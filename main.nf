params.fastqs
params.fasta
params.outdir
// bed intervals of interest
// params.regions

outdir = file(params.outdir)
// regions = file(params.regions)

if( !params.fasta ) { exit 1, "--fasta is not defined" }
// assumes .alt, .amb, .ann, .bwt, .pac, and .sa are present
fasta = file(params.fasta)

Channel
    .value(file("${params.fasta}.{amb,ann,bwt,pac,sa}"))
    .set { bwaidx_ch }

// emits [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
Channel
    .fromFilePairs(params.fastqs, flat: true)
    .set { fastq_ch }


// qc the fastqs
process fastp {
    publishDir path: "$outdir/json", pattern: "*.json"
    publishDir path: "$outdir/html", pattern: "*.html"

    input:
    set sample_id, file(r1), file(r2) from fastq_ch

    output:
    set sample_id, file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") into checked_fastq_ch
    file("${sample_id}.fastp.json") into fastp_report_ch
    file("${sample_id}.fastp.html")

    script:
    """
    fastp --thread ${task.cpus} --in1 $r1 --out1 ${sample_id}_R1.fastq.gz \
        --in2 $r2 --out2 ${sample_id}_R2.fastq.gz -y \
        --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html
    """
}


process bwamem {
    input:
    set sample_id, file(r1), file(r2) from checked_fastq_ch
    file(fasta)
    file(bwaidx) from bwaidx_ch

    output:
    set sample_id, file("${sample_id}.bam") into bwa_ch

    script:
    rg = "@RG\\tID:${sample_id}\\tPU:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:illumina"
    """
    bwa mem -K 100000000 -R \"${rg}\" -t ${task.cpus} -M ${fasta} $r1 $r2 \
        | samtools sort -n --threads ${task.cpus} -m 2G --output-fmt BAM -o ${sample_id}.bam
    """
}


process markduplicates {
    tag "$sample_id"
    publishDir path: "$outdir/alignments", overwrite: true

    input:
    set sample_id, file(bam) from bwa_ch

    output:
    set sample_id, file("${sample_id}.md.bam"), file("${sample_id}.md.bam.bai") into md_ch

    script:
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g MarkDuplicatesSpark \
        --input $bam \
        --output ${sample_id}.md.bam \
        --tmp-dir . \
        --spark-master \'local[*]\'
    gatk --java-options -Xmx${task.memory.toGiga()}g BuildBamIndex \
        --INPUT ${sample_id}.md.bam \
        --OUTPUT ${sample_id}.md.bam.bai \
        --TMP_DIR .
    """
}


process alignstats {
    publishDir path: "$outdir/json", pattern: "*.json"

    input:
    set sample_id, file(bam), file(bai) from md_ch

    output:
    file("*.json")

    script:
    """
    alignstats -C -P ${task.cpus} -i $bam -j bam -o ${sample_id}.alignstats.json
    """
}
