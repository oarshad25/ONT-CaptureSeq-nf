/*
* Align reads to the reference with minimap2 and directly stream SAM output
* into samtools sort to produce sorted BAM, BAM index and flagstat metrics.
*/

process MINIMAP2_SAMTOOLS {
    tag "${meta.id}"

    label 'medium'

    publishDir "${params.outdir}/bam/unfiltered/", mode: "link"

    conda "bioconda::minimap2=2.30 bioconda::samtools=1.23"

    container "${workflow.containerEngine == 'apptainer'
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools:a81ff6397062f3f9'
        : 'community.wave.seqera.io/library/minimap2_samtools:f2e55a1cd407fcaf'}"

    input:
    // input reads to align
    tuple val(meta), path(fastq)

    // genome fasta file
    path genome

    // junction bed file, optional
    path junc_bed

    // preset, -x argument to Minimap2, e.g. "splice"
    val x

    // kmer, -k argument
    val k

    // -u argument
    val u

    // -G argument
    val G

    // -I argument
    val I

    //--cs argument
    val cs

    //--junc-bonus argument
    val junc_bonus

    // any additional options
    val extra_opts

    output:
    // sorted and indexed bam
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bambai

    // flagstats
    tuple val(meta), path("${meta.id}.flagstat.txt"), emit: flagstat

    script:
    def preset_flag = x ? "-x ${x}" : "-x splice"
    def kmer_flag = k ? "-k ${k}" : ""
    def strand_flag = u ? "-u ${u}" : ""
    def max_intron_length_flag = G ? "-G ${G}" : ""
    def indexing_batchsize_flag = I ? "-I ${I}" : ""
    def cs_flag = cs ? "--cs=${cs}" : ""
    def junc_bed_flag = junc_bed.name != 'NO_FILE' ? "--junc-bed ${junc_bed}" : ''
    def junc_bonus_flag = junc_bonus ? "--junc-bonus ${junc_bonus}" : ""

    """
    minimap2 \\
        -t ${task.cpus} \\
        -a \\
        ${preset_flag} \\
        ${kmer_flag} \\
        ${strand_flag} \\
        ${max_intron_length_flag} \\
        ${indexing_batchsize_flag} \\
        ${cs_flag} \\
        --secondary=yes \\
        ${extra_opts} \\
        ${junc_bed_flag} \\
        ${junc_bonus_flag} \\
        ${genome} \\
        ${fastq} | \\
    # convert SAM to BAM
    samtools sort -@ ${task.cpus} -O bam -o ${meta.id}.bam -

    # index BAM
    samtools index -@ ${task.cpus} ${meta.id}.bam

    # post-mapping QC
    samtools flagstat -@ ${task.cpus} ${meta.id}.bam > ${meta.id}.flagstat.txt
    """
}
