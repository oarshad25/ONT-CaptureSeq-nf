/*
Convert annotation gtf to bed12 using minimap2 script paftools

Minimap2 can optionally take input annotation and prioritise on annotated splice junctions.
To use this feature we need to convert the annotation GTF to BED12.
This module uses Minimap's paftools scripts gff2bed utility to do the conversion.
*/

process MINIMAP2_PAFTOOLS_GFF2BED {
    label 'low'

    conda "bioconda::minimap2=2.30 bioconda::samtools=1.23"
    container "${workflow.containerEngine == 'apptainer'
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools:a81ff6397062f3f9'
        : 'community.wave.seqera.io/library/minimap2_samtools:f2e55a1cd407fcaf'}"

    input:
    // path to annotation gtf
    path annotation

    output:
    // annotation in bed12 format
    path "*.bed", emit: bed12

    script:
    """
    paftools.js gff2bed ${annotation} > ${annotation.baseName}.bed
    """
}
