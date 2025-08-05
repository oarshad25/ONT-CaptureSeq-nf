/*
Convert annotation gtf to bed12 using minimap2 script paftools

Minimap2 can optionally take input annotation and prioritise on annotated splice junctions.
To use this feature we need to convert the annotation GTF to BED12.
This module uses Minimap's paftools scripts gff2bed utility to do the conversion.
*/

process MINIMAP2_PAFTOOLS_GFF2BED {
    label 'low'

    conda "bioconda::minimap2=2.28"
    container "${workflow.containerEngine == 'apptainer'
        ? 'https://depot.galaxyproject.org/singularity/minimap2:2.28--h577a1d6_4'
        : 'quay.io/biocontainers/minimap2:2.28--h577a1d6_4'}"

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
