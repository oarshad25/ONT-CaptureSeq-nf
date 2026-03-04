#!/usr/bin/env nextflow

/*
================================================================================
   Utility functions
================================================================================
*/

/*
* Helper functions
*/


// Function to find all files with extensions in 'extns' in directory 'inputdir'.
// Returns list of all files found.
// If no matching files found and flag 'check' is set, abort with error

def findFilesWithExtns(inputdir, extns) {
    // whether to error out if no files matching extension found in inputdir
    def check = true

    // https://github.com/asadprodhan/How-to-channel-sequencing-reads-from-multiple-subdirectories-into-nextflow-pipeline
    def all_files = file(inputdir).listFiles()

    // get files with extensions in fastq_extns
    def files_with_extns = all_files.findAll { fn -> extns.find { fn.name.endsWith(it) } }

    if (!files_with_extns) {
        if (check) {
            error("No files with extensions ${extns} found in ${inputdir}")
            System.exit(1)
        }
        else {
            log.warn("No files with extensions ${extns} found in ${inputdir}")
            return files_with_extns
        }
    }
    else {
        return files_with_extns
    }
}

// find fastq files in fastqdir
def findFastqFiles(fastqdir) {
    def fastq_extns = ['.fq', '.fq.gz', '.fastq', '.fastq.gz']

    findFilesWithExtns(fastqdir, fastq_extns)
}

/**
 * Resolves the featurecounts strandedness parameter (-s flag).
 *
 * If `params.featurecounts_s` is explicitly provided, it is validated and used directly.
 * If not provided (null or empty), the value is inferred from `params.run_cdna_qc`:
 *   - true  → 1 (stranded, for cDNA QC workflows)
 *   - false → 0 (unstranded)
 *
 * @param params  The Nextflow params object.
 * @return        An integer: 0 (unstranded), 1 (stranded), or 2 (reverse stranded).
 * @throws        Exception if `params.featurecounts_s` is provided but not a valid value (0, 1, or 2).
 */
def resolveFeatureCountsStrand(params) {
    if (params.featurecounts_s != null && params.featurecounts_s != "") {
        // normalise to inetger if not already
        def value = params.featurecounts_s instanceof Integer
            ? params.featurecounts_s
            : params.featurecounts_s.toInteger()

        // validate provided parameter value if provided and error out if invalid
        if (!(value in [0, 1, 2])) {
            error(
                "Invalid value for params.featurecounts_s: '${params.featurecounts_s}'. " + "Must be 0 (unstranded), 1 (stranded), or 2 (reverse stranded)."
            )
        }
        return value
    }

    // if `params.featurecounts_s` not provided, set starnd based on `params.run_cdna_qc`
    return params.run_cdna_qc ? 1 : 0
}
