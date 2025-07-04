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
