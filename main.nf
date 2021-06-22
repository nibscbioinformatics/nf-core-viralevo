#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core-viralevo
========================================================================================
     nf-core-viralevo Analysis Pipeline
     #### Homepage / Documentation
     Github: https://github.com/nibscbioinformatics/nf-core-viralevo
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
           GENOME PARAMETER VALUES                                                 
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gff   = WorkflowMain.getGenomeAttribute(params, 'gff')

/*
========================================================================================
    VALIDATE AND PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow NFCORE_VIRALEVO {

    //
    // WORKFLOW: Run main nf-core-viralevo analysis pipeline
    //
    include { VIRALEVO } from './workflows/viralevo'
    VIRALEVO ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_VIRALEVO ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
