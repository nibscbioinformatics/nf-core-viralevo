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

include { VIRALEVO } from './workflows/viralevo'
    
//
// WORKFLOW: Run main nf-core-viralevo analysis pipeline
//
workflow NFCORE_VIRALEVO {
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
