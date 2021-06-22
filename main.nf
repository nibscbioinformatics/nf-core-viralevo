#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core-viralevo
========================================================================================
     nf-core-viralevo Analysis Pipeline
     #### Homepage / Documentation
     https://github.com/kaurravneet4123
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

log.info Utils.logo(workflow, params.monochrome_logs)

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/conva --input samplesheet.csv -profile conda"
    log.info NfcoreSchema.paramsHelp(workflow, params, json_schema, command)
    exit 0
}

========================================================================================
/* --        FASTA PARAMETER VALUE                                                 -- */
========================================================================================


params.genome = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE AND PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params, json_schema)
log.info NfcoreSchema.paramsSummaryLog(workflow, params, json_schema)
log.info Workflow.citation(workflow)
log.info Utils.dashedLine(params.monochrome_logs)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

workflow NFCORE_VIRALEVO {

    //
    // WORKFLOW: Run main nf-core-viralevo analysis pipeline
    //
    include { VIRALEVO } from './workflows/viralevo' addParams( summary_params: summary_params )
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
