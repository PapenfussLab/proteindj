#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RunRFDiffusion } from '../modules/rfdiffusion.nf'

workflow RFDiffusionWorkflow {
    take:
    rfdCommand // The prepared RFDiffusion command
    numDesigns // Number of designs to generate
    batchSize // Size of each batch
    mode // Mode for RFDiffusion
    inputFiles // Input files for RFDiffusion

    main:
    // Create the channel for RFDiffusion
    rf_ch = Channel
        .fromList((0..<numDesigns).collate(batchSize))
        .map { batch ->
            def batchId = batch.isEmpty() ? 0 : (batch[0] / batchSize).intValue()
            def designStartnum = batch.min()
            tuple(
                rfdCommand,
                batchId,
                batchSize,
                designStartnum,
                mode,
                inputFiles,
            )
        }

    // Run RFDiffusion with the generated channel
    RunRFDiffusion(rf_ch)

    emit:
    pdbs = RunRFDiffusion.out.pdbs
    pdbs_jsons = RunRFDiffusion.out.pdbs_jsons
}
