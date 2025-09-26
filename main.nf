#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RFDiffusionWorkflow } from './workflows/rfdiffusion.nf'
include { FilterRFD ; RunRFDiffusion } from './modules/rfdiffusion.nf'
include { PrepFAMPNN ; FilterFAMPNN; RunFAMPNN } from './modules/fampnn.nf'
include { FilterMPNN; PrepMPNN ; RunMPNN } from './modules/proteinmpnn.nf'
include { AlignAF2; FilterAF2; RunAF2 } from './modules/af2.nf'
include { AnalyseBestDesigns } from './modules/analysis.nf'
include { PublishResults } from './modules/publish.nf'
include { AlignBoltz ; FilterBoltz; PrepBoltz ; RunBoltz } from './modules/boltz.nf'
include { CombineMetadata } from './modules/combine_metadata.nf'
include { Compress as CompressRFD } from './modules/compress'
include { Compress as CompressMPNN } from './modules/compress'
include { Compress as CompressFAMPNN } from './modules/compress'
include { Compress as CompressAF2 } from './modules/compress'
include { Compress as CompressBoltz } from './modules/compress'
include { MergeUncroppedTarget } from './modules/merge_uncropped_target.nf'

workflow {
    // Permit use of topic channels in Nextflow v24 by enabling preview features
    try {
        nextflow.preview.topic = true
    } catch (Exception _e) {
        // Silently continue - topic channels work by default in v25+
    }
    
    def outputDirectory = params.out_dir

    if (params.run_rfd_only && (params.skip_rfd_seq || params.skip_rfd_seq_pred )) {
        error("Cannot use --run_rfd_only with skip flags --skip_rfd_seq or --skip_rfd_seq_pred. These options are contradictory.")
    }
    if (params.run_rfd_only && params.skip_rfd) {
        error("Cannot use --run_rfd_only with --skip_rfd. These options are contradictory.")
    }

    // Calculate batch size based on maximum GPUs
    def num_batches = Math.min(params.gpus, params.rfd_num_designs).intValue()
    def batch_size = Math.ceil(params.rfd_num_designs / num_batches).intValue()
    def num_designs = num_batches * batch_size
    
    println("***********************************************************************")
    println("██████╗ ██████╗  ██████╗ ████████╗███████╗██╗███╗   ██╗██████╗      ██╗")
    println("██╔══██╗██╔══██╗██╔═══██╗╚══██╔══╝██╔════╝██║████╗  ██║██╔══██╗     ██║")
    println("██████╔╝██████╔╝██║   ██║   ██║   █████╗  ██║██╔██╗ ██║██║  ██║     ██║")
    println("██╔═══╝ ██╔══██╗██║   ██║   ██║   ██╔══╝  ██║██║╚██╗██║██║  ██║██   ██║")
    println("██║     ██║  ██║╚██████╔╝   ██║   ███████╗██║██║ ╚████║██████╔╝╚█████╔╝")
    println("╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚══════╝╚═╝╚═╝  ╚═══╝╚═════╝  ╚════╝ ")
    println("                   ProteinDJ Protein Design Pipeline                   ")
    println("          Developers: Dylan Silke, Josh Hardy, Julie Iskander          ")
    println("***********************************************************************")
    println("* Pipeline Mode: ${params.rfd_mode}")
    println("* Number of RFdiffusion designs: ${num_designs}")
    println("* Number of sequences for each design: ${params.seqs_per_design}")
    println("* Output Directory: ${outputDirectory}")
    println("***********************************************************************\n")

    // Create output directory for copy of config files used in run
    def configDir = file("${outputDirectory}/configs")
    configDir.mkdirs()
    workflow.configFiles.each { configFile ->
        configFile.copyTo("${configDir}/${configFile.getName()}")
    }

    // Create output directory for copy of input files used in run
    def inputsDir = file("${outputDirectory}/inputs")
    inputsDir.mkdirs()

    ///////////////////////
    // FOLD DESIGN STAGE //
    ///////////////////////

    // Run RFdiffusion if not skipped
    if (!params.skip_rfd & !params.skip_rfd_seq & !params.skip_rfd_seq_pred) {
        // Check if num_designs has been provided
        if (!params.rfd_num_designs) {
            error("Please provide the number of designs for RFdiffusion to generate")
        }
        def rfdParams = new RFDiffusionParams(params)
        // Generate the command string
        def rfdCommand = rfdParams.generateCommandString()
        log.info("RFdiffusion command: ${rfdCommand} inference.num_designs=${batch_size}")

        // Collect input files
        def inputFiles = collectInputFiles(params)

        // Copy input files to output directory
        inputFiles.each { inputFile ->
            "rsync -r ${inputFile} ${inputsDir}/.".execute()
        }

        // Launch RFDiffusion Workflow
        RFDiffusionWorkflow(
            rfdCommand,
            params.rfd_num_designs,
            batch_size,
            params.rfd_mode,
            inputFiles,
        )

        RFDiffusionWorkflow.out.pdbs_jsons.set { rfd_pdbs_jsons }
        // Compress output files
        CompressRFD("rfd", rfd_pdbs_jsons.flatten().collect())

        // Batch RFD PDBs and JSONS for CPU tasks
        Utils
            .rebatchTuples(rfd_pdbs_jsons, 200)
            .set { rfd_tuples }
        // RFdiffusion filtering - secondary structure and radius of gyration
        FilterRFD(rfd_tuples)

        // If Running RFD only these are the final pdbs
        if (params.run_rfd_only) {
            FilterRFD.out.pdbs_jsons
                .flatten()
                .collect()
                .ifEmpty(file("${projectDir}/lib/placeholder.pdb"))
                .set { final_pdbs }
        }
        else {
            FilterRFD.out.pdbs_jsons.set { filt_rfd_pdbs_jsons }
        }
    }
    else if (params.skip_rfd & !params.skip_rfd_seq & !params.skip_rfd_seq_pred) {
        // Skip RFDiffusion and use existing PDBs and JSONs from specified directory
        println("Skipping RFDiffusion stage as skip_rfd=true.")
        println("Running Sequence Design, Prediction, and Analysis stages only.")
        println("Looking for PDBs and JSONs in: ${params.skip_input_dir}")
        // Check if directory exists
        if (!file(params.skip_input_dir).exists()) {
            throw new FileNotFoundException("Skip input file directory not found at: ${params.skip_input_dir}. Please ensure the path is correct.")
        }
        def previous_pdbs = file("${params.skip_input_dir}").listFiles().findAll { it.name.endsWith('.pdb') }
        def previous_jsons = file("${params.skip_input_dir}").listFiles().findAll { it.name.endsWith('.json') }
        // Error handling for missing files
        if (previous_pdbs.isEmpty()) {
            throw new FileNotFoundException("No PDB files found in directory: ${params.skip_input_dir}. Please provide PDB files to proceed with the workflow.")
        }
        if (previous_jsons.isEmpty()) {
            throw new FileNotFoundException("No JSON files found in directory: ${params.skip_input_dir}. Please provide JSON files to proceed with the workflow.")
        }
        println("Found ${previous_pdbs.size()} PDB files")
        println("Found ${previous_jsons.size()} JSON files\n")

        // Copy PDB and JSON files from the previous results directory to inputs directory
        previous_pdbs.each { pdbFile ->
            pdbFile.copyTo("${inputsDir}/${pdbFile.getName()}")
        }
        previous_jsons.each { jsonFile ->
            jsonFile.copyTo("${inputsDir}/${jsonFile.getName()}")
        }

        // Create channel with PDB-JSON tuples from specified directory
        Channel
            .of([previous_pdbs, previous_jsons])
            .set { rfd_pdbs_jsons }
        // Batch RFD PDBs and JSONS for CPU tasks
        Utils
            .rebatchTuples(rfd_pdbs_jsons, 200)
            .set { filt_rfd_pdbs_jsons }
    }
    else {
        println("Skipping RFDiffusion stage as skip_rfd_seq=true or skip_rfd_seq_pred=true.")
    }
    ///////////////////////////
    // SEQUENCE DESIGN STAGE //
    ///////////////////////////
    // Run Sequence Design if not skipped
    if (!params.skip_rfd_seq & !params.skip_rfd_seq_pred & !params.run_rfd_only) {
        // Sequence design (either MPNN or FAMPNN)
        if (params.seq_method == "mpnn") {
            // Add FIXED labels to PDBs for target residues so the sequence does not change
            PrepMPNN(filt_rfd_pdbs_jsons)
            
            // Method specific batching
            if (params.mpnn_relax_max_cycles > 0) {
                // use smaller batches for fast relax (slow)
                PrepMPNN.out.pdbs
                    .collect()
                    .flatten()
                    .buffer( size: 2, remainder: true )
                    .set { seq_input_pdbs }
            }
            else {
                // use larger batches without fast relax
                PrepMPNN.out.pdbs
                    .collect()
                    .flatten()
                    .buffer( size: 10, remainder: true )
                    .set { seq_input_pdbs }
            }

            // Launch ProteinMPNN
            RunMPNN(seq_input_pdbs)

            // Compress output files
            CompressMPNN("mpnn", RunMPNN.out.pdbs_jsons.flatten().collect())

            // Rebatch sequence assignment files for CPU Filtering Step
            Utils
                .rebatchTuples(RunMPNN.out.pdbs_jsons, 200)
                .set { seq_tuple }

            // Filter designs by sequence score
            FilterMPNN(seq_tuple)
            FilterMPNN.out.pdbs
                .flatten()
                .collect()
                .set { filt_seq_pdbs }
        }
        else if (params.seq_method == "fampnn") {
            // FAMPNN path
            // Rebatch files for Prep Step
            Utils
                .rebatchTuples(filt_rfd_pdbs_jsons, 10)
                .set { fampnn_prep_input_tuple }
            
            // Restore side-chains to RFD output and prepare CSV file with fixed residues
            PrepFAMPNN(fampnn_prep_input_tuple)
            PrepFAMPNN.out.csv
                .collectFile(name: 'merged_results.csv', keepHeader: true)
                .set { mega_csv }

            // GPU-aware batching for RunFAMPNN
            Utils
                .rebatchGPU(PrepFAMPNN.out.pdbs, params.gpus)
                .set { fampnn_pdbs }

            // Add CSV path to PDB channel
            fampnn_pdbs
                .combine(mega_csv)
                .set { fampnn_input }

            if (params.rfd_mode in ['binder_denovo', 'binder_foldconditioning', 'binder_motifscaffolding', 'binder_partialdiffusion']) {
                // Perform design and scoring on binder (chain A)
                RunFAMPNN(fampnn_input, 'A')
            }
            else {
                // Perform design and scoring on all chains
                RunFAMPNN(fampnn_input, 'all_chains')
            }

            // Compress output files
            CompressFAMPNN("fampnn", RunFAMPNN.out.pdbs_jsons.flatten().collect())

            // Rebatch sequence assignment files for CPU Filtering Step
            Utils
                .rebatchTuples(RunFAMPNN.out.pdbs_jsons, 200)
                .set { seq_tuple }

            // Filter designs by sequence score
            FilterFAMPNN(seq_tuple)
            FilterFAMPNN.out.pdbs
                .flatten()
                .collect()
                .set { filt_seq_pdbs }
        }
        else {
            error("Not a valid sequence assignment method")
        }
    }
    else if (!params.skip_rfd_seq_pred & !params.run_rfd_only) {
        // Skip sequence design and run prediction using existing PDBs from specified directory
        println("Skipping Sequence Design stage as skip_rfd_seq=true.")
        println("Running Prediction and Analysis stages only.")
        println("Looking for PDBs in: ${params.skip_input_dir}")
        // Check if directory exists
        if (!file(params.skip_input_dir).exists()) {
            throw new FileNotFoundException("Skip input file directory not found at: ${params.skip_input_dir}. Please ensure the path is correct.")
        }
        def pdbs_for_pred = file("${params.skip_input_dir}").listFiles().findAll { it.name.endsWith('.pdb') }
        if (pdbs_for_pred.isEmpty()) {
            throw new FileNotFoundException("No PDB files found in directory: ${params.skip_input_dir}. Please provide PDB files to proceed with the workflow.")
        }
        println("Found ${pdbs_for_pred.size()} PDB files")

        // Copy PDB files from the previous results directory to inputs directory
        pdbs_for_pred.each { pdbFile ->
            pdbFile.copyTo("${inputsDir}/${pdbFile.getName()}")
        }

        // Create channel with PDBs from specified directory
        Channel
            .of(pdbs_for_pred)
            .set { filt_seq_pdbs }
    }
    else if (params.skip_rfd_seq_pred) {
        println("Skipping Sequence Design stage as skip_rfd_seq_pred=true.")
    }
    else {
        println("Skipping Sequence Design stage as run_rfd_only=true.")
    }
    ////////////////////////////////
    // STRUCTURE PREDICTION STAGE //
    ////////////////////////////////
    // Run Structure Prediction if not skipped
    if (!params.skip_rfd_seq_pred & !params.run_rfd_only) {
        // Optional uncropped target PDB merge for binder design
        if (params.rfd_mode in ['binder_denovo', 'binder_foldconditioning', 'binder_motifscaffolding', 'binder_partialdiffusion']) {
            // if uncropped target PDB file is provided, merge with designs
            if (params.uncropped_target_pdb) {
                def uncroppedPDBfile = file(params.uncropped_target_pdb)
                if (!uncroppedPDBfile.exists()) {
                    throw new FileNotFoundException("Uncropped target PDB file not found at path: ${params.uncropped_target_pdb}. Please ensure the file exists and the path is correct.")
                }
                MergeUncroppedTarget(filt_seq_pdbs, uncroppedPDBfile).set { pred_input_pdbs }
            }
            else {
                filt_seq_pdbs.set { pred_input_pdbs }
            }
        } else {
            filt_seq_pdbs.set { pred_input_pdbs }
        }
        // Structure Prediction (either AlphaFold2 Initial-Guess or Boltz-2)
        if (params.pred_method == "af2") {
          
            // reallocate batching for GPU
            Utils
                .rebatchGPUByNumRes(pred_input_pdbs, params.gpus)
                .set { pred_input_tuple }
            
            // AlphaFold2-Initial Guess
            RunAF2(pred_input_tuple)

            // Compress output files
            CompressAF2("af2", RunAF2.out.pdbs_jsons.flatten().collect())

            // Batch files for CPUs
            Utils
                .rebatchTuples(RunAF2.out.pdbs_jsons, 200)
                .set { af2_tuple }

            // Filtering of AF2 results
            FilterAF2(af2_tuple)

            if (params.rfd_mode in ['binder_denovo', 'binder_foldconditioning', 'binder_motifscaffolding', 'binder_partialdiffusion']) {
                // Alignment of PDBs to target chain(s). Only need one reference file
                AlignAF2(FilterAF2.out.pdbs.flatten().collect(), pred_input_pdbs.flatten().last())
                AlignAF2.out.pdbs
                    .flatten()
                    .collect()
                    .set { analysis_input_pdbs }
            } else {
                FilterAF2.out.pdbs
                    .flatten()
                    .collect()
                    .set { analysis_input_pdbs }
            }
        }
        else if (params.pred_method == "boltz") {
            // Prep yaml files for Boltz-2
            PrepBoltz(pred_input_pdbs)

            // reallocate batching for GPU
            Utils
                .rebatchGPU(PrepBoltz.out.yamls, params.gpus)
                .set { pred_input_tuple }
            
            // Perform prediction of designs using Boltz-2
            RunBoltz(pred_input_tuple)

            // Batch files for CPUs
            Utils
                .rebatchTuples(RunBoltz.out.pdbs_jsons, 200)
                .set { boltz_tuple }

            // Align Boltz Predictions to FAMPNN output and calculate RMSD
            if (params.rfd_mode in ['binder_denovo', 'binder_foldconditioning', 'binder_motifscaffolding', 'binder_partialdiffusion']) {
                AlignBoltz(boltz_tuple, filt_seq_pdbs, 'binder')
            }
            else {
                AlignBoltz(boltz_tuple, filt_seq_pdbs, 'monomer')
            }
            // Compress output files
            CompressBoltz("boltz", AlignBoltz.out.pdbs_jsons.flatten().collect())

            // Filtering of Boltz-2 results
            FilterBoltz(AlignBoltz.out.pdbs_jsons)
            FilterBoltz.out.pdbs
                .flatten()
                .collect()
                .set { analysis_input_pdbs }
        }
        else {
            error("Not a valid structure prediction method")
        }
    }
    else if (!params.run_rfd_only) {
        // Skip prediction and run analysis only using existing PDBs from specified directory
        println("Skipping Structure Prediction stage as skip_rfd_seq_pred=true.")
        println("Running Analysis Stage only")
        println("Looking for PDBs in: ${params.skip_input_dir}")
        if (!file(params.skip_input_dir).exists()) {
            throw new FileNotFoundException("Skip input file directory not found at: ${params.skip_input_dir}. Please ensure the path is correct.")
        }
        def pdbs_for_analysis = file("${params.skip_input_dir}").listFiles().findAll { it.name.endsWith('.pdb') }
        if (pdbs_for_analysis.isEmpty()) {
            throw new FileNotFoundException("No PDB files found in directory: ${params.skip_input_dir}. Please provide PDB files to proceed with the workflow.")
        }
        println("Found ${pdbs_for_analysis.size()} PDB files")

        // Copy PDB files from the previous results directory to inputs directory
        pdbs_for_analysis.each { pdbFile ->
            pdbFile.copyTo("${inputsDir}/${pdbFile.getName()}")
        }

        // Create channel with PDBs from specified directory
        Channel
            .of(pdbs_for_analysis)
            .set { analysis_input_pdbs }
    }
    else {
        println("Skipping Structure Prediction stage as run_rfd_only=true.")
    }
    ////////////////////
    // ANALYSIS STAGE //
    ////////////////////
    if (!params.run_rfd_only) {
        // Analysis of PDBs to generate additional metrics 
        AnalyseBestDesigns(analysis_input_pdbs)
        // Use placeholder PDB file if no designs survive filtering
        analysis_input_pdbs
            .flatten()
            .collect()
            .ifEmpty(file("${projectDir}/lib/placeholder.pdb"))
            .set { final_pdbs }
    }
    else {
        println("Skipping Analysis stage as run_rfd_only=true.")
    }

    // Open topic channels to collect metadata for all designs. 
    // Channel for metadata with only fold_id and not seq_id
    channel
        .topic('metadata_ch_fold')
        .flatten()
        .collectFile(name: "metadata_fold.jsonl", newLine: true)
        .ifEmpty { file("${projectDir}/lib/empty-meta.jsonl") }
        .set { metadata_fold }
    // Channel for metadata with both fold_id and seq_id
    channel
        .topic('metadata_ch_fold_seq')
        .flatten()
        .collectFile(name: "metadata_fold_seq.jsonl", newLine: true)
        .ifEmpty { file("${projectDir}/lib/empty-meta.jsonl") }
        .set { metadata_fold_seq }

    // Combine Metadata into CSV
    CombineMetadata(metadata_fold, metadata_fold_seq).csv.collectFile(name: "all_designs.csv").set { all_designs_metadata }

    // Count outputs
    if (params.run_rfd_only) {
        Utils.countPdbFiles(rfd_tuples).set { rfd_count }
        Utils.countPdbFiles(final_pdbs).set { filter_rfd_count }
        seq_count = 0
        filter_seq_count = 0
        filter_pred_count = 0
    }
    else if (params.skip_rfd_seq_pred) {
        rfd_count = 0
        filter_rfd_count = 0
        seq_count = 0
        filter_seq_count = 0
        filter_pred_count = 0
    }
    else if (params.skip_rfd_seq) {
        rfd_count = 0
        filter_rfd_count = 0
        seq_count = 0
        filter_seq_count = 0
        Utils.countPdbFiles(analysis_input_pdbs).set { filter_pred_count }
    }
    else if (params.skip_rfd) {
        rfd_count = 0
        filter_rfd_count = 0
        Utils.countPdbFiles(seq_tuple).set { seq_count }
        Utils.countPdbFiles(filt_seq_pdbs).set { filter_seq_count }
        Utils.countPdbFiles(analysis_input_pdbs).set { filter_pred_count }
    }
    else {
        Utils.countPdbFiles(rfd_tuples).set { rfd_count }
        Utils.countPdbFiles(filt_rfd_pdbs_jsons).set { filter_rfd_count }
        Utils.countPdbFiles(seq_tuple).set { seq_count }
        Utils.countPdbFiles(filt_seq_pdbs).set { filter_seq_count }
        Utils.countPdbFiles(analysis_input_pdbs).set { filter_pred_count }
    }

    // Generate report and statistics of run
    PublishResults(
        final_pdbs,
        all_designs_metadata,
        rfd_count,
        filter_rfd_count,
        seq_count,
        filter_seq_count,
        filter_pred_count
    )
    
    // Save log file on completion
    workflow.onComplete {
        def logFile = file('.nextflow.log')
        def outputDir = file(outputDirectory)
        if (logFile.exists()) {
            logFile.copyTo(outputDir.resolve('nextflow.log'))
        }
    }
}

// Collect required input files for RFdiffusion
def collectInputFiles(params) {
    def inputs = []

    // Add required input files
    if (params.rfd_mode in [
        'binder_denovo',
        'binder_foldcond',
        'binder_motifscaff',
        'binder_partialdiff',
        'monomer_motifscaff',
        'monomer_partialdiff',
    ]) {
        if (params.rfd_input_pdb) {
            inputs << file(params.rfd_input_pdb)
        }
    }
    if (params.rfd_mode in ['monomer_denovo', 'monomer_foldcond']) {
        // Add 'placeholder' PDB file, since RFdiffusion requires xyz coordinates
        inputs << file("${projectDir}/lib/placeholder.pdb")
    }
    if (params.rfd_mode in ['binder_foldcond', 'monomer_foldcond']) {
        if (params.rfd_scaffold_dir) {
            // Add scaffolds_dir and contents
            inputs << file(params.rfd_scaffold_dir)
        }
    }
    if (params.rfd_mode == 'binder_foldcond') {
        if (params.rfd_target_ss) {
            inputs << file(params.rfd_target_ss)
        }
        if (params.rfd_target_adj) {
            inputs << file(params.rfd_target_adj)
        }
    }

    return inputs
}
