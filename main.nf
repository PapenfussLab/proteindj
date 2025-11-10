#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GenerateRFDContigs; GenerateRFDfoldcond; FilterFold ; RunRFD } from './modules/rfdiffusion.nf'
include { AnalyseBC; PrepBC; RunBC } from './modules/bindcraft.nf'
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
include { Compress as CompressBC } from './modules/compress'
include { MergeUncroppedTarget } from './modules/merge_uncropped_target.nf'

workflow {
    // Permit use of topic channels in Nextflow v24 by enabling preview features
    try {
        nextflow.preview.topic = true
    } catch (Exception _e) {
        // Silently continue - topic channels work by default in v25+
    }
    
    def outputDirectory = params.out_dir

    if (params.run_fold_only && (params.skip_fold_seq || params.skip_fold_seq_pred )) {
        error("Cannot use --run_fold_only with skip flags --skip_fold_seq or --skip_fold_seq_pred. These options are contradictory.")
    }
    if (params.run_fold_only && params.skip_fold) {
        error("Cannot use --run_fold_only with --skip_fold. These options are contradictory.")
    }

    // Calculate batch size based on maximum GPUs
    def num_batches = Math.min(params.gpus, params.num_designs).intValue()
    def batch_size = Math.ceil(params.num_designs / num_batches).intValue()
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
    println("* Pipeline Mode: ${params.design_mode}")
    println("* Number of designs: ${num_designs}")
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

    // Run Fold Design if not skipped
    if (!params.skip_fold & !params.skip_fold_seq & !params.skip_fold_seq_pred) {
        // Check if num_designs has been provided
        if (!params.num_designs) {
            error("Please provide the number of designs to generate")
        }
        // Validate input PDB file
        if (params.design_mode in ['bindcraft','binder_denovo', 'binder_motifscaff', 'binder_partialdiff', 'binder_foldcond', 'monomer_motifscaff', 'monomer_partialdiff']) {
            if (!params.input_pdb) {
                throw new IllegalArgumentException("Please provide input PDB file path required by $params.design_mode mode")
            }
            inputFile = new File(params.input_pdb)
            if (!inputFile.exists()) {
                throw new FileNotFoundException("Input PDB file not found at path: ${params.input_pdb}. Please ensure the file exists and the path is correct.")
            }
        }
        // Validate design length
        if (params.design_mode in ['bindcraft','binder_denovo', 'monomer_denovo']){
            validateDesignLength(params.design_length)
        }
        
        
        if (params.design_mode=="bindcraft"){
            // Use Bindcraft for fold design
            validateBindCraftParams(
                params.bc_chains,
                params.hotspot_residues,
                params.design_length,
                params.num_designs,
                params.input_pdb,
                params.bc_advanced_json)
            log.info("Using BindCraft to hallucinate binders with the following design parameters:")
            log.info("* Design length = ${params.design_length}")
            if (params.bc_chains){
                log.info("* Target chains = ${params.bc_chains}")
            }
            if (params.hotspot_residues){
                log.info("* Target hotspots = ${params.hotspot_residues}")
            }
            // Select advanced settings json file    
            def bc_advanced_json 
            if(!params.bc_advanced_json){
                bc_advanced_json = getAdvancedSettingsPath(
                    params.bc_design_protocol,
                    params.bc_template_protocol
                )
            } else {
                bc_advanced_json = file(params.bc_advanced_json)
            }

            // Get path of filter settings json file
            def bc_filters_json = file("${projectDir}/lib/bindcraft/settings_filters/no_filters.json")

            // Collect input files
            def inputFiles = collectInputFiles(params)

            // Copy input files to output directory
            inputFiles.each { inputFile ->
                "rsync -r ${inputFile} ${inputsDir}/.".execute()
            }
            
            // Create channel with items for requested designs
            bc_ch = Channel
                .fromList((0..<num_batches))

            PrepBC(bc_ch,batch_size,file(params.input_pdb),bc_advanced_json)

            // Run BindCraft for each batch
            RunBC(PrepBC.out, bc_filters_json, file(params.input_pdb))

            // Collect batches and run analysis and conversion of BindCraft outputs
            AnalyseBC(RunBC.out.pdbs_csvs.flatten().collect()) 
            AnalyseBC.out.pdbs_jsons.set { bc_pdbs_jsons }
            
            // Compress output files
            CompressBC("bc", bc_pdbs_jsons.flatten().collect())

            // Batch PDBs and JSONS for CPU tasks
            Utils
                .rebatchTuples(bc_pdbs_jsons, 200)
                .set { fold_tuples }

        } else { // Use RFdiffusion for fold design
            validateRFDParameters(params)
            // Check for user-provided contigs or whether to automatically generate them
            if (params.design_mode in ['binder_foldcond', 'monomer_foldcond']){
               Channel.value('NoContigsNeededForFoldConditioning').set{rfdContigs}
            } else if (params.rfd_contigs && params.design_mode in ['binder_denovo', 'binder_partialdiff', 'binder_motifscaff', 'monomer_denovo','monomer_motifscaff', 'monomer_partialdiff']){
               // Use provided value
               description=generateContigDescription("$params.rfd_contigs")
               println description
               Channel.value(params.rfd_contigs).set{rfdContigs}
            } else if (!params.rfd_contigs && params.design_mode in ['binder_motifscaff', 'monomer_motifscaff']){
                error("rfd_contigs is required for $params.design_mode mode.")
            } else if (!params.rfd_contigs && params.design_mode in ['binder_denovo', 'binder_partialdiff', 'monomer_partialdiff']){
                // Auto-generate contigs for RFdiffusion if not provided
                println("Automatically generating RFdiffusion contigs from input PDB. Will include all residues.")
                GenerateRFDContigs(file(params.input_pdb),params.design_mode)
                GenerateRFDContigs.out.view({ contigs -> "Generated the RFdiffusion contigs: $contigs" }).set{rfdContigs}
            } else if (!params.rfd_contigs && params.design_mode == 'monomer_denovo'){
                // Contigs for monomer_denovo are equivalent to design_length
                Channel.value("[$params.design_length]").set{rfdContigs}
            }

            if(params.design_mode == 'binder_foldcond'){
                GenerateRFDfoldcond(file(params.input_pdb))
                GenerateRFDfoldcond.out.target_adj.set{target_adj}
                GenerateRFDfoldcond.out.target_ss.set{target_ss}
            } else {
                Channel.value(file("${projectDir}/lib/NO_FILE")).set{target_adj}
                Channel.value(file("${projectDir}/lib/NO_FILE1")).set{target_ss}
            }

            // Collect input files
            def inputFiles = collectInputFiles(params)

            // Copy input files to output directory
            inputFiles.each { inputFile ->
                "rsync -r ${inputFile} ${inputsDir}/.".execute()
            }
            // Create the channel for RFdiffusion
            rf_ch = Channel
                .fromList((0..<num_designs).collate(batch_size))
                .map { batch ->
                    def batchId = batch.isEmpty() ? 0 : (batch[0] / batch_size).intValue()
                    def designStartnum = batch.min()
                    tuple(
                        batchId,
                        batch_size,
                        designStartnum,
                        params.design_mode,
                        inputFiles,
                    )
                }
                .combine(rfdContigs)
                .combine(target_adj)
                .combine(target_ss)
            
            // Run RFdiffusion with the generated channel
            RunRFD(rf_ch)

            RunRFD.out.pdbs_jsons.set { rfd_pdbs_jsons }
            // Compress output files
            CompressRFD("rfd", rfd_pdbs_jsons.flatten().collect())

            // Batch RFD PDBs and JSONS for CPU tasks
            Utils
                .rebatchTuples(rfd_pdbs_jsons, 200)
                .set { fold_tuples }
        }

        // Fold filtering - secondary structure and radius of gyration
        FilterFold(fold_tuples)

        // If Running Fold Design only these are the final pdbs
        if (params.run_fold_only) {
            FilterFold.out.pdbs_jsons
                .flatten()
                .collect()
                .ifEmpty(file("${projectDir}/lib/placeholder.pdb"))
                .set { final_pdbs }
        }
        else {
            FilterFold.out.pdbs_jsons.set { filt_fold_pdbs_jsons }
        }
    }
    else if (params.skip_fold & !params.skip_fold_seq & !params.skip_fold_seq_pred) {
        // Skip Fold Design and use existing PDBs and JSONs from specified directory
        println("Skipping Fold Design stage as skip_fold=true.")
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
            .set { filt_fold_pdbs_jsons }
    }
    else {
        println("Skipping Fold Design stage as skip_fold_seq=true or skip_fold_seq_pred=true.")
    }
    ///////////////////////////
    // SEQUENCE DESIGN STAGE //
    ///////////////////////////
    // Run Sequence Design if not skipped
    if (!params.skip_fold_seq & !params.skip_fold_seq_pred & !params.run_fold_only) {
        // Sequence design (either MPNN or FAMPNN)
        if (params.seq_method == "mpnn") {
            // Add FIXED labels to PDBs for target residues so the sequence does not change
            PrepMPNN(filt_fold_pdbs_jsons)
            
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
                .rebatchTuples(filt_fold_pdbs_jsons, 10)
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

            if (params.design_mode in ['binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff']) {
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
    else if (!params.skip_fold_seq_pred & !params.run_fold_only) {
        // Skip sequence design and run prediction using existing PDBs from specified directory
        println("Skipping Sequence Design stage as skip_fold_seq=true.")
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
    else if (params.skip_fold_seq_pred) {
        println("Skipping Sequence Design stage as skip_fold_seq_pred=true.")
    }
    else {
        println("Skipping Sequence Design stage as run_fold_only=true.")
    }
    ////////////////////////////////
    // STRUCTURE PREDICTION STAGE //
    ////////////////////////////////
    // Run Structure Prediction if not skipped
    if (!params.skip_fold_seq_pred & !params.run_fold_only) {
        // Optional uncropped target PDB merge for binder design
        if (params.design_mode in ['bindcraft', 'binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff']) {
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

            if (params.design_mode in ['binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff']) {
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
            
            // Handle templates - use empty channel if not present
            PrepBoltz.out.templates
                .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                .set { templates_ch }

            // Perform prediction of designs using Boltz-2
            RunBoltz(pred_input_tuple, templates_ch)

            // Batch files for CPUs
            Utils
                .rebatchTuples(RunBoltz.out.pdbs_jsons, 200)
                .set { boltz_tuple }

            // Align Boltz Predictions to FAMPNN output and calculate RMSD
            if (params.design_mode in ['binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff']) {
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
    else if (!params.run_fold_only) {
        // Skip prediction and run analysis only using existing PDBs from specified directory
        println("Skipping Structure Prediction stage as skip_fold_seq_pred=true.")
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
        println("Skipping Structure Prediction stage as run_fold_only=true.")
    }
    ////////////////////
    // ANALYSIS STAGE //
    ////////////////////
    if (!params.run_fold_only) {
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
        println("Skipping Analysis stage as run_fold_only=true.")
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
    if (params.run_fold_only) {
        Utils.countPdbFiles(fold_tuples).set { fold_count }
        Utils.countPdbFiles(final_pdbs).set { filter_fold_count }
        seq_count = 0
        filter_seq_count = 0
        filter_pred_count = 0
    }
    else if (params.skip_fold_seq_pred) {
        fold_count = 0
        filter_fold_count = 0
        seq_count = 0
        filter_seq_count = 0
        filter_pred_count = 0
    }
    else if (params.skip_fold_seq) {
        fold_count = 0
        filter_fold_count = 0
        seq_count = 0
        filter_seq_count = 0
        Utils.countPdbFiles(analysis_input_pdbs).set { filter_pred_count }
    }
    else if (params.skip_fold) {
        fold_count = 0
        filter_fold_count = 0
        Utils.countPdbFiles(seq_tuple).set { seq_count }
        Utils.countPdbFiles(filt_seq_pdbs).set { filter_seq_count }
        Utils.countPdbFiles(analysis_input_pdbs).set { filter_pred_count }
    }
    else {
        Utils.countPdbFiles(fold_tuples).set { fold_count }
        Utils.countPdbFiles(filt_fold_pdbs_jsons).set { filter_fold_count }
        Utils.countPdbFiles(seq_tuple).set { seq_count }
        Utils.countPdbFiles(filt_seq_pdbs).set { filter_seq_count }
        Utils.countPdbFiles(analysis_input_pdbs).set { filter_pred_count }
    }

    // Generate report and statistics of run
    PublishResults(
        final_pdbs,
        all_designs_metadata,
        fold_count,
        filter_fold_count,
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
///////////////////////
// HELPER FUNCTIONS //
//////////////////////
def validateDesignLength(design_length){
    // Validate design length (required, one or two comma-separated integers, min<=max)
    if (!design_length) {
        throw new IllegalArgumentException("Please provide a value for design_length, e.g. '65', '65-150'.")
    }
    def designLengthVals = design_length.split('-')
    if (designLengthVals.size() > 2 || !designLengthVals.every { it.isInteger() }) {
        throw new IllegalArgumentException("design_length parameter must contain one or two integers (dash-separated) e.g. '65' '65-150'.")
    }
    if (designLengthVals.size() == 2){
        def minLength = Integer.parseInt(designLengthVals[0])
        def maxLength = Integer.parseInt(designLengthVals[1])
        if (minLength > maxLength || minLength < 1) {
            throw new IllegalArgumentException("design_length values must be valid: min ≤ max and min ≥ 1.")
        }
    }
}

def validateRFDParameters(params) {
    // Validate Parameters for RFdiffusion including design_mode and mode required parameters
    VALID_MODES = ['binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff', 'monomer_denovo', 'monomer_foldcond', 'monomer_motifscaff', 'monomer_partialdiff']
    if (!(params.design_mode in VALID_MODES)) {
        throw new IllegalArgumentException("Invalid design mode: ${params.design_mode}. Must be one of: ${VALID_MODES.join(', ')}")
    }

    switch (params.design_mode) {
        case 'binder_partialdiff':
            if (!params.rfd_partial_diffusion_timesteps) {
                throw new IllegalArgumentException("rfd_partial_diffusion_timesteps is required when mode is 'binder_partialdiff'")
            }
            break            
        case 'binder_foldcond':
            // Validate scaffold directory and contents
            if (!params.rfd_scaffold_dir) {
                throw new IllegalArgumentException("Please provide path to directory containing scaffold files for fold conditioning (rfd_scaffold_dir)")
            }
            def scaffoldsDir = new File(params.rfd_scaffold_dir)
            if (!scaffoldsDir.exists() || !scaffoldsDir.isDirectory()) {
                throw new IllegalArgumentException("rfd_scaffold_dir does not exist or is not a directory")
            }
            
            def ss_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_ss.pt') }
            def adj_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_adj.pt') }
            if (!ss_files || !adj_files) {
                throw new IllegalArgumentException("rfd_scaffold_dir does not contain required _ss.pt and _adj.pt files")
            }
            break
        case 'monomer_foldcond':
            if (!params.rfd_scaffold_dir) {
                throw new IllegalArgumentException("Please provide path to directory containing scaffold files for fold conditioning (rfd_scaffold_dir)")
            }

            def scaffoldsDir = new File(params.rfd_scaffold_dir)
            if (!scaffoldsDir.exists() || !scaffoldsDir.isDirectory()) {
                throw new IllegalArgumentException("rfd_scaffold_dir does not exist or is not a directory")
            }
            
            def ss_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_ss.pt') }
            def adj_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_adj.pt') }
            if (!ss_files || !adj_files) {
                throw new IllegalArgumentException("rfd_scaffold_dir does not contain required _ss.pt and _adj.pt files")
            }
            break                
        case 'monomer_partialdiff':
            if (!params.rfd_partial_diffusion_timesteps) {
                throw new IllegalArgumentException("rfd_partial_diffusion_timesteps is required when mode is 'monomer_partialdiff'")
            }
            break                
    }
}

// Collect required input files
def collectInputFiles(params) {
    def inputs = []

    // Add required input files
    if (params.design_mode in [
        'bindcraft',
        'binder_denovo',
        'binder_foldcond',
        'binder_motifscaff',
        'binder_partialdiff',
        'monomer_motifscaff',
        'monomer_partialdiff',
    ]) {
        if (params.input_pdb) {
            inputs << file(params.input_pdb)
        }
    }

    if(params.design_mode == 'bindcraft' & params.bc_advanced_json){
        inputs << file(params.bc_advanced_json)
    }

    if (params.design_mode in ['monomer_denovo', 'monomer_foldcond']) {
        // Add 'placeholder' PDB file, since RFdiffusion requires xyz coordinates
        inputs << file("${projectDir}/lib/placeholder.pdb")
    }
    if (params.design_mode in ['binder_foldcond', 'monomer_foldcond']) {
        if (params.rfd_scaffold_dir) {
            // Add scaffolds_dir and contents
            inputs << file(params.rfd_scaffold_dir)
        }
    }

    return inputs
}

def generateContigDescription(String contigs) {
    // Function to generate human-readable contig descriptions
    def cleaned = contigs.replaceAll(/[\[\]]/, '').trim()
    def chainParts = cleaned.split(/\s+/).findAll { it }
    def description = ["The contigs for RFdiffusion ${contigs} include ${chainParts.size()} chain${chainParts.size() > 1 ? 's' : ''}. RFdiffusion will:"]
    
    chainParts.each { part ->
        def segments = part.split('/')
        def chainId = null
        
        // Determine chain ID from first segment with chain designation
        segments.find { seg ->
            def matcher = seg =~ /^([A-Za-z])?(\d+)-(\d+)$/
            if (matcher.matches() && matcher.group(1)) {
                chainId = matcher.group(1)
                return true
            }
            return false
        }
        
        segments.each { seg ->
            switch (seg) {
                case '0':
                    description << "* Insert a chainbreak ${chainId ? "after chain ${chainId}" : ""}"
                    break
                default:
                    def matcher = seg =~ /^([A-Za-z])?(\d+)-(\d+)$/
                    if (matcher.matches()) {
                        def (segChain, start, end) = [matcher.group(1), matcher.group(2), matcher.group(3)]
                        if (segChain) {
                            description << "* Keep residues ${start}-${end} of chain ${segChain}"
                        } else {
                            // No chain ID - this is a partial diffusion or motifscaffolding mode
                            def count = end.toInteger() - start.toInteger() + 1
                            if (start == end) {
                                description << "* Diffuse ${start} residues${chainId ? " for chain ${chainId}" : " for a new chain"}"
                            } else {
                                description << "* Diffuse ${start}-${end} residues${chainId ? " for chain ${chainId}" : " for a new chain"}"
                            }
                        }
                    }
            }
        }
    }
    
    return description.join('\n')
}

def validateBindCraftParams(bc_chains,hotspot_residues,design_length,num_designs,input_pdb,bc_advanced_json) {

    // Validate PDB chains to target (required, comma-separated, only letters)
    if (bc_chains) {
        if (!bc_chains.matches(/^([A-Za-z]+)(,[A-Za-z]+)*$/)) {
            throw new IllegalArgumentException("bc_chains parameter must be a comma-separated list of chain identifiers, e.g. 'A,B' or 'A'.")
        }
    }

    // Validate hotspot residues (optional, allow null, otherwise enforce expected format)
    def hotspot = hotspot_residues
    if (hotspot && !hotspot.matches(/^([A-Za-z]*[0-9]+([\-][0-9]+)?)(,[A-Za-z]*[0-9]+([\-][0-9]+)?)*$|^([A-Za-z]+,?)*$/)) {
        throw new IllegalArgumentException("hotspot_residues format invalid. Acceptable: '1,2-10', 'A10,A12,B2-10', or chains 'A,B'.")
    }

    // Validate optional BindCraft advanced settings JSON
    if (bc_advanced_json) {
    def advancedFile = new File(bc_advanced_json)
    if (!advancedFile.exists()) {
        throw new FileNotFoundException("Advanced settings JSON file not found at path: ${bc_advanced_json}. Please ensure the file exists and the path is correct.")
    }
    }
}

def getAdvancedSettingsPath(bc_design_protocol, bc_template_protocol) {
    // Generate the advanced settings path for BindCraft according to protocol parameters.
    // Design protocol tag
    def designProtocolTag
    switch(bc_design_protocol) {
        case "default":
            designProtocolTag = "default_4stage_multimer"
            break
        case "beta-sheet":
            designProtocolTag = "betasheet_4stage_multimer"
            break
        case "peptide":
            designProtocolTag = "peptide_3stage_multimer"
            break
        default:
            throw new IllegalArgumentException("Unsupported BindCraft design protocol: ${bc_design_protocol}")
    }

    // Template protocol tag
    def templateProtocolTag
    switch(bc_template_protocol) {
        case "default":
            templateProtocolTag = ""
            break
        case "flexible":
            templateProtocolTag = "_flexible"
            break
        default:
            throw new IllegalArgumentException("Unsupported BindCraft template protocol: ${bc_template_protocol}")
    }

    // Compose the path
    def advancedSettingsPath = file("${projectDir}/lib/bindcraft/settings_advanced/" +
        designProtocolTag +
        templateProtocolTag +
        ".json"
    )

    log.info "Selected advanced design settings file: ${advancedSettingsPath}\n"
    return advancedSettingsPath
}
