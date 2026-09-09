#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GenerateRFDContigs; GenerateRFDFoldCond; FilterFold ; RunRFD } from './modules/rfdiffusion.nf'
include { AnalyseBC; PrepBC; RunBC } from './modules/bindcraft.nf'
include { AnalyseBG; PrepBG; RunBG } from './modules/boltzgen.nf'
include { PrepFAMPNN ; RunFAMPNN } from './modules/fampnn.nf'
include { FilterSeq; PrepMPNN ; RunMPNN } from './modules/proteinmpnn.nf'
include { AlignAF2; FilterAF2; RunAF2 } from './modules/af2.nf'
include { AnalysePredictions; FilterAnalysis } from './modules/analysis.nf'
include { PublishResults } from './modules/publish.nf'
include { AlignBoltz ; FilterBoltz; PrepBoltz ; RunBoltz; AnalyseBoltz } from './modules/boltz.nf'
include { PrepBoltz as PrepBoltzUnbound ; RunBoltz as RunBoltzUnbound ; AlignBoltz as AlignBoltzUnbound } from './modules/boltz.nf'
include { CombineMetadata } from './modules/combine_metadata.nf'
include { Compress as CompressRFD } from './modules/compress'
include { Compress as CompressMPNN } from './modules/compress'
include { Compress as CompressFAMPNN } from './modules/compress'
include { Compress as CompressAF2 } from './modules/compress'
include { Compress as CompressBoltz } from './modules/compress'
include { Compress as CompressBoltzUnbound } from './modules/compress'
include { Compress as CompressBC } from './modules/compress'
include { Compress as CompressBG } from './modules/compress'
include { MergeUncroppedTarget } from './modules/merge_uncropped_target.nf'

workflow {
    // Permit use of topic channels in Nextflow v24 by enabling preview features
    try {
        nextflow.preview.topic = true
    } catch (Exception _e) {
        // Silently continue - topic channels work by default in v25+
    }
    
    if (params.cpus < 2){
        throw new IllegalArgumentException("--cpus must be >= 2")
    }
    if (params.cpus_per_gpu < 2){
        throw new IllegalArgumentException("--cpus_per_gpu must be >= 2")
    }
    if (params.gpus < 1){
        throw new IllegalArgumentException("--gpus must be >= 1")
    }

    // Auto-create the AF2 JAX compilation cache dir since apptainer/singularity fail to bind-mount a missing path
    if (params.af2_jax_compilation_cache_dir != null) {
        def jaxCacheDir = launchDir.resolve(params.af2_jax_compilation_cache_dir.toString()).normalize()
        if (jaxCacheDir.exists() && !jaxCacheDir.isDirectory()) {
            throw new IllegalArgumentException("af2_jax_compilation_cache_dir (${jaxCacheDir}) exists but is not a directory")
        }
        jaxCacheDir.mkdirs()
        if (!jaxCacheDir.exists()) {
            throw new IllegalArgumentException("Failed to create af2_jax_compilation_cache_dir (${jaxCacheDir}). Check parent directory permissions.")
        }
    }

    def outputDirectory = params.out_dir

    VALID_MODES = ['bindcraft_denovo', 'boltzgen_denovo', 'boltzgen_motifscaff', 'rfd_denovo', 'rfd_foldcond', 'rfd_motifscaff', 'rfd_partialdiff']
    if (!(params.design_mode in VALID_MODES)) {
        throw new IllegalArgumentException("Invalid design mode: ${params.design_mode}. Must be one of: ${VALID_MODES.join(', ')}")
    }

    // Auto-detect monomer vs binder behavior for RFdiffusion/BoltzGen modes (computed once, used throughout)
    def running_fold_design = !params.skip_fold && !params.skip_fold_seq && !params.skip_fold_seq_pred
    def is_binder_mode
    if (params.design_mode == 'bindcraft_denovo') {
        is_binder_mode = true
    } else if (running_fold_design) {
        // Fold design will run this invocation - detect from RFdiffusion/BoltzGen input parameters
        is_binder_mode = detectIsBinderModeFromParams(params.design_mode, params)
    } else {
        // Fold design is being skipped this invocation - derive binder/monomer status from
        // the chain count of an actual resumed PDB in params.skip_input_dir.
        is_binder_mode = detectIsBinderModeFromResumedPdb(params.skip_input_dir)
    }

    if (params.run_fold_only && (params.skip_fold_seq || params.skip_fold_seq_pred )) {
        error("Cannot use --run_fold_only with skip flags --skip_fold_seq or --skip_fold_seq_pred. These options are contradictory.")
    }
    if (params.run_fold_only && params.skip_fold) {
        error("Cannot use --run_fold_only with --skip_fold. These options are contradictory.")
    }

    if (params.run_fold_only && params.rank_designs) {
        error("--rank_designs cannot be used with --run_fold_only since no prediction metrics are generated to rank on.")
    }

    if (params.skip_fold_seq_pred && params.rank_designs) {
        error("--rank_designs cannot be used with --skip_fold_seq_pred since no prediction metrics are generated to rank on.")
    }

    // Validate ranking metric matches prediction method
    def ranking_metric = null
    if (params.rank_designs && params.ranking_metric) {
        ranking_metric = validateranking_metric(params.ranking_metric, params.pred_method)
    } else if (params.rank_designs && !params.ranking_metric) {
        // Use default ranking metrics based on mode and prediction method
        // 'af2_boltz' defaults to boltz_* metrics since Boltz is the final/most-refined stage
        if (!is_binder_mode) {
            // For monomer modes, use overall quality metrics (no interface)
            ranking_metric = params.pred_method in ['boltz', 'af2_boltz'] ? 'boltz_ptm' : 'af2_plddt_overall'
        } else {
            // For binder modes, use interface-specific metrics
            ranking_metric = params.pred_method in ['boltz', 'af2_boltz'] ? 'boltz_ipSAE_min' : 'af2_pae_interaction'
        }
    }

    // Calculate batch size based on maximum GPUs
    def num_batches = Math.min(params.gpus, params.num_designs).intValue()
    def batch_size = Math.ceil(params.num_designs / num_batches).intValue()
    def num_designs = num_batches * batch_size

    // Pre-compute filter parameter strings here in the workflow body (to enable caching and resume)
    def fold_filter_params = Utils.formatFilterParams(params, "fold", ["min_ss", "max_ss", "min_helices", "max_helices", "min_strands", "max_strands", "min_rog", "max_rog"])
    def seq_max_score = params.seq_method == 'mpnn' ? params.mpnn_max_score : params.fampnn_max_psce
    def seq_max_score_param = seq_max_score != null ? "--max-score ${seq_max_score}" : ''
    def seq_filter_params = Utils.formatFilterParams(params, "seq", ["min_ext_coef", "max_ext_coef", "min_pi", "max_pi"])
    def af2_filter_params = Utils.formatFilterParams(params, "af2", ["min_iptm", "max_pae_interaction", "max_pae_overall", "max_pae_binder", "max_pae_target", "min_plddt_overall", "min_plddt_binder", "min_plddt_target", "max_rmsd_overall", "max_rmsd_binder_bndaln", "max_rmsd_binder_tgtaln", "max_rmsd_target"])
    def boltz_filter_params = Utils.formatFilterParams(params, "boltz", ["max_rmsd_overall", "max_rmsd_binder", "max_rmsd_target", "min_conf_score", "min_ptm", "min_ptm_binder", "min_ptm_target", "min_iptm", "min_plddt", "min_iplddt", "max_pde", "max_ipde", "min_ipSAE_min", "min_LIS", "min_pDockQ2_min", "max_pae_interaction", "max_unbound_rmsd", "min_unbound_conf_score", "min_unbound_ptm", "min_unbound_plddt", "max_unbound_pde"])
    def analysis_filter_params = Utils.formatFilterParams(params, "pr", ["min_helices", "max_helices", "min_strands", "max_strands", "min_total_ss", "max_total_ss", "min_rog", "max_rog", "min_intface_bsa", "min_intface_shpcomp", "min_intface_hbonds", "max_intface_unsat_hbonds", "max_intface_deltag", "max_intface_deltagtobsa", "max_surfhphobics", "max_sap", "max_sap_complex"])

    println("***********************************************************************")
    println("██████╗ ██████╗  ██████╗ ████████╗███████╗██╗███╗   ██╗██████╗      ██╗")
    println("██╔══██╗██╔══██╗██╔═══██╗╚══██╔══╝██╔════╝██║████╗  ██║██╔══██╗     ██║")
    println("██████╔╝██████╔╝██║   ██║   ██║   █████╗  ██║██╔██╗ ██║██║  ██║     ██║")
    println("██╔═══╝ ██╔══██╗██║   ██║   ██║   ██╔══╝  ██║██║╚██╗██║██║  ██║██   ██║")
    println("██║     ██║  ██║╚██████╔╝   ██║   ███████╗██║██║ ╚████║██████╔╝╚█████╔╝")
    println("╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚══════╝╚═╝╚═╝  ╚═══╝╚═════╝  ╚════╝ ")
    println("                   ProteinDJ Protein Design Pipeline                   ")
    println("         Developers: Dylan Silke, Josh Hardy, Julie Iskander,          ")
    println("                  Junqi Pan, Thomas Litfin, David Ladd                 ")
    println("***********************************************************************")
    println("* Pipeline Mode: ${params.design_mode}")
    println("* Number of designs: ${num_designs}")
    println("* Number of sequences for each design: ${params.seqs_per_design}")
    if (params.rank_designs) {
        println("* Ranking Metric: ${ranking_metric}")
        if (params.max_designs != null) {
            println("* Maximum designs output after ranking: ${params.max_designs}")
        }
    }
    println("* Output Directory: ${outputDirectory}")
    println("***********************************************************************\n")

    // Create output directory for copy of config files used in run
    def outputDir = file(outputDirectory)
    def configDir = outputDir.resolve('configs')
    configDir.mkdirs()
    workflow.configFiles.each { configFile ->
        def configPath = file(configFile)
        configPath.copyTo(configDir.resolve(configPath.name))
    }

    // Create output directory for copy of input files used in run
    def inputsDir = outputDir.resolve('inputs')
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

        // Dispatch on design_mode: validation, startup logging, and input-file collection are
        // handled by the engine strategy object; channel/process wiring stays below per-engine.
        def engine = DesignEngine.forMode(params, is_binder_mode)

        // Validate input PDB file
        if (engine.requiresInputPdb()) {
            if (!params.input_pdb) {
                throw new IllegalArgumentException("Please provide input PDB file path required by $params.design_mode mode")
            }
            def inputFile = file(params.input_pdb)
            if (!inputFile.exists()) {
                throw new FileNotFoundException("Input PDB file not found at path: ${params.input_pdb}. Please ensure the file exists and the path is correct.")
            }
        }

        engine.validateParams()
        engine.startupMessages().each { log.info(it) }

        if (params.design_mode=="bindcraft_denovo"){
            // Use Bindcraft for fold design

            // Select advanced settings json file    
            def bc_advanced_json 
            if(!params.bc_advanced_json){
                bc_advanced_json = file("${projectDir}/lib/bindcraft/settings_advanced/${engine.advancedSettingsRelativePath()}")
                log.info "Selected advanced design settings file: ${bc_advanced_json}\n"
            } else {
                bc_advanced_json = file(params.bc_advanced_json)
            }

            // Get path of filter settings json file
            def bc_filters_json = file("${projectDir}/lib/bindcraft/settings_filters/no_filters.json")

            // Collect input files
            def inputFiles = engine.collectInputFiles().collect { file(it) }

            // Copy input files to output directory
            inputFiles.each { inputFile ->
                inputFile.copyTo(inputsDir.resolve(inputFile.name))
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

        } else if (params.design_mode in ['boltzgen_denovo','boltzgen_motifscaff']) {
            // Use BoltzGen for fold design

            // Collect input files
            def inputFiles = engine.collectInputFiles().collect { file(it) }

            // Copy input files to output directory
            inputFiles.each { inputFile ->
                inputFile.copyTo(inputsDir.resolve(inputFile.name))
            }
            
            // Create channel with items for requested designs
            bg_ch = Channel
                .fromList((0..<num_batches))

            // No target for monomer denovo design - use the NO_FILE placeholder
            def bg_input_pdb = params.input_pdb ? file(params.input_pdb) : file("${projectDir}/lib/NO_FILE")
            PrepBG(bg_ch, bg_input_pdb, params.design_mode)

            // Run BoltzGen for each batch (PrepBG emits a cleaned copy of input_pdb alongside the design spec)
            RunBG(PrepBG.out, batch_size)

            // Collect batches and run analysis and conversion of BoltzGen outputs
            AnalyseBG(RunBG.out.cifs_npzs.flatten().collect(), params.design_mode)
            AnalyseBG.out.pdbs_jsons.set { bg_pdbs_jsons }

            // Compress output files
            CompressBG("bg", bg_pdbs_jsons.flatten().collect())

            // Batch PDBs and JSONS for CPU tasks
            Utils
                .rebatchTuples(bg_pdbs_jsons, 200)
                .set { fold_tuples }

        } else { // Use RFdiffusion for fold design
            // Check which contigs auto-generation strategy applies for this mode
            if (params.design_mode == 'rfd_foldcond'){
               Channel.value('NoContigsNeededForFoldConditioning').set{rfdContigs}
            } else if (params.design_mode == 'rfd_motifscaff'){
                // Auto-generate contigs from motifscaff_spec/target-chain detection, resolved in
                // Groovy by RFDiffusionParams at command-generation time (no PDB parsing here).
                // engine.validateParams() has already confirmed at least one of motifscaff_spec/
                // motifscaff_inpaint_seq/flexible_residues is set.
                println("Automatically generating RFdiffusion contigs from motifscaff_spec/motifscaff_inpaint_seq/flexible_residues.")
                Channel.value(RFDiffusionParams.AUTO_CONTIGS_SENTINEL).set{rfdContigs}
            } else if (params.design_mode == 'rfd_partialdiff' && params.rfd_partialdiff_spec){
                // Auto-generate contigs from rfd_partialdiff_spec/target-chain detection, resolved in
                // Groovy by RFDiffusionParams at command-generation time (no PDB parsing here).
                println("Automatically generating RFdiffusion contigs from rfd_partialdiff_spec.")
                Channel.value(RFDiffusionParams.AUTO_PARTIALDIFF_SENTINEL).set{rfdContigs}
            } else if (params.design_mode == 'rfd_denovo' && !is_binder_mode){
                // Contigs for monomer rfd_denovo are equivalent to design_length. Expand a bare
                // single-value design_length (e.g. '65') to 'N-N'
                Channel.value("[${Utils.designLengthContigRange(params.design_length)}]").set{rfdContigs}
            } else {
                // Auto-generate contigs for RFdiffusion if not provided (rfd_denovo binder, or rfd_partialdiff monomer/binder without rfd_partialdiff_spec)
                println("Automatically generating RFdiffusion contigs from input PDB. Will include all residues.")
                GenerateRFDContigs(file(params.input_pdb), params.design_mode, is_binder_mode)
                GenerateRFDContigs.out.view({ contigs -> "Generated the RFdiffusion contigs: $contigs" }).set{rfdContigs}
            }

            if(params.design_mode == 'rfd_foldcond' && is_binder_mode){
                GenerateRFDFoldCond(file(params.input_pdb))
                GenerateRFDFoldCond.out.target_adj.set{target_adj}
                GenerateRFDFoldCond.out.target_ss.set{target_ss}
            } else {
                Channel.value(file("${projectDir}/lib/NO_FILE")).set{target_adj}
                Channel.value(file("${projectDir}/lib/NO_FILE1")).set{target_ss}
            }

            // Collect input files
            def inputFiles = engine.collectInputFiles().collect { file(it) }
            if (engine.needsPlaceholderPdb()) {
                inputFiles << file("${projectDir}/lib/placeholder.pdb")
            }

            // Copy input files to output directory
            inputFiles.each { inputFile ->
                inputFile.copyTo(inputsDir.resolve(inputFile.name))
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
                .map { batchId, batchSizeVal, designStartnum, mode, files, contigs, adj, ss ->
                    def rfdParams = new RFDiffusionParams(params + [is_binder_mode: is_binder_mode])
                    def rfdCommand = rfdParams.generateCommandString(contigs)
                    tuple(batchId, batchSizeVal, designStartnum, mode, files, rfdCommand, adj, ss)
                }
            
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
        FilterFold(fold_tuples, fold_filter_params)

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
        def previous_pdbs = files(file(params.skip_input_dir).resolve('*.pdb'))
        def previous_jsons = files(file(params.skip_input_dir).resolve('*.json'))
        // Error handling for missing files
        if (previous_pdbs.isEmpty()) {
            throw new FileNotFoundException("No PDB files found in directory: ${params.skip_input_dir}. Please provide PDB files to proceed with the workflow.")
        }
        if (previous_jsons.isEmpty()) {
            throw new FileNotFoundException("No JSON files found in directory: ${params.skip_input_dir}. Please provide JSON files to proceed with the workflow.")
        }
        println("Found ${previous_pdbs.size()} PDB files")
        println("Found ${previous_jsons.size()} JSON files\n")

        // Validate naming convention
        validateFileNaming(previous_pdbs, previous_jsons, 'fold_only')

        // Copy PDB and JSON files from the previous results directory to inputs directory
        previous_pdbs.each { pdbFile ->
            pdbFile.copyTo(inputsDir.resolve(pdbFile.name))
        }
        previous_jsons.each { jsonFile ->
            jsonFile.copyTo(inputsDir.resolve(jsonFile.name))
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
            
            // GPU-aware batching for RunMPNN
            Utils
                .rebatchGPU(PrepMPNN.out.pdbs, params.gpus)
                .set { seq_input_pdbs }

            // Launch ProteinMPNN
            RunMPNN(seq_input_pdbs)

            // Compress output files
            CompressMPNN("mpnn", RunMPNN.out.pdbs_jsons.flatten().collect())

            // Rebatch sequence assignment files for CPU Filtering Step
            Utils
                .rebatchTuples(RunMPNN.out.pdbs_jsons, 200)
                .set { seq_tuple }
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

            if (is_binder_mode) {
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
        }
        else {
            error("Not a valid sequence assignment method")
        }

        // Filter designs by sequence score
        FilterSeq(seq_tuple, seq_max_score_param, seq_filter_params)
        FilterSeq.out.pdbs
            .flatten()
            .collect()
            .set { filt_seq_pdbs }
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
        def pdbs_for_pred = files(file(params.skip_input_dir).resolve('*.pdb'))
        if (pdbs_for_pred.isEmpty()) {
            throw new FileNotFoundException("No PDB files found in directory: ${params.skip_input_dir}. Please provide PDB files to proceed with the workflow.")
        }
        println("Found ${pdbs_for_pred.size()} PDB files")

        // Validate naming convention
        validateFileNaming(pdbs_for_pred, null, 'fold_seq')

        // Copy PDB files from the previous results directory to inputs directory
        pdbs_for_pred.each { pdbFile ->
            pdbFile.copyTo(inputsDir.resolve(pdbFile.name))
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
        if (is_binder_mode) {
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
                .set { pred_tuple }

            // Filtering of AF2 results
            FilterAF2(pred_tuple, af2_filter_params)

            if (is_binder_mode) {
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
            // Handle MSA file input
            if (params.boltz_input_msa) {
                msa_input = file(params.boltz_input_msa, checkIfExists: true)
            } else {
                msa_input = file("${projectDir}/lib/NO_FILE")
            }

            // Prep yaml files for Boltz-2
            PrepBoltz(pred_input_pdbs, msa_input, is_binder_mode, false)

            // Handle templates - use empty channel if not present
            PrepBoltz.out.templates
                .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                .set { templates_ch }

            // Handle MSA file - use empty channel if not present
            PrepBoltz.out.msa_file
                .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                .set { msa_ch }

            // reallocate batching for GPU
            Utils
                .rebatchGPU(PrepBoltz.out.yamls, params.gpus)
                .combine(templates_ch)
                .combine(msa_ch)
                .set { pred_input_tuple }

            // Perform prediction of designs using Boltz-2
            RunBoltz(pred_input_tuple)

            // Batch files for CPUs
            Utils
                .rebatchTuples(RunBoltz.out.pdbs_jsons, 200)
                .set { pred_tuple }

            // Convert npz files to value channel for reuse across all batches
            RunBoltz.out.npzs.collect().set { npz_files_for_analysis }

            // Convert pred_input_pdbs to value channel for reuse across all batches
            pred_input_pdbs.collect().set { designs_for_alignment }

            // Calculate Boltz-2 interface scores for binders only
            if (is_binder_mode) {
                AnalyseBoltz(pred_tuple, npz_files_for_analysis)
                AnalyseBoltz.out.pdbs_jsons
                    .set { boltz_with_metrics }
            } else{
                pred_tuple.set { boltz_with_metrics }
            }

            // Align Boltz Predictions to FAMPNN output and calculate RMSD
            if (is_binder_mode) {
                AlignBoltz(boltz_with_metrics, designs_for_alignment, 'binder')
            }
            else {
                AlignBoltz(boltz_with_metrics, designs_for_alignment, 'monomer')
            }

            // Optional unbound (target-free) binder-only prediction, for QC of fold consistency
            // in the absence of the target. Runs as a fully separate Boltz-2 prediction of chain A only.
            def unbound_json_files
            if (is_binder_mode && params.boltz_predict_unbound_binder) {
                PrepBoltzUnbound(pred_input_pdbs, file("${projectDir}/lib/NO_FILE"), is_binder_mode, true)

                PrepBoltzUnbound.out.templates
                    .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                    .set { unbound_templates_ch }

                PrepBoltzUnbound.out.msa_file
                    .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                    .set { unbound_msa_ch }

                Utils
                    .rebatchGPU(PrepBoltzUnbound.out.yamls, params.gpus)
                    .combine(unbound_templates_ch)
                    .combine(unbound_msa_ch)
                    .set { unbound_pred_input_tuple }

                RunBoltzUnbound(unbound_pred_input_tuple)

                Utils
                    .rebatchTuples(RunBoltzUnbound.out.pdbs_jsons, 200)
                    .set { unbound_pred_tuple }

                // Align to the original (binder+target) design, comparing chain A only
                AlignBoltzUnbound(unbound_pred_tuple, designs_for_alignment, 'unbound_binder')

                CompressBoltzUnbound("boltz_unbound", AlignBoltzUnbound.out.pdbs_jsons.flatten().collect())

                AlignBoltzUnbound.out.pdbs_jsons
                    .map { pdbs, jsons -> jsons }
                    .flatten()
                    .collect()
                    .set { unbound_json_files_ch }
                unbound_json_files = unbound_json_files_ch
            } else {
                unbound_json_files = file("${projectDir}/lib/NO_FILE")
            }

            // Compress output files
            CompressBoltz("boltz", AlignBoltz.out.pdbs_jsons.flatten().collect())

            // Filtering of Boltz-2 results
            FilterBoltz(AlignBoltz.out.pdbs_jsons, boltz_filter_params, unbound_json_files)
            FilterBoltz.out.pdbs
                .flatten()
                .collect()
                .set { analysis_input_pdbs }
        }
        else if (params.pred_method == "af2_boltz") {
            // Cascaded prediction: run AF2 first (fast pre-filter), then re-feed the
            // original design PDBs of AF2 survivors into Boltz-2

            // --- Stage 1: AF2 predict + filter ---

            // reallocate batching for GPU
            Utils
                .rebatchGPUByNumRes(pred_input_pdbs, params.gpus)
                .set { af2c_pred_input_tuple }

            // AlphaFold2-Initial Guess
            RunAF2(af2c_pred_input_tuple)

            // Compress output files
            CompressAF2("af2", RunAF2.out.pdbs_jsons.flatten().collect())

            // Batch files for CPUs
            Utils
                .rebatchTuples(RunAF2.out.pdbs_jsons, 200)
                .set { af2c_pred_tuple }

            // Filtering of AF2 results
            FilterAF2(af2c_pred_tuple, af2_filter_params)

            if (is_binder_mode) {
                // Alignment of PDBs to target chain(s). Only need one reference file
                AlignAF2(FilterAF2.out.pdbs.flatten().collect(), pred_input_pdbs.flatten().last())
                AlignAF2.out.pdbs
                    .flatten()
                    .collect()
                    .set { af2c_survivor_pdbs }
            } else {
                FilterAF2.out.pdbs
                    .flatten()
                    .collect()
                    .set { af2c_survivor_pdbs }
            }

            // --- Join AF2 survivors back to their original (pre-prediction) design PDBs ---
            // AF2 renames outputs by appending '_af2pred' to the basename (e.g. 'fold_0_seq_1.pdb'
            // -> 'fold_0_seq_1_af2pred.pdb'), so strip that suffix before matching on filename.
            pred_input_pdbs
                .flatten()
                .map { pdb -> tuple(pdb.getName(), pdb) }
                .join(
                    af2c_survivor_pdbs
                        .flatten()
                        .map { pdb -> tuple(pdb.getName().replaceFirst(/_af2pred\.pdb$/, '.pdb'), pdb) }
                )
                .map { name, orig_pdb, survivor_pdb -> orig_pdb }
                .collect()
                .set { boltz_input_pdbs }

            // --- Stage 2: Boltz-2 predict + filter on AF2 survivors ---

            // Handle MSA file input
            if (params.boltz_input_msa) {
                msa_input = file(params.boltz_input_msa, checkIfExists: true)
            } else {
                msa_input = file("${projectDir}/lib/NO_FILE")
            }

            // Prep yaml files for Boltz-2
            PrepBoltz(boltz_input_pdbs, msa_input, is_binder_mode, false)

            // Handle templates - use empty channel if not present
            PrepBoltz.out.templates
                .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                .set { templates_ch }

            // Handle MSA file - use empty channel if not present
            PrepBoltz.out.msa_file
                .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                .set { msa_ch }

            // reallocate batching for GPU
            Utils
                .rebatchGPU(PrepBoltz.out.yamls, params.gpus)
                .combine(templates_ch)
                .combine(msa_ch)
                .set { pred_input_tuple }

            // Perform prediction of designs using Boltz-2
            RunBoltz(pred_input_tuple)

            // Batch files for CPUs
            Utils
                .rebatchTuples(RunBoltz.out.pdbs_jsons, 200)
                .set { pred_tuple }

            // Convert npz files to value channel for reuse across all batches
            RunBoltz.out.npzs.collect().set { npz_files_for_analysis }

            // Convert boltz_input_pdbs to value channel for reuse across all batches
            boltz_input_pdbs.collect().set { designs_for_alignment }

            // Calculate Boltz-2 interface scores for binders only
            if (is_binder_mode) {
                AnalyseBoltz(pred_tuple, npz_files_for_analysis)
                AnalyseBoltz.out.pdbs_jsons
                    .set { boltz_with_metrics }
            } else{
                pred_tuple.set { boltz_with_metrics }
            }

            // Align Boltz Predictions to FAMPNN output and calculate RMSD
            if (is_binder_mode) {
                AlignBoltz(boltz_with_metrics, designs_for_alignment, 'binder')
            }
            else {
                AlignBoltz(boltz_with_metrics, designs_for_alignment, 'monomer')
            }

            // Optional unbound (target-free) binder-only prediction, for QC of fold consistency
            // in the absence of the target. Runs as a fully separate Boltz-2 prediction of chain A only.
            def unbound_json_files
            if (is_binder_mode && params.boltz_predict_unbound_binder) {
                PrepBoltzUnbound(boltz_input_pdbs, file("${projectDir}/lib/NO_FILE"), is_binder_mode, true)

                PrepBoltzUnbound.out.templates
                    .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                    .set { unbound_templates_ch }

                PrepBoltzUnbound.out.msa_file
                    .ifEmpty(file("${projectDir}/lib/NO_FILE"))
                    .set { unbound_msa_ch }

                Utils
                    .rebatchGPU(PrepBoltzUnbound.out.yamls, params.gpus)
                    .combine(unbound_templates_ch)
                    .combine(unbound_msa_ch)
                    .set { unbound_pred_input_tuple }

                RunBoltzUnbound(unbound_pred_input_tuple)

                Utils
                    .rebatchTuples(RunBoltzUnbound.out.pdbs_jsons, 200)
                    .set { unbound_pred_tuple }

                // Align to the original (binder+target) design, comparing chain A only
                AlignBoltzUnbound(unbound_pred_tuple, designs_for_alignment, 'unbound_binder')

                CompressBoltzUnbound("boltz_unbound", AlignBoltzUnbound.out.pdbs_jsons.flatten().collect())

                AlignBoltzUnbound.out.pdbs_jsons
                    .map { pdbs, jsons -> jsons }
                    .flatten()
                    .collect()
                    .set { unbound_json_files_ch }
                unbound_json_files = unbound_json_files_ch
            } else {
                unbound_json_files = file("${projectDir}/lib/NO_FILE")
            }

            // Compress output files
            CompressBoltz("boltz", AlignBoltz.out.pdbs_jsons.flatten().collect())

            // Filtering of Boltz-2 results
            FilterBoltz(AlignBoltz.out.pdbs_jsons, boltz_filter_params, unbound_json_files)
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
        def pdbs_for_analysis = files(file(params.skip_input_dir).resolve('*.pdb'))
        if (pdbs_for_analysis.isEmpty()) {
            throw new FileNotFoundException("No PDB files found in directory: ${params.skip_input_dir}. Please provide PDB files to proceed with the workflow.")
        }
        println("Found ${pdbs_for_analysis.size()} PDB files")

        // Validate naming convention
        validateFileNaming(pdbs_for_analysis, null, 'fold_seq_pred')

        // Copy PDB files from the previous results directory to inputs directory
        pdbs_for_analysis.each { pdbFile ->
            pdbFile.copyTo(inputsDir.resolve(pdbFile.name))
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
        AnalysePredictions(analysis_input_pdbs)

        // Filtering of analysis results
        FilterAnalysis(AnalysePredictions.out.jsonl, AnalysePredictions.out.relaxed_pdbs, analysis_filter_params)

        // Use placeholder PDB file if no designs survive filtering
        FilterAnalysis.out.pdbs
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

    // Count outputs. Each guard below mirrors the entry condition of the stage that
    // populated the corresponding channel earlier in the workflow, so referencing a
    // not-yet-set channel is avoided by ternary short-circuiting rather than branching.
    fold_count = running_fold_design ? Utils.countPdbFiles(fold_tuples) : 0
    filter_fold_count = running_fold_design
        ? (params.run_fold_only ? Utils.countPdbFiles(final_pdbs) : Utils.countPdbFiles(filt_fold_pdbs_jsons))
        : 0

    def ran_seq = !params.skip_fold_seq && !params.skip_fold_seq_pred && !params.run_fold_only
    seq_count = ran_seq ? Utils.countPdbFiles(seq_tuple) : 0
    filter_seq_count = ran_seq ? Utils.countPdbFiles(filt_seq_pdbs) : 0

    def ran_pred = !params.skip_fold_seq_pred && !params.run_fold_only
    pred_count = ran_pred ? Utils.countPdbFiles(pred_tuple) : 0

    filter_pred_count = !params.run_fold_only ? Utils.countPdbFiles(analysis_input_pdbs) : 0
    filter_analysis_count = !params.run_fold_only ? Utils.countPdbFiles(FilterAnalysis.out.pdbs) : 0

    def is_af2_boltz_cascade = ran_pred && params.pred_method == 'af2_boltz'
    af2_count = is_af2_boltz_cascade ? Utils.countPdbFiles(af2c_pred_tuple) : 0
    af2_filter_count = is_af2_boltz_cascade ? Utils.countPdbFiles(af2c_survivor_pdbs) : 0
    boltz_count = is_af2_boltz_cascade ? Utils.countPdbFiles(pred_tuple) : 0
    boltz_filter_count = is_af2_boltz_cascade ? Utils.countPdbFiles(analysis_input_pdbs) : 0

    // Generate report and statistics of run
    PublishResults(
        final_pdbs,
        all_designs_metadata,
        fold_count,
        filter_fold_count,
        seq_count,
        filter_seq_count,
        pred_count,
        filter_pred_count,
        filter_analysis_count,
        af2_count,
        af2_filter_count,
        boltz_count,
        boltz_filter_count,
        is_binder_mode
    )

    // Real post-ranking count, only known once rank_designs.py has actually run
    PublishResults.out.final_designs_count
        .view { count -> "* Final designs kept after ranking: ${count}\n" }
    
    // Save log file on completion
    workflow.onComplete {
        def logFile = file('.nextflow.log')
        if (logFile.exists()) {
            logFile.copyTo(outputDir.resolve('nextflow.log'))
        }
    }
}
///////////////////////
// HELPER FUNCTIONS //
//////////////////////
def validateranking_metric(ranking_metric, pred_method) {
    // Validate that ranking_metric matches the prediction method
    if (pred_method == 'af2' && !ranking_metric.startsWith('af2_')) {
        throw new IllegalArgumentException(
            "Ranking metric '${ranking_metric}' does not match prediction method '${pred_method}'. " +
            "For AlphaFold2 predictions, use metrics with 'af2_' prefix (e.g., 'af2_pae_interaction', 'af2_plddt_overall')."
        )
    }
    if (pred_method == 'boltz' && !ranking_metric.startsWith('boltz_')) {
        throw new IllegalArgumentException(
            "Ranking metric '${ranking_metric}' does not match prediction method '${pred_method}'. " +
            "For Boltz-2 predictions, use metrics with 'boltz_' prefix (e.g., 'boltz_iptm', 'boltz_ipSAE_min', 'boltz_LIS')."
        )
    }
    if (pred_method == 'af2_boltz' && !(ranking_metric.startsWith('af2_') || ranking_metric.startsWith('boltz_'))) {
        throw new IllegalArgumentException(
            "Ranking metric '${ranking_metric}' does not match prediction method '${pred_method}'. " +
            "For the combined AF2 + Boltz-2 cascade, use metrics with either the 'af2_' prefix " +
            "(e.g., 'af2_pae_interaction') or the 'boltz_' prefix (e.g., 'boltz_ipSAE_min')."
        )
    }
    return ranking_metric
}

// Auto-detect binder vs monomer behavior for RFdiffusion/BoltzGen modes when fold design will
// actually run this invocation. Uses params.input_pdb to determine chain count.
def detectIsBinderModeFromParams(design_mode, params) {
   if (design_mode in ['rfd_denovo', 'rfd_foldcond', 'boltzgen_denovo']) {
        return params.input_pdb as boolean
    }
    // rfd_motifscaff / rfd_partialdiff / boltzgen_motifscaff always require an input PDB
    if (!params.input_pdb) {
        throw new IllegalArgumentException("input_pdb is required for '${design_mode}' mode.")
    }
    def inputFile = file(params.input_pdb)
    if (!inputFile.exists()) {
        throw new FileNotFoundException("Input PDB file not found at path: ${params.input_pdb}. Please ensure the file exists and the path is correct.")
    }
    Set chainIds = Utils.getPdbChainIds(inputFile)
    if (chainIds.isEmpty()) {
        throw new IllegalArgumentException("Could not determine chain(s) from input_pdb for '${design_mode}' mode.")
    }
    return chainIds.size() >= 2
}

// Auto-detect binder vs monomer behavior when fold design is being skipped this invocation.
// Derived from the chain count of an actual resumed PDB in params.skip_input_dir, without
// requiring/validating RFdiffusion-specific input parameters (input_pdb).
def detectIsBinderModeFromResumedPdb(skip_input_dir) {
    if (!skip_input_dir || !file(skip_input_dir).exists()) {
        throw new FileNotFoundException("skip_input_dir not found at path: ${skip_input_dir}. Please ensure the path is correct.")
    }
    def pdbs = files(file(skip_input_dir).resolve('*.pdb'))
    if (!pdbs) {
        throw new FileNotFoundException("No PDB files found in directory: ${skip_input_dir}. Cannot auto-detect monomer/binder design mode.")
    }
    def firstPdb = pdbs instanceof List ? pdbs[0] : pdbs
    return Utils.getPdbChainIds(firstPdb).size() >= 2
}

def validateFileNaming(pdbFiles, jsonFiles = null, validationType = 'fold_seq') {
    /**
     * Validates PDB and optional JSON files against naming conventions
     * 
     * @param pdbFiles List of PDB files to validate
     * @param jsonFiles List of JSON files to validate (optional, for 'fold_only' type)
     * @param validationType One of: 'fold_only', 'fold_seq', 'fold_seq_pred'
     * @return Map with validated files or throws IllegalArgumentException
     */
    
    def validationResults = [:]
    def invalidPdbs = []
    def invalidJsons = []
    def pdbPattern
    def jsonPattern
    def pdbIndices = [] as Set
    def jsonIndices = [] as Set
    
    // Define patterns based on validation type
    if (validationType == 'fold_only') {
        // Pattern: fold_x.pdb and fold_x.json
        pdbPattern = ~/^fold_(\d+)\.pdb$/
        jsonPattern = ~/^fold_(\d+)\.json$/

        // Validate PDBs
        pdbFiles.each { pdbFile ->
            def matcher = pdbFile.name =~ pdbPattern
            if (matcher.matches()) {
                pdbIndices.add(matcher[0][1] as Integer)
            } else {
                invalidPdbs.add(pdbFile.name)
            }
        }

        // Validate JSONs
        if (jsonFiles) {
            jsonFiles.each { jsonFile ->
                def matcher = jsonFile.name =~ jsonPattern
                if (matcher.matches()) {
                    jsonIndices.add(matcher[0][1] as Integer)
                } else {
                    invalidJsons.add(jsonFile.name)
                }
            }
        }

        // Report PDB errors
        if (!invalidPdbs.isEmpty()) {
            def errorMsg = "Invalid PDB filename(s) detected. Files must follow the naming convention 'fold_x.pdb' where x is an integer.\n"
            errorMsg += "Invalid files found:"
            invalidPdbs.each { errorMsg += " ${it}" }
            errorMsg += "\nExample valid names: fold_0.pdb, fold_1.pdb, fold_10.pdb"
            throw new IllegalArgumentException(errorMsg)
        }

        // Report JSON errors
        if (!invalidJsons.isEmpty()) {
            def errorMsg = "Invalid JSON filename(s) detected. Files must follow the naming convention 'fold_x.json' where x is an integer.\n"
            errorMsg += "Invalid files found: "
            invalidJsons.each { errorMsg += " ${it}" }
            errorMsg += "\nExample valid names: fold_0.json, fold_1.json, fold_10.json"
            throw new IllegalArgumentException(errorMsg)
        }

        // Validate pairing
        def missingJsons = pdbIndices - jsonIndices
        def missingPdbs = jsonIndices - pdbIndices

        if (!missingJsons.isEmpty() || !missingPdbs.isEmpty()) {
            def errorMsg = "Mismatch between PDB and JSON files. Each fold_x.pdb must have a corresponding fold_x.json file.\n"
            if (!missingJsons.isEmpty()) {
                errorMsg += "PDB files missing corresponding JSON files:\n"
                missingJsons.sort().each { errorMsg += "  - fold_${it}.pdb (missing fold_${it}.json)\n" }
            }
            if (!missingPdbs.isEmpty()) {
                errorMsg += "JSON files missing corresponding PDB files:\n"
                missingPdbs.sort().each { errorMsg += "  - fold_${it}.json (missing fold_${it}.pdb)\n" }
            }
            throw new IllegalArgumentException(errorMsg)
        }

        println("All PDB and JSON files passed naming validation")
        println("Found ${pdbIndices.size()} properly paired fold files\n")

    } else if (validationType == 'fold_seq') {
        // Pattern: fold_x_seq_y.pdb
        pdbPattern = ~/^fold_\d+_seq_\d+\.pdb$/

        pdbFiles.each { pdbFile ->
            if (!(pdbFile.name =~ pdbPattern)) {
                invalidPdbs.add(pdbFile.name)
            }
        }

        if (!invalidPdbs.isEmpty()) {
            def errorMsg = "Invalid PDB filename(s) detected. Files must follow the naming convention 'fold_x_seq_y.pdb' where x and y are integers.\n"
            errorMsg += "Invalid files found:"
            invalidPdbs.each { errorMsg += " ${it}" }
            errorMsg += "\nExample valid names: fold_0_seq_1.pdb, fold_10_seq_25.pdb"
            throw new IllegalArgumentException(errorMsg)
        }

        println("All PDB files passed naming validation (fold_x_seq_y.pdb)")

    } else if (validationType == 'fold_seq_pred') {
        // Pattern: fold_x_seq_y_*.pdb (with any suffix after the last underscore)
        pdbPattern = ~/^fold_\d+_seq_\d+_.+\.pdb$/

        pdbFiles.each { pdbFile ->
            if (!(pdbFile.name =~ pdbPattern)) {
                invalidPdbs.add(pdbFile.name)
            }
        }

        if (!invalidPdbs.isEmpty()) {
            def errorMsg = "Invalid PDB filename(s) detected. Files must follow the naming convention 'fold_x_seq_y_*.pdb' where x and y are integers and * is any suffix.\n"
            errorMsg += "Invalid files found:"
            invalidPdbs.each { errorMsg += " ${it}" }
            errorMsg += "\nExample valid names: fold_0_seq_1_af2pred.pdb, fold_10_seq_25_boltzpred.pdb"
            throw new IllegalArgumentException(errorMsg)
        }

        println("All PDB files passed naming validation (fold_x_seq_y_*.pdb)")

    } else {
        throw new IllegalArgumentException("Invalid validation type: ${validationType}. Must be one of: 'fold_only', 'fold_seq', 'fold_seq_pred'")
    }

    return validationResults
}
