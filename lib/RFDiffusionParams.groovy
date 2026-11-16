import java.nio.file.Paths

public class RFDiffusionParams extends HashMap<String, Object> {
    static final String RFD_SCRIPT_PATH = "/app/RFdiffusion/scripts/run_inference.py"

    // Sentinel passed as the contigs channel value for rfd_motifscaff when contigs should be
    // auto-generated from motifscaff_spec/motifscaff_inpaint_seq/flexible_residues (rather than a real
    // pre-computed contig string), since RFDiffusionParams has direct access to input_pdb/params
    // and can resolve them itself at command-generation time - mirrors 'NoContigsNeededForFoldConditioning'.
    static final String AUTO_CONTIGS_SENTINEL = 'AutoGenerateFromMotifScaffSpec'

    // Sentinel passed as the contigs channel value for rfd_partialdiff when contigs should be
    // auto-generated from rfd_partialdiff_spec (rather than a real pre-computed contig string).
    static final String AUTO_PARTIALDIFF_SENTINEL = 'AutoGenerateFromPartialDiffSpec'
    
    static final String MODEL_BASE_PATH = "/app/RFdiffusion/models/"
    static final Map<String, String> MODEL_NAMES = [
        'base': 'Base',
        'complex_base': 'Complex_base',
        'complex_fold_base': 'Complex_Fold_base',
        'inpaint_seq': 'InpaintSeq',
        'inpaint_seq_fold': 'InpaintSeq_Fold',
        'active_site': 'ActiveSite',
        'base_epoch8': 'Base_epoch8',
        'complex_beta': 'Complex_beta'
    ]   

    RFDiffusionParams(Map params) {
        super(params)
    }

    // Validate hotspot input
    def validateHotspots(String contigs, String hotspots) {
        // Parse contigs into chain:ranges map
        def chainRanges = [:]
        contigs.replaceAll(/[\[\]]/, '').split().each { seg ->
            seg.split('/').each { part ->
                def match = (part =~ /^([A-Z])?(\d+)-(\d+)$/)
                if (match) {
                    def chain = match[0][1]
                    def (start, end) = [match[0][2], match[0][3]]*.toInteger()
                    chainRanges[chain] = chainRanges.getOrDefault(chain, []) << [start, end]
                }
            }
        }
        // Hotspot validation
        hotspots.replaceAll(/[\[\]]/, '').split(/,\s*/).each { hs ->
            def match = (hs =~ /^([A-Z])?(\d+)$/)
            def chain = match[0][1]
            def res = match[0][2].toInteger()
            // First check if chain exists
            if (!chainRanges.containsKey(chain)) {
                throw new IllegalArgumentException("Hotspot $hs invalid: Chain $chain is not in contigs")
            }
            // Then check residue range
            if (!chainRanges[chain].any { res in it[0]..it[1] }) {
                throw new IllegalArgumentException("Hotspot $hs invalid: Residue $res not in chain $chain ranges ${chainRanges[chain]}")
            }
        }
    }
    
    // Expand a hotspot_residues spec (individual residues, ranges, and/or bare chain IDs) into a
    // flat, comma-separated list of individual chain+residue tokens, as required by RFdiffusion's
    // ppi.hotspot_res. A chain identifier is required per token (consistent with BindCraft/BoltzGen).
    private String expandHotspots(String hotspots) {
        if (!hotspots.matches(Utils.RESIDUE_SPEC_REGEX)) {
            throw new IllegalArgumentException("hotspot_residues format invalid. Acceptable: 'A56,A115-120,B' (chain identifier is required).")
        }
        def pdbFile = new File(this.input_pdb.toString())
        def expanded = []
        hotspots.split(',').each { token ->
            def matcher = (token =~ /^([A-Za-z]+)(\d+)?(?:-(\d+))?$/)
            matcher.matches()
            def (chain, start, end) = [matcher.group(1), matcher.group(2), matcher.group(3)]
            if (!start) {
                // Bare chain ID - expand to every residue present for that chain in the input PDB
                Utils.getPdbChainResidueNumbers(pdbFile, chain).each { resNum -> expanded << "${chain}${resNum}" }
            } else if (end) {
                // Range - expand to every residue number in the range (inclusive)
                (start.toInteger()..end.toInteger()).each { resNum -> expanded << "${chain}${resNum}" }
            } else {
                // Individual residue - keep as-is
                expanded << token
            }
        }
        return expanded.join(',')
    }

    private String getSimpleFileName(String filePath) {
        return Paths.get(filePath).getFileName().toString()
    }

    // Format a chain's residues as RFdiffusion full-chain contig notation with chain-letter prefix
    // and breaks at gaps, e.g. 'A1-40/A45-88' (single residues without a dash, e.g. 'A50').
    private String formatChainRanges(File pdbFile, String chainId) {
        def resNums = Utils.getPdbChainResidueNumbers(pdbFile, chainId)
        def ranges = []
        def rangeStart = null
        def rangeEnd = null
        resNums.each { resNum ->
            if (rangeStart == null) {
                rangeStart = resNum
                rangeEnd = resNum
            } else if (resNum == rangeEnd + 1) {
                rangeEnd = resNum
            } else {
                ranges << "${chainId}${rangeStart}-${rangeEnd}"
                rangeStart = resNum
                rangeEnd = resNum
            }
        }
        if (rangeStart != null) {
            ranges << "${chainId}${rangeStart}-${rangeEnd}"
        }
        return ranges.join('/')
    }

    // Translate a comma-separated keep/diffuse spec for the design chain into RFdiffusion's '/'-joined
    // segments. In RFdiffusion contigs, '/' joins segments within ONE chain while whitespace separates
    // chains, so commas must become '/' (never spaces). RFdiffusion also requires every fixed residue
    // segment to be an explicit range, so a single residue 'A60' is written 'A60-60'.
    private String specToContigSegments(String spec) {
        return spec.split(',').collect { token ->
            token =~ /^[A-Za-z]\d+$/ ? "${token}-${token.substring(1)}" : token
        }.join('/')
    }

    // Auto-generate an rfd_motifscaff contigmap.contigs string from motifscaff_spec (BoltzGen-style
    // keep/insert grammar, comma-separated) plus the input PDB's chains. Target chain(s), in binder mode, are
    // auto-detected and appended with full-range notation joined by '/0', exactly as the existing
    // rfd_denovo/rfd_partialdiff binder auto-contig generation does.
    private String resolveMotifScaffContigs() {
        def pdbFile = new File(this.input_pdb.toString())
        def designChain = this.is_binder_mode ? 'A' : Utils.getPdbChainIds(pdbFile).sort().first()
        def chainASegment = this.motifscaff_spec ? specToContigSegments(this.motifscaff_spec) : formatChainRanges(pdbFile, designChain)
        if (!this.is_binder_mode) {
            return "[${chainASegment}]"
        }
        def targetSegments = Utils.getPdbChainIds(pdbFile).findAll { it != designChain }.sort().collect { chain ->
            formatChainRanges(pdbFile, chain)
        }
        return "[${([chainASegment] + targetSegments).join('/0 ')}]"
    }

    // Auto-generate an rfd_partialdiff contigmap.contigs string from rfd_partialdiff_spec (same
    // keep/diffuse grammar as motifscaff_spec, but diffuse tokens are fixed counts only since partial
    // diffusion cannot change chain A's length) plus the input PDB's chains. Target chain(s), in binder
    // mode, are auto-detected and appended with full-range notation joined by '/0', exactly as
    // resolveMotifScaffContigs does.
    private String resolvePartialDiffContigs() {
        def pdbFile = new File(this.input_pdb.toString())
        def designChain = this.is_binder_mode ? 'A' : Utils.getPdbChainIds(pdbFile).sort().first()
        def chainASegment = specToContigSegments(this.rfd_partialdiff_spec)
        if (!this.is_binder_mode) {
            return "[${chainASegment}]"
        }
        def targetSegments = Utils.getPdbChainIds(pdbFile).findAll { it != designChain }.sort().collect { chain ->
            formatChainRanges(pdbFile, chain)
        }
        return "[${([chainASegment] + targetSegments).join('/0 ')}]"
    }

    // Translate a motifscaff_inpaint_seq/flexible_residues comma-token list (RESIDUE_SPEC_REGEX or
    // DESIGN_CHAIN_SPEC_REGEX grammar) into RFdiffusion's bracket/slash contigmap.inpaint_seq /
    // contigmap.inpaint_str format, e.g. 'A10-50,A60' -> '[A10-50/A60]'. Bare chain-ID tokens
    // (whole-chain shorthand) are expanded to that chain's full range in the input PDB.
    private String resolveInpaintBracket(String spec) {
        if (!spec) {
            return null
        }
        def pdbFile = new File(this.input_pdb.toString())
        def tokens = spec.split(',').collect { token ->
            def matcher = (token =~ /^([A-Za-z]+)(\d+)?(?:-(\d+))?$/)
            matcher.matches()
            def (chain, start) = [matcher.group(1), matcher.group(2)]
            start ? token : formatChainRanges(pdbFile, chain)
        }
        return "[${tokens.join('/')}]"
    }

    // Translate a chain-qualified residue-spec string (rfd_partialdiff_fixed_seq) into RFdiffusion's
    // flat 0-indexed contigmap.provide_seq ranges. Only valid for the design chain (chain A in
    // binder mode, the sole chain in monomer mode), since that chain is always the first/bare-digit
    // segment of a partial-diffusion contig by ProteinDJ convention - its rank within the chain's own
    // sorted residue numbers therefore equals its flat 0-based position in the assembled structure.
    private String resolveProvideSeq(String fixedSeqSpec) {
        def pdbFile = new File(this.input_pdb.toString())
        def designChain = this.is_binder_mode ? 'A' : Utils.getPdbChainIds(pdbFile).sort().first()
        def flatIndices = [] as SortedSet
        fixedSeqSpec.split(',').each { token ->
            def matcher = (token =~ /^([A-Za-z]+)(\d+)?(?:-(\d+))?$/)
            matcher.matches()
            def (chain, start, end) = [matcher.group(1), matcher.group(2), matcher.group(3)]
            if (chain != designChain) {
                throw new IllegalArgumentException("rfd_partialdiff_fixed_seq token '${token}' must reference the design chain ('${designChain}') only.")
            }
            def chainResNums = Utils.getPdbChainResidueNumbers(pdbFile, chain)
            def targetNums = !start ? chainResNums : (end ? (start.toInteger()..end.toInteger()) : [start.toInteger()])
            targetNums.each { resNum ->
                def rank = chainResNums.indexOf(resNum)
                if (rank == -1) {
                    throw new IllegalArgumentException("rfd_partialdiff_fixed_seq residue ${chain}${resNum} not found in input_pdb chain ${chain}.")
                }
                flatIndices << rank
            }
        }
        def ranges = []
        def rangeStart = null
        def rangeEnd = null
        flatIndices.each { idx ->
            if (rangeStart == null) {
                rangeStart = idx
                rangeEnd = idx
            } else if (idx == rangeEnd + 1) {
                rangeEnd = idx
            } else {
                ranges << "${rangeStart}-${rangeEnd}".toString()
                rangeStart = idx
                rangeEnd = idx
            }
        }
        if (rangeStart != null) {
            ranges << "${rangeStart}-${rangeEnd}".toString()
        }
        return ranges.join(',')
    }

    private void addCommonParameters(List<String> cmd) {
        cmd << "inference.write_trajectory=False"
        cmd << "inference.output_prefix=./rfd_results/fold"

        def variant = this.design_mode - 'rfd_'  // 'denovo' | 'foldcond' | 'motifscaff' | 'partialdiff'
        boolean isBinder = this.is_binder_mode

        // Contigs apply to every variant except foldcond (which uses scaffoldguided.* instead)
        if (variant != 'foldcond' && this.resolvedContigs) {
            cmd << "\'contigmap.contigs=${this.resolvedContigs}\'"
        }

        // input_pdb is always required for motifscaff/partialdiff; only for denovo/foldcond in binder mode
        if ((isBinder || variant in ['motifscaff', 'partialdiff']) && this.input_pdb) {
            cmd << "inference.input_pdb=${getSimpleFileName(this.input_pdb)}"
        }
        
        // Add model parameter with automatic "_ckpt.pt" suffix
        if (this.rfd_ckpt_override) {
            String modelName = MODEL_NAMES[this.rfd_ckpt_override.toLowerCase()]
            if (modelName) {
                cmd << "inference.ckpt_override_path=${MODEL_BASE_PATH}${modelName}_ckpt.pt"
            } else {
                // If not in our mapping, use directly with _ckpt.pt suffix
                cmd << "inference.ckpt_override_path=${MODEL_BASE_PATH}${this.rfd_ckpt_override}_ckpt.pt"
            }
        } else if (variant == 'motifscaff' && isBinder && !this.motifscaff_inpaint_seq && !this.flexible_residues) {
            // rfd_motifscaff never sets ppi.hotspot_res, so without this RFdiffusion's own
            // checkpoint auto-selection would fall back to the monomer-only Base_ckpt.pt even
            // though a target chain is present. Skipped when inpaint_seq/inpaint_str is used
            // since there is no checkpoint trained for complex-conditioning + inpainting together,
            // so RFdiffusion's InpaintSeq_ckpt.pt auto-selection is left in charge instead.
            cmd << "inference.ckpt_override_path=${MODEL_BASE_PATH}Complex_base_ckpt.pt"
        }
        
        // Add noise scale parameters
        if (this.rfd_noise_scale != null) {
            cmd << "denoiser.noise_scale_ca=${this.rfd_noise_scale} denoiser.noise_scale_frame=${this.rfd_noise_scale}"
        }

        // Add any essential extra configurations
        if (this.rfd_extra_config) {
            cmd << "${this.rfd_extra_config}"
        }
    }
    
    private void addBinderDenovoParameters(List<String> cmd) {
        // Add hotspots validation and parameter in rfd_denovo binder mode
        if (this.hotspot_residues) {
            def expandedHotspots = expandHotspots(this.hotspot_residues)
            // Validate hotspots against the auto-generated contigs before adding to command
            if (this.resolvedContigs) {
                validateHotspots(this.resolvedContigs, expandedHotspots)
            }
            cmd << "\'ppi.hotspot_res=[${expandedHotspots}]\'"
        }
        def inpaintStr = resolveInpaintBracket(this.flexible_residues)
        if (inpaintStr) {
            cmd << "contigmap.inpaint_str=${inpaintStr}"
        }
    }
    
    private void addBinderFoldConditioningParameters(List<String> cmd) {
        cmd << "scaffoldguided.scaffoldguided=True"
        cmd << "scaffoldguided.target_pdb=True"

        // scaffoldguided.mask_loops defaults to True. Override if false.
        if (this.rfd_foldcond_mask_loops == false) {
            cmd << "scaffoldguided.mask_loops=False"
        }
        
        // Use simple filenames
        // Note inference.input_pdb is also passed to avoid an error although target_path is the value that is actually used
        cmd << "scaffoldguided.target_path=${getSimpleFileName(this.input_pdb)}"
        cmd << "scaffoldguided.target_ss=target_ss.pt"
        cmd << "scaffoldguided.target_adj=target_adj.pt"
        
        // For scaffolds_dir, use a relative path
        if (this.rfd_foldcond_scaffold_dir) {
            cmd << "scaffoldguided.scaffold_dir=${getSimpleFileName(this.rfd_foldcond_scaffold_dir)}"
        }
        
        // Add hotspots for rfd_foldcond binder mode
        if (this.hotspot_residues) {
            cmd << "\'ppi.hotspot_res=[${expandHotspots(this.hotspot_residues)}]\'"
        }
        def inpaintStr = resolveInpaintBracket(this.flexible_residues)
        if (inpaintStr) {
            cmd << "contigmap.inpaint_str=${inpaintStr}"
        }
    }

    private void addBinderMotifScaffoldingParameters(List<String> cmd) {
        if (this.rfd_motifscaff_length) {
            cmd << "contigmap.length=${this.rfd_motifscaff_length}"  
        } 
        def inpaintSeq = resolveInpaintBracket(this.motifscaff_inpaint_seq)
        if (inpaintSeq) {
            cmd << "contigmap.inpaint_seq=${inpaintSeq}" 
        }
        def inpaintStr = resolveInpaintBracket(this.flexible_residues)
        if (inpaintStr) {
            cmd << "contigmap.inpaint_str=${inpaintStr}"
        }
    }

    private void addBinderPartialDiffusionParameters(List<String> cmd) {
        cmd << "diffuser.partial_T=${this.rfd_partialdiff_timesteps}"
        if (this.rfd_partialdiff_fixed_seq) {
            cmd << "\'contigmap.provide_seq=[${resolveProvideSeq(this.rfd_partialdiff_fixed_seq)}]\'"
        }
        def inpaintStr = resolveInpaintBracket(this.flexible_residues)
        if (inpaintStr) {
            cmd << "contigmap.inpaint_str=${inpaintStr}"
        }
    }

    private void addMonomerDenovoParameters(List<String> cmd) {
        // Add 'placeholder' PDB file, since RFdiffusion requires xyz coordinates
        cmd << "inference.input_pdb=placeholder.pdb"
    }

    private void addMonomerFoldConditioningParameters(List<String> cmd) {
        cmd << "scaffoldguided.scaffoldguided=True"

        // Add 'placeholder' PDB file, since RFdiffusion requires xyz coordinates
        cmd << "inference.input_pdb=placeholder.pdb"

        // scaffoldguided.mask_loops defaults to True. Override if false.
        if (this.rfd_foldcond_mask_loops == false) {
            cmd << "scaffoldguided.mask_loops=False"
        }

        // For scaffolds_dir, use a relative path
        if (this.rfd_foldcond_scaffold_dir) {
            cmd << "scaffoldguided.scaffold_dir=${getSimpleFileName(this.rfd_foldcond_scaffold_dir)}"
        }

    }

    private void addMonomerMotifScaffoldingParameters(List<String> cmd) {
        def inpaintSeq = resolveInpaintBracket(this.motifscaff_inpaint_seq)
        if (inpaintSeq) {
            cmd << "contigmap.inpaint_seq=${inpaintSeq}" 
        }
        def inpaintStr = resolveInpaintBracket(this.flexible_residues)
        if (inpaintStr) {
            cmd << "contigmap.inpaint_str=${inpaintStr}"
        }
        if (this.rfd_motifscaff_length) {
            cmd << "contigmap.length=${this.rfd_motifscaff_length}"  
        } 
    }

    private void addMonomerPartialDiffusionParameters(List<String> cmd) {
        cmd << "diffuser.partial_T=${this.rfd_partialdiff_timesteps}"
        if (this.rfd_partialdiff_fixed_seq) {
            cmd << "\'contigmap.provide_seq=[${resolveProvideSeq(this.rfd_partialdiff_fixed_seq)}]\'"
        }
        def inpaintStr = resolveInpaintBracket(this.flexible_residues)
        if (inpaintStr) {
            cmd << "contigmap.inpaint_str=${inpaintStr}"
        }
    }


    String generateCommandString(String contigsOverride = null) {
        // Contigs are always supplied by main.nf's channel wiring: either a concrete contig
        // string (rfd_denovo monomer/foldcond) or an auto-generation sentinel resolved below.
        def effectiveContigs = contigsOverride

        // Auto-generate rfd_motifscaff contigs from motifscaff_spec/target-chain detection when the
        // contigs channel carries the auto-generation sentinel instead of a real contig string.
        if (effectiveContigs == AUTO_CONTIGS_SENTINEL) {
            effectiveContigs = resolveMotifScaffContigs()
        } else if (effectiveContigs == AUTO_PARTIALDIFF_SENTINEL) {
            effectiveContigs = resolvePartialDiffContigs()
        }

        // Store the resolved contigs for use by addCommonParameters/addBinderDenovoParameters
        this.resolvedContigs = effectiveContigs
        
        def cmd = [RFD_SCRIPT_PATH]
        
        // Add common parameters
        addCommonParameters(cmd)
        
        // Add mode-specific parameters, dispatching on variant + auto-detected monomer/binder status
        def variant = this.design_mode - 'rfd_'
        switch (variant) {
            case 'denovo':
                this.is_binder_mode ? addBinderDenovoParameters(cmd) : addMonomerDenovoParameters(cmd)
                break
            case 'foldcond':
                this.is_binder_mode ? addBinderFoldConditioningParameters(cmd) : addMonomerFoldConditioningParameters(cmd)
                break
            case 'motifscaff':
                this.is_binder_mode ? addBinderMotifScaffoldingParameters(cmd) : addMonomerMotifScaffoldingParameters(cmd)
                break
            case 'partialdiff':
                this.is_binder_mode ? addBinderPartialDiffusionParameters(cmd) : addMonomerPartialDiffusionParameters(cmd)
                break
        }
        return cmd.join(' ')
    }
}

