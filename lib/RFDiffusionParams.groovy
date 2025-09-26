import java.nio.file.Paths

public class RFDiffusionParams extends HashMap<String, Object> {
    static final Set<String> VALID_MODES = ['binder_denovo', 'binder_foldcond', 'binder_motifscaff', 'binder_partialdiff', 'monomer_denovo', 'monomer_foldcond', 'monomer_motifscaff', 'monomer_partialdiff'] as Set
    static final String RFD_SCRIPT_PATH = "/app/RFdiffusion/scripts/run_inference.py"
    
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
                    def chain = match[0][1] ?: 'A'
                    def (start, end) = [match[0][2], match[0][3]]*.toInteger()
                    chainRanges[chain] = chainRanges.getOrDefault(chain, []) << [start, end]
                }
            }
        }
        // Hotspot validation
        hotspots.replaceAll(/[\[\]]/, '').split(/,\s*/).each { hs ->
            def match = (hs =~ /^([A-Z])?(\d+)$/)
            def chain = match[0][1] ?: 'A'
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
    
    // Function to generate human-readable contig descriptions
    def generateContigDescription(String contigs) {
        // Clean and split the contig string into chain parts
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
            // Process each segment in current chain part
            segments.each { seg ->
                switch (seg) {
                    case '0':
                        description << "* Insert a chainbreak after chain ${chainId}"
                        break
                    default:
                        def matcher = seg =~ /^([A-Za-z])?(\d+)-(\d+)$/
                        if (matcher.matches()) {
                            def (segChain, start, end) = [matcher.group(1), matcher.group(2), matcher.group(3)]
                            if (segChain) {
                                description << "* Keep residues ${start}-${end} of chain ${segChain}"
                            } else {
                                if (chainId) {
                                    description << "* Diffuse ${start}-${end} residues for chain ${chainId}"
                                } else {
                                    description << "* Diffuse ${start}-${end} residues for a new chain"
                                }
                            }
                        }
                }
            }
        }
        return description.join('\n')
    }
   
    private void validateParameters() {
        // Validate mode
        if (!(this.rfd_mode in VALID_MODES)) {
            throw new IllegalArgumentException("Invalid mode: ${this.rfd_mode}. Must be one of: ${VALID_MODES.join(', ')}")
        }
       
        // Validate mode-specific parameters
        validateModeSpecificParameters()
    }
   
    private void validateModeSpecificParameters() {
        // Common parameters for multiple modes
        if (this.rfd_mode in ['binder_denovo', 'binder_motifscaff', 'binder_partialdiff', 'monomer_denovo', 'monomer_motifscaff', 'monomer_partialdiff']) {
            validateContigs()
        }
       
        if (this.rfd_mode in ['binder_denovo', 'binder_motifscaff', 'binder_partialdiff', 'binder_foldcond', 'monomer_motifscaff', 'monomer_partialdiff']) {
            validateInputPdb()
        }
       
        // Mode-specific parameters
        switch (this.rfd_mode) {
            case 'binder_partialdiff':
                if (!this.rfd_partial_diffusion_timesteps) {
                    throw new IllegalArgumentException("rfd_partial_diffusion_timesteps is required when mode is 'binder_partialdiff'")
                }
                break            
            case 'binder_foldcond':
                validateBinderFoldConditioningParameters()
                break
            case 'monomer_foldcond':
                validateMonomerFoldConditioningParameters()
                break                
            case 'monomer_partialdiff':
                if (!this.rfd_partial_diffusion_timesteps) {
                    throw new IllegalArgumentException("rfd_partial_diffusion_timesteps is required when mode is 'monomer_partialdiff'")
                }
                break                
        }
    }
   
    private void validateContigs() {
        if (!this.rfd_contigs) {
            throw new IllegalArgumentException("Please provide valid contigs for RFdiffusion")
        }
        // Print descriptions of contigs
        println generateContigDescription(this.rfd_contigs)
        println ""
    }
   
    private void validateInputPdb() {
        if (!this.rfd_input_pdb) {
            throw new IllegalArgumentException("Please provide input PDB file path for RFdiffusion")
        }
       
        def inputFile = new File(this.rfd_input_pdb)
        if (!inputFile.exists()) {
            throw new FileNotFoundException("Input PDB file not found at path: ${this.rfd_input_pdb}. Please ensure the file exists and the path is correct.")
        }
    }
   
    private void validateBinderFoldConditioningParameters() {
        // Validate scaffold directory and contents
        if (!this.rfd_scaffold_dir) {
            throw new IllegalArgumentException("Please provide path to directory containing scaffold files for fold conditioning (rfd_scaffold_dir)")
        }
        def scaffoldsDir = new File(this.rfd_scaffold_dir)
        if (!scaffoldsDir.exists() || !scaffoldsDir.isDirectory()) {
            throw new IllegalArgumentException("rfd_scaffold_dir does not exist or is not a directory")
        }
        
        def ss_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_ss.pt') }
        def adj_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_adj.pt') }
        if (!ss_files || !adj_files) {
            throw new IllegalArgumentException("rfd_scaffold_dir does not contain required _ss.pt and _adj.pt files")
        }
        
        // Validate scaffold target files
        if (this.rfd_target_ss) {
            def ssFile = new File(this.rfd_target_ss)
            if (!ssFile.exists()) {
                throw new FileNotFoundException("*_ss.pt file for scaffold-guided binder design not found")
            }
        } else {
            throw new IllegalArgumentException("rfd_target_ss is required for binder_foldcond mode")
        }
        
        if (this.rfd_target_adj) {
            def adjFile = new File(this.rfd_target_adj)
            if (!adjFile.exists()) {
                throw new FileNotFoundException("*_adj.pt file for scaffold-guided binder design not found")
            }
        } else {
            throw new IllegalArgumentException("rfd_target_adj is required for binder_foldcond mode")
        }
    }

    private void validateMonomerFoldConditioningParameters() {
        if (!this.rfd_scaffold_dir) {
            throw new IllegalArgumentException("Please provide path to directory containing scaffold files for fold conditioning (rfd_scaffold_dir)")
        }

        def scaffoldsDir = new File(this.rfd_scaffold_dir)
        if (!scaffoldsDir.exists() || !scaffoldsDir.isDirectory()) {
            throw new IllegalArgumentException("rfd_scaffold_dir does not exist or is not a directory")
        }
        
        def ss_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_ss.pt') }
        def adj_files = scaffoldsDir.listFiles().findAll { it.name.endsWith('_adj.pt') }
        if (!ss_files || !adj_files) {
            throw new IllegalArgumentException("rfd_scaffold_dir does not contain required _ss.pt and _adj.pt files")
        }

    }
    
    private String getSimpleFileName(String filePath) {
        return Paths.get(filePath).getFileName().toString()
    }

    private void addCommonParameters(List<String> cmd) {
        cmd << "inference.write_trajectory=False"
        cmd << "inference.output_prefix=./rfd_results/fold"
        
        // Add contigs if applicable to this mode
        if (this.rfd_mode in ['binder_denovo', 'binder_partialdiff', 'binder_motifscaff', 'monomer_denovo', 'monomer_motifscaff', 'monomer_partialdiff'] && this.rfd_contigs) {
            cmd << "\'contigmap.contigs=${this.rfd_contigs}\'"
        }
        
        // Use just filename, will be in the process working directory
        if (this.rfd_mode in ['binder_denovo', 'binder_foldcond', 'binder_partialdiff', 'binder_motifscaff', 'monomer_motifscaff', 'monomer_partialdiff'] && this.rfd_input_pdb) {
            cmd << "inference.input_pdb=${getSimpleFileName(this.rfd_input_pdb)}"
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
        // Add hotspots validation and parameter in binder_denovo mode
        if (this.rfd_hotspots) {
            // Validate hotspots before adding to command
            if (this.rfd_contigs) {
                validateHotspots(this.rfd_contigs, this.rfd_hotspots)
            }
            cmd << "\'ppi.hotspot_res=${this.rfd_hotspots}\'"
        }
    }
    
    private void addBinderFoldConditioningParameters(List<String> cmd) {
        cmd << "scaffoldguided.scaffoldguided=True"
        cmd << "scaffoldguided.target_pdb=True"

        // scaffoldguided.mask_loops defaults to True. Override if false.
        if (this.rfd_mask_loops == false) {
            cmd << "scaffoldguided.mask_loops=False"
        }
        
        // Use simple filenames
        // Note inference.input_pdb is also passed to avoid an error although target_path is the value that is actually used
        cmd << "scaffoldguided.target_path=${getSimpleFileName(this.rfd_input_pdb)}"
        cmd << "scaffoldguided.target_ss=${getSimpleFileName(this.rfd_target_ss)}"
        cmd << "scaffoldguided.target_adj=${getSimpleFileName(this.rfd_target_adj)}"
        
        // For scaffolds_dir, use a relative path
        if (this.rfd_scaffold_dir) {
            cmd << "scaffoldguided.scaffold_dir=${getSimpleFileName(this.rfd_scaffold_dir)}"
        }
        
        // Add hotspots for binder_foldcond mode
        if (this.rfd_hotspots) {
            cmd << "\'ppi.hotspot_res=${this.rfd_hotspots}\'"
        }
    }

    private void addBinderMotifScaffoldingParameters(List<String> cmd) {
        if (this.rfd_length) {
            cmd << "contigmap.length=${this.rfd_length}"  
        } 
        if (this.rfd_inpaint_seq) {
            cmd << "contigmap.inpaint_seq=${this.rfd_inpaint_seq}" 
        }
    }

    private void addBinderPartialDiffusionParameters(List<String> cmd) {
        cmd << "diffuser.partial_T=${this.rfd_partial_diffusion_timesteps}"
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
        if (this.rfd_mask_loops == false) {
            cmd << "scaffoldguided.mask_loops=False"
        }

        // For scaffolds_dir, use a relative path
        if (this.rfd_scaffold_dir) {
            cmd << "scaffoldguided.scaffold_dir=${getSimpleFileName(this.rfd_scaffold_dir)}"
        }

    }

    private void addMonomerMotifScaffoldingParameters(List<String> cmd) {
        if (this.rfd_inpaint_seq) {
            cmd << "contigmap.inpaint_seq=${this.rfd_inpaint_seq}" 
        }
        if (this.rfd_length) {
            cmd << "contigmap.length=${this.rfd_length}"  
        } 
    }

    private void addMonomerPartialDiffusionParameters(List<String> cmd) {
        cmd << "diffuser.partial_T=${this.rfd_partial_diffusion_timesteps}"
    }


    String generateCommandString() {
        validateParameters()
        
        def cmd = [RFD_SCRIPT_PATH]
        
        // Add common parameters
        addCommonParameters(cmd)
        
        // Add mode-specific parameters
        switch (this.rfd_mode) {
            case 'binder_denovo':
                addBinderDenovoParameters(cmd)
                break
            case 'binder_foldcond':
                addBinderFoldConditioningParameters(cmd)
                break
            case 'binder_motifscaff':
                addBinderMotifScaffoldingParameters(cmd)
                break
            case 'binder_partialdiff':
                addBinderPartialDiffusionParameters(cmd)
                break
            case 'monomer_denovo':
                addMonomerDenovoParameters(cmd)
                break  
            case 'monomer_foldcond':
                addMonomerFoldConditioningParameters(cmd)
                break  
            case 'monomer_motifscaff':
                addMonomerMotifScaffoldingParameters(cmd)
                break  
            case 'monomer_partialdiff':
                addMonomerPartialDiffusionParameters(cmd)
                break                
        }
        return cmd.join(' ')
    }
}

