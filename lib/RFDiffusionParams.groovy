import java.nio.file.Paths

public class RFDiffusionParams extends HashMap<String, Object> {
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
    
    private String getSimpleFileName(String filePath) {
        return Paths.get(filePath).getFileName().toString()
    }

    private void addCommonParameters(List<String> cmd) {
        cmd << "inference.write_trajectory=False"
        cmd << "inference.output_prefix=./rfd_results/fold"
        
        // Add contigs if applicable to this mode
        if (this.design_mode in ['binder_denovo', 'binder_partialdiff', 'binder_motifscaff', 'monomer_denovo', 'monomer_motifscaff', 'monomer_partialdiff'] && this.rfd_contigs) {
            cmd << "\'contigmap.contigs=${this.rfd_contigs}\'"
        }
        
        // Use just filename, will be in the process working directory
        if (this.design_mode in ['binder_denovo', 'binder_foldcond', 'binder_partialdiff', 'binder_motifscaff', 'monomer_motifscaff', 'monomer_partialdiff'] && this.input_pdb) {
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
        if (this.hotspot_residues) {
            // Validate hotspots before adding to command
            if (this.rfd_contigs) {
                validateHotspots(this.rfd_contigs, this.hotspot_residues)
            }
            cmd << "\'ppi.hotspot_res=[${this.hotspot_residues}]\'"
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
        cmd << "scaffoldguided.target_path=${getSimpleFileName(this.input_pdb)}"
        cmd << "scaffoldguided.target_ss=target_ss.pt"
        cmd << "scaffoldguided.target_adj=target_adj.pt"
        
        // For scaffolds_dir, use a relative path
        if (this.rfd_scaffold_dir) {
            cmd << "scaffoldguided.scaffold_dir=${getSimpleFileName(this.rfd_scaffold_dir)}"
        }
        
        // Add hotspots for binder_foldcond mode
        if (this.hotspot_residues) {
            cmd << "\'ppi.hotspot_res=[${this.hotspot_residues}]\'"
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


    String generateCommandString(String contigsOverride = null) {
        // Use override if provided, otherwise use params value
        def effectiveContigs = contigsOverride ?: this.rfd_contigs
        
        // Temporarily set contigs for validation
        this.rfd_contigs = effectiveContigs
        
        def cmd = [RFD_SCRIPT_PATH]
        
        // Add common parameters
        addCommonParameters(cmd)
        
        // Add mode-specific parameters
        switch (this.design_mode) {
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

