// Strategy object for the pipeline's 3 fold-design engines (BindCraft, BoltzGen, RFdiffusion).
// Encapsulates engine-specific parameter validation, startup logging, and required input file
// collection, so main.nf can dispatch on design_mode once instead of repeating per-engine checks
// inline. Channel/process wiring (which needs Nextflow's workflow scope) stays in main.nf.
abstract class DesignEngine {
    Map params
    boolean isBinderMode

    DesignEngine(Map params, boolean isBinderMode) {
        this.params = params
        this.isBinderMode = isBinderMode
    }

    // Whether this mode/monomer-binder combination requires params.input_pdb to be set.
    // Defaults to true; override for engines/modes that support monomer design without a target.
    boolean requiresInputPdb() {
        return true
    }

    // Validate engine-specific parameters, throwing IllegalArgumentException/FileNotFoundException on failure
    abstract void validateParams()

    // Human-readable summary lines of the design parameters that will be used (caller does the logging)
    abstract List<String> startupMessages()

    // Paths (as Strings, exactly as supplied in params) of files that must be copied into the
    // run's inputs/ directory. Engine-agnostic extras (e.g. RFdiffusion's placeholder.pdb, which
    // lives under projectDir) are added by the caller since projectDir isn't available here.
    abstract List<String> collectInputFiles()

    static DesignEngine forMode(Map params, boolean isBinderMode) {
        if (params.design_mode == 'bindcraft_denovo') {
            return new BindCraftEngine(params, isBinderMode)
        }
        if (params.design_mode in ['boltzgen_denovo', 'boltzgen_motifscaff']) {
            return new BoltzGenEngine(params, isBinderMode)
        }
        if (params.design_mode == 'rfd4_denovo') {
            return new RFD4Engine(params, isBinderMode)
        }
        return new RFDEngine(params, isBinderMode)
    }
}
