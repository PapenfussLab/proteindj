class RFD4Engine extends DesignEngine {
    RFD4Engine(Map params, boolean isBinderMode) {
        super(params, isBinderMode)
    }

    void validateParams() {
        Utils.validateDesignLength(params.design_length)
        if (!params.input_pdb) {
            throw new IllegalArgumentException('rfd4_denovo requires input_pdb (the target protein)')
        }
        if (params.hotspot_residues && !params.hotspot_residues.matches(Utils.RESIDUE_SPEC_REGEX)) {
            throw new IllegalArgumentException("hotspot_residues format invalid. Acceptable: 'A56,A115-120,B'.")
        }
        if (params.flexible_residues && !params.flexible_residues.matches(Utils.RESIDUE_SPEC_REGEX)) {
            throw new IllegalArgumentException("flexible_residues format invalid. Acceptable: 'A10-13,A16,B'.")
        }
        if (params.seqs_per_design != 1) {
            throw new IllegalArgumentException('rfd4_denovo generates one sequence per structure; set seqs_per_design=1')
        }
        def inferenceBatch = params.rfd4_inference_batch_size?.toString()
        if (!inferenceBatch?.isInteger() || inferenceBatch.toInteger() < 1) {
            throw new IllegalArgumentException('rfd4_inference_batch_size must be a positive integer')
        }
        if (!params.rfd4_models) {
            throw new IllegalArgumentException('rfd4_models is required')
        }
        def models = new File(params.rfd4_models.toString())
        if (!models.isAbsolute() || !models.isDirectory()) {
            throw new IllegalArgumentException('rfd4_models must be an existing absolute directory')
        }
        if (!new File(models, 'config.json').isFile() || !new File(models, 'model.safetensors').isFile()) {
            throw new IllegalArgumentException('rfd4_models must contain config.json and model.safetensors')
        }
    }

    List<String> startupMessages() {
        def messages = ["Using RFD4-Proteina to co-design binder sequence and structure; refolding with ${params.pred_method}."]
        if (params.flexible_residues) {
            messages << "* Flexible target residues = ${params.flexible_residues} (structure unconditioned; sequence fixed)"
        }
        return messages
    }

    List<String> collectInputFiles() { return [params.input_pdb] }
}
