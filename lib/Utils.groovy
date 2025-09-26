class Utils {
    // Method to rebatch channels of tuples of PDBs and JSON files
    static def rebatchTuples(input_channel, batch_size = 50) {
        return input_channel 
            .transpose()
            .buffer(size: batch_size, remainder: true)
            .map { pairs -> 
                def first_elements = pairs.collect { it[0] }
                def second_elements = pairs.collect { it[1] }
                return [first_elements, second_elements]
            }
    }

    // Method to rebatch channels of PDBs by number of GPUs
    static def rebatchGPU(input_channel, gpus) {
        return input_channel
            .collect()
            .flatMap { all_pdbs ->
                def total_size = all_pdbs.size()
                def nbatches = Math.min(gpus, total_size)
                def bsize = (total_size / nbatches).doubleValue()
                def idx = 0
                all_pdbs.collect { pdb ->
                    def batch_id = Math.floor(idx++ / bsize).intValue()
                    [batch_id, pdb]
                }
            }
            .groupTuple()
    }
    
    // Method to rebatch channels of PDBs by number of GPUs, sorted by residue count
    static def rebatchGPUByNumRes(input_channel, gpus) {
        return input_channel
            .collect()
            .flatMap { all_pdbs ->
                // Sort PDB files by number of residues
                def sorted_pdbs = all_pdbs.sort { pdb ->
                    return countResidues(pdb)
                }
                
                def total_size = sorted_pdbs.size()
                def nbatches = Math.min(gpus, total_size)
                def bsize = (total_size / nbatches).doubleValue()
                def idx = 0
                sorted_pdbs.collect { pdb ->
                    def batch_id = Math.floor(idx++ / bsize).intValue()
                    [batch_id, pdb]
                }
            }
            .groupTuple()
    }

    // Helper function to count unique residues in a PDB file
    static def countResidues(pdb_file) {
        def residue_set = new HashSet()
        
        pdb_file.eachLine { line ->
            // Only process ATOM and HETATM lines
            if (line.startsWith("ATOM  ") || line.startsWith("HETATM")) {
                // Extract chain ID (column 22) and residue number (columns 23-26)
                if (line.length() >= 26) {
                    def chain_id = line.substring(21, 22)
                    def res_num = line.substring(22, 26).trim()
                    def res_name = line.substring(17, 20).trim()
                    
                    // Create unique identifier: chain + residue_number + residue_name
                    def residue_key = "${chain_id}_${res_num}_${res_name}"
                    residue_set.add(residue_key)
                }
            }
        }
        
        return residue_set.size()
    }
    
    /**
     * Format filter parameters for command-line arguments
     * @param params The Nextflow params object
     * @param paramPrefix The prefix for parameter names
     * @param paramNames List of parameter names to process
     * @return Formatted parameter string for command line
     */
    static def formatFilterParams(params, paramPrefix, paramNames) {
        return paramNames.collect { name ->
            def paramValue = params["${paramPrefix}_${name}"]
            if (paramValue != null) {
                def cmdParam = name.replaceAll('_', '-')
                return "--${paramPrefix}-${cmdParam} ${paramValue}"
            } else {
                return ""
            }
        }.findAll { it != "" }.join(' ')
    }

    /**
     * Count PDB files from a channel
     * @param channel The input channel containing files
     * @return A value channel with the count of PDB files
     */
    static def countPdbFiles(channel) {
        return channel
            .flatten()
            .collect()
            .map { files -> 
                files.findAll { file -> 
                    file.toString().endsWith('.pdb') 
                }.size() 
            }
            .ifEmpty(0)
    }
}
