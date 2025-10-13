[🏠 ProteinDJ](../README.md) > Parameter Guide

# ProteinDJ Parameter Guide

This guide summarizes the key parameters used in ProteinDJ's Nextflow pipeline configuration. Parameters are essential for controlling the design mode, input/output locations, model choices, advanced options, filtering thresholds, and cluster resource allocation.

---

## Essential Parameters

These parameters are required for ProteinDJ and are used by every mode.

| Parameter         | Default         | Description                                                                                                                                                                                     |
| ----------------- | --------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `design_mode`        | null            | Pipeline mode. Choose from: `monomer_denovo`, `monomer_foldcond`, `monomer_motifscaff`, `monomer_partialdiff`, `binder_denovo`, `binder_foldcond`, `binder_motifscaff`, `binder_partialdiff`, or `bindcraft` |
| `num_designs` | 8               | Number of designs to generate using RFdiffusion or Bindcraft                                                                                                                                                 |
| `seqs_per_design` | 8               | Number of sequences to generate per design                                                                                                                                          |
| `out_dir`         | `./pdj_results` | Output directory for results. Existing results will be overwritten                                                                                                                              |

---

## Mode-Specific Parameters

These parameters are used for some of the ProteinDJ modes.

| Parameter                         | Default | Description                                                                                                                                                        | Required for Modes                                                                                                         |
| --------------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------- |
| `rfd_contigs`                     | null    | Contigs specification strings, e.g. `[A17-145/0 50-100]`. See docs for examples. Required for multiple modes.                                                      | `binder_denovo`, `binder_motifscaff`, `binder_partialdiff`, `monomer_denovo`, `monomer_motifscaff`, `monomer_partialdiff`  |
| `input_pdb`                   | null    | Path to input PDB file (e.g., target for binders, `'./target.pdb'`). Required for several modes.                                                                   | `binder_denovo`, `binder_foldcond`, `binder_motifscaff`, `binder_partialdiff`, `monomer_motifscaff`, `monomer_partialdiff` |
| `hotspot_residues`                    | null    | Hotspot residues for binder design, e.g. `[A56,A115,A123]`. Optional for `binder_denovo` and `binder_foldcond`.                                                    |                                                                                                                            |
| `rfd_scaffold_dir`                | null    | Directory containing scaffold secondary structure and block adjacency files (e.g. `'./binderscaffolds/scaffolds_assorted'`). Required for fold conditioning modes. | `binder_foldcond`, `monomer_foldcond`                                                                                      |
| `rfd_target_ss`                   | null    | Fold conditioning secondary structure file for target (e.g. `'./target_ss.pt'`). Required for `binder_foldcond`.                                                   | `binder_foldcond`                                                                                                          |
| `rfd_target_adj`                  | null    | Fold conditioning block adjacency file for target (e.g. `'./target_adj.pt'`). Required for `binder_foldcond`.                                                      | `binder_foldcond`                                                                                                          |
| `rfd_mask_loops`                  | true    | Whether to ignore loops during scaffold secondary structure constraints. Optional for `binder_foldcond` and `monomer_foldcond` modes.                              |                                                                                                                            |
| `rfd_inpaint_seq`                 | null    | Residues with masked sequence during diffusion, e.g. `[A17-19/A143-145]`. Optional for motif scaffolding modes `monomer_motifscaff`, `monomer_foldcond` .          |                                                                                                                            |
| `rfd_length`                      | null    | Length of diffused chain during motif scaffolding, e.g. `100-100` or `100-120`. Optional for motif scaffolding modes `binder_motifscaff`, `monomer_motifscaff`.    |                                                                                                                            |
| `rfd_partial_diffusion_timesteps` | null    | Number of timesteps for partial diffusion (1-50) e.g. 20. Required for partial diffusion modes.                                                                    | `binder_partialdiff`, `monomer_partialdiff`                                                                                |

### Parameter Requirements by Mode

`-` = ignored

| Parameter                       | monomer denovo | monomer foldcond | monomer motifscaff | monomer partialdiff | binder denovo | binder foldcond | binder motifscaff | binder partialdiff |
| ------------------------------- | -------------- | ---------------- | ------------------ | ------------------- | ------------- | --------------- | ----------------- | ------------------ |
| rfd_contigs                     | Required       | -                | Required           | Required            | Required      | -               | Required          | Required           |
| input_pdb                   | -              | -                | Required           | Required            | Required      | Required        | Required          | Required           |
| hotspot_residues                    | -              | -                | -                  | -                   | _Optional_    | _Optional_      | -                 | -                  |
| rfd_scaffold_dir                | -              | Required         | -                  | -                   | -             | Required        | -                 | -                  |
| rfd_target_adj                  | -              | -                | -                  | -                   | -             | Required        | -                 | -                  |
| rfd_target_ss                   | -              | -                | -                  | -                   | -             | Required        | -                 | -                  |
| rfd_mask_loops                  | -              | _Optional_       | -                  | -                   | -             | _Optional_      | -                 | -                  |
| rfd_inpaint_seq                 | -              | -                | _Optional_         | -                   | -             | -               | _Optional_        | -                  |
| rfd_length                      | -              | -                | _Optional_         | -                   | -             | -               | _Optional_        | -                  |
| rfd_partial_diffusion_timesteps | -              | -                | -                  | Required            | -             | -               | -                 | Required           |
| rfd_ckpt_override               | _Optional_     | _Optional_       | _Optional_         | _Optional_          | _Optional_    | _Optional_      | _Optional_        | _Optional_         |
| rfd_noise_scale                 | _Optional_     | _Optional_       | _Optional_         | _Optional_          | _Optional_    | _Optional_      | _Optional_        | _Optional_         |

---

## Workflow Advanced Parameters

These parameters control the workflow of ProteinDJ.

| Parameter           | Default | Description                                                                                                                                                                                                              |
| ------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `seq_method`        | 'mpnn'  | Sequence design method. Options: `'mpnn'` (ProteinMPNN Fast-Relax) or `'fampnn'` (Full-Atom MPNN).                                                                                                                       |
| `pred_method`       | 'af2'   | Structure prediction method. Options: `'af2'` (AlphaFold2 Initial-Guess) or `'boltz'` (Boltz-2).                                                                                                                         |
| `zip_pdbs`      | true   | Whether to compress output final designs in a tar.gz archive. If false, results will be output as uncompressed PDB files. |
| `run_rfd_only`      | false   | Whether to run only RFdiffusion and skip sequence design, prediction, and analysis.                                                                                                                                      |
| `skip_rfd`          | false   | Skip RFdiffusion, run sequence design, prediction, and analysis only. Requires valid `skip_input_dir` containing PDB and JSON files with metadata. Binder design PDBs must have binder as chain A and target as chain B. |
| `skip_rfd_seq`      | false   | Skip RFdiffusion and sequence design, run structure prediction and analysis only. Requires valid `skip_input_dir` containing PDB files. Binder design PDBs must have binder as chain A and target as chain B.            |
| `skip_rfd_seq_pred` | false   | Skip RFdiffusion, sequence design, and prediction, run analysis only. Requires valid `skip_input_dir` containing PDB files. Binder design PDBs must have binder as chain A and target as chain B.                        |
| `skip_input_dir`    | null    | Directory path for files when skipping stages (e.g. `'./rfd_results'`).                                                                                                                                                  |

---

## RFdiffusion Advanced Parameters

Advanced parameters to control the behaviour of RFdiffusion. Can be used with any mode.

| Parameter           | Default | Description                                                                                                                                                                                                                                                                         |
| ------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `rfd_ckpt_override` | null    | Overrides the diffusion model checkpoint. Options include `'active_site'`, `'base'`, `'base_epoch8'`, `'complex_base'`, `'complex_beta'`, `'complex_fold_base'`, `'inpaint_seq'`, `'inpaint_seq_fold'`. Not generally recommended except for `binder_denovo` with `'complex_beta'`. |
| `rfd_noise_scale`   | null    | Scale of noise applied to translations and rotations. Default is 1 for monomer modes (increases diversity), recommended 0-0.5 for binders (increases success rates).                                                                                                                |
| `rfd_extra_config`  | null    | Additional raw configuration passed to RFdiffusion not covered by Nextflow parameters e.g. `'contigmap.inpaint_str=[B165-178]'`.                                                                                                                                                    |

---

## ProteinMPNN-FastRelax Advanced Parameters

Advanced parameters to control the behaviour of ProteinMPNN-FastRelax.

| Parameter                           | Default      | Description                                                                                                                                                                             |
| ----------------------------------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mpnn_omitAAs`                      | `'CX'`       | Residue types to exclude during design (one-letter code, case insensitive).                                                                                                             |
| `mpnn_temperature`                  | 0.1          | Temperature for sequence sampling; higher values increase diversity. Recommended to lower significantly for binders to improve success (e.g., 0.0001).                                  |
| `mpnn_checkpoint_type`              | `'soluble'`  | Checkpoint selection: `'vanilla'` or `'soluble'`.                                                                                                                                       |
| `mpnn_checkpoint_model`             | `'v_48_020'` | Checkpoint model variant indicating backbone noise level used during training. e.g. 'v_48_020' noised with 0.20Å Gaussian noise. Choose from 'v_48_002','v_48_010,'v_48_020','v_48_030' |
| `mpnn_backbone_noise`               | 0            | Std dev Gaussian noise added to backbone atoms. 0 = none, 0.1-0.3 = mild perturbation.                                                                                                  |
| `mpnn_relax_max_cycles`             | 0            | Max fast relaxation cycles per sequence; 0 disables FastRelax functionality.                                                                                                            |
| `mpnn_relax_seqs_per_cycle`         | 1            | Number of sequences generated between relaxation steps; best scored sequence is kept.                                                                                                   |
| `mpnn_relax_output`                 | false        | Whether to run relaxation before saving output.                                                                                                                                         |
| `mpnn_relax_convergence_rmsd`       | 0.2          | RMSD early convergence threshold for relaxation cycles. Design is considered converged if the C-alpha RMSD (Å) between cycles is <= this threshold.                                     |
| `mpnn_relax_convergence_score`      | 0.1          | Score early convergence threshold for relaxation cycles. Design is considered converged if the improvement in score between cycles is <= this threshold.                                |
| `mpnn_relax_convergence_max_cycles` | 1            | Design is considered converged if it meets both convergence criteria for n consecutive cycles.                                                                                          |
| `mpnn_extra_config`                 | null         | Additional raw configuration passed to ProteinMPNN-FastRelax. e.g. `-protein_features="full"'`                                                                                          |

---

## Full-Atom MPNN (FAMPNN) Advanced Parameters

Advanced parameters to control the behaviour of Full-Atom MPNN

| Parameter               | Default | Description                                                                                                  |
| ----------------------- | ------- | ------------------------------------------------------------------------------------------------------------ |
| `fampnn__fix_target_sidechains` | false    | Fix target side-chain positions when performing binder sequence design (default: false) |
| `fampnn_psce_threshold` | 0.3     | Will only keep sidechains below this PSCE threshold during design. Null means no filtering.                  |
| `fampnn_temperature`    | 0.1     | Temperature for sampling; higher increases sequence diversity. Recommended lower for binders (e.g., 0.0001). |
| `fampnn_exclude_cys`    | true    | Exclude cysteine residues from design.                                                                       |
| `fampnn_extra_config`   | null    | Additional raw configuration passed to Full-Atom MPNN. e.g. 'timestep_schedule.num_steps=100 seq_only=true'. |

---

## Prediction Advanced Parameters

| Parameter                 | Default | Description                                                                               |
| ------------------------- | ------- | ----------------------------------------------------------------------------------------- |
| `uncropped_target_pdb`    | null    | Path to uncropped target PDB for prediction (e.g., full complex).                         |
| `af2_initial_guess`       | true    | Use an initial guess for target chains in AlphaFold2 structure predictions.               |
| `af2_extra_config`        | null    | Additional raw parameters passed to AlphaFold2 Initial-Guess. e.g. '-recycle 3'           |
| `boltz_recycling_steps`   | 3       | Number of recycling steps in Boltz-2 predictions.                                         |
| `boltz_diffusion_samples` | 1       | Number of diffusion samples in Boltz-2 predictions.                                       |
| `boltz_sampling_steps`    | 200     | Number of sampling steps in Boltz-2 predictions.                                          |
| `boltz_use_potentials`    | false   | Use physics-based potentials during inference to improve physical plausibility of predictions (also known as Boltz-2x). Increases prediction time.                                          |
| `boltz_extra_config`      | null    | Additional raw parameters for Boltz-2 predictions. e.g. '--msa_pairing_strategy complete' |

---

## Filtering Parameters <a name="params-filter"></a>

Due to the inherently stochastic nature of protein design, often we see problematic results during the pipeline. It can save computation time to discard designs mid-pipeline that fail to meet success criteria. We have implemented three filtering stages that can be used to reject poor designs:

- Fold Filtering - Filters designs according to the number of secondary structure elements and radius of gyration.
- Sequence Filtering - Filters designs according to the score of the generated sequence
- AlphaFold2/Boltz-2 Filtering - Filters designs according to the quality of the structure prediction

The most powerful predictors of experimental success are structure prediction metrics, but some metrics are more effective than others. Here are some recommended filters for binder design from the literature and their corresponding parameters in ProteinDJ:

| Parameter                  | RFdiffusion paper<sup>1</sup> | AlphaProteo whitepaper<sup>2</sup> |
| -------------------------- | --------------------- | ---------------------- |
| af2_max_pae_interaction    | 10                    | 7                      |
| af2_min_plddt_overall      | 80                    | 90                     |
| af2_max_rmsd_binder_bndaln | 1                     | 1.5                    |

<sup> 1. Watson, J.L. et al. Nature 620, 1089–1100 (2023). https://doi.org/10.1038/s41586-023-06415-8; 2. Zambaldi, V. et al. arXiv (2024). https://doi.org/10.48550/arXiv.2409.08022
</sup>

We recommend disabling other filters for small-scale and pilot experiments, and using these results to decide on values to use for filtering large-scale runs.

#### RFdiffusion Filtering Parameters

RFdiffusion Filtering Parameters. Metrics are calculated on the binder chain only in binder design modes, otherwise all chains are used in calculations.

| Parameter         | Description                                                             |
| ----------------- | ----------------------------------------------------------------------- |
| `rfd_min_helices` | Minimum number of alpha-helices required.                               |
| `rfd_max_helices` | Maximum number of alpha-helices allowed.                                |
| `rfd_min_strands` | Minimum number of beta-strands required.                                |
| `rfd_max_strands` | Maximum number of beta-strands allowed.                                 |
| `rfd_min_ss`      | Minimum number of secondary structure elements (α-helices + β-strands). |
| `rfd_max_ss`      | Maximum number of secondary structure elements (α-helices + β-strands). |
| `rfd_min_rog`     | Minimum radius of gyration (Å).                                         |
| `rfd_max_rog`     | Maximum radius of gyration (Å).                                         |

#### BindCraft Filtering Parameters

There are no filtering paramters for BindCraft as it has built-in filtering of designs and will automatically reject designs that meet any of the following criteria:
- Low confidence (pLDDT < 0.7)
- Severe clashes (clashes detected between C-alpha atoms)
- Insufficient contact between binder and target (less than three residues contacting the target)

Note that if a design fails, BindCraft will generate a new design until it finds one that meets the criteria. This can lead to long run times compared to RFdiffusion but results in binder designs that are more likely to succeed in the subsequent Structure Prediction stage. 

#### Sequence Filtering Parameters

Sequence Filtering Parameters for ProteinMPNN and Full-Atom MPNN. Recommended to disable unless you know what you are doing.

| Parameter         | Suggested Value | Description                                                                                                      |
| ----------------- | --------------- | ---------------------------------------------------------------------------------------------------------------- |
| `mpnn_max_score`  | `null`          | Maximum MPNN score (negative log likelihood). A lower score means ProteinMPNN is more confident in the sequence. |
| `fampnn_max_psce` | `null`          | Max PSCE score for designed side-chains. A lower score means FAMPNN is more confident in the sequence.           |

#### AlphaFold2 Filtering Parameters

AlphaFold2 Initial-Guess Filtering Parameters.

| Parameter                    | Suggested Value | Description                                            |
| ---------------------------- | --------------- | ------------------------------------------------------ |
| `af2_max_pae_interaction`    | `10`            | Max predicted aligned error for interactions           |
| `af2_max_pae_overall`        | `5`             | Max predicted aligned error for all chains             |
| `af2_max_pae_binder`         | `5`             | Max predicted aligned error for binder                 |
| `af2_max_pae_target`         | `5`             | Max predicted aligned error for target                 |
| `af2_max_rmsd_overall`       | `2`             | Max C-alpha RMSD between AF2 prediction and input design when all chains are aligned           |
| `af2_max_rmsd_binder_bndaln` | `1`             | Max binder C-alpha RMSD between AF2 prediction and input design when binder chains are aligned |
| `af2_max_rmsd_binder_tgtaln` | `2`             | Max binder C-alpha RMSD between AF2 prediction and input design when target chains are aligned |
| `af2_max_rmsd_target`        | `1`             | Max target C-alpha RMSD between AF2 prediction and input design when target chains are aligned |
| `af2_min_plddt_overall`      | `90`            | Min average pLDDT score for all chains                 |
| `af2_min_plddt_binder`       | `85`            | Min pLDDT score required for binder                    |
| `af2_min_plddt_target`       | `90`            | Min pLDDT score required for target                    |

#### Boltz-2 Filtering Parameters

Boltz-2 Filtering Parameters.

| Parameter                   | Suggested Value | Description                                                                                     |
| --------------------------- | --------------- | ----------------------------------------------------------------------------------------------- |
| `boltz_max_overall_rmsd`    | `null`          | Max C-alpha RMSD between all chains of Boltz-2 prediction and input design                        |
| `boltz_max_binder_rmsd`     | `null`          | Max C-alpha RMSD between binder chains of Boltz-2 prediction and input design. Binder modes only. |
| `boltz_max_target_rmsd`     | `null`          | Max C-alpha RMSD between target chains of Boltz-2 prediction and input design. Binder modes only. |
| `boltz_min_conf_score`      | `null`          | Minimum confidence score of the prediction                                                      |
| `boltz_min_ptm`             | `null`          | Minimum predicted template modelling score of the prediction                                    |
| `boltz_min_ptm_interface`   | `null`          | Minimum predicted template modelling score of the prediction interface                          |
| `boltz_min_plddt`           | `null`          | Minimum pLDDT score of the prediction                                                           |
| `boltz_min_plddt_interface` | `null`          | Minimum pLDDT score of the prediction interface                                                 |
| `boltz_max_pde`             | `null`          | Maximum predicted distance error of the prediction                                              |
| `boltz_max_pde_interface`   | `null`          | Maximum predicted distance error for the prediction interface                                   |

## Cluster Parameters

The cluster parameters may need adjusting depending on your HPC setup and available hardware. You must ensure that the paths to the containers and models for RFdiffusion, AlphaFold2 and Boltz-2 are valid and contain the expected files (see ['Installation Guide'](installation.md)).

| Parameter       | Example Value                               | Description                                        |
| --------------- | ------------------------------------------- | -------------------------------------------------- |
| `container_dir` | `'./containers'`                            | Path to pipeline containers directory or cloud URI |
| `rfd_models`    | `"${projectDir}/models/rfd"`                | Path to the RFdiffusion model checkpoints.         |
| `af2_models `   | `"${projectDir}/models/af2"`                | Path to the AlphaFold2 models.                     |
| `boltz_models`  | `"${projectDir}/models/boltz"`              | Path to the Boltz-2 models.                        |
| `gpu_model`     | `'A30'`                                     | GPU model to request, e.g., 'A30'.                 |
| `gpus`          | Number of GPUs to request                   | `1`, `2`, `4`, `8`                                 |
| `cpus_per_gpu`  | Number of CPUs to request per GPU           | `8`, `12`                                          |
| `memory_gpu`    | Memory to request for GPU jobs              | `'24GB'`, `'48GB'`                                 |
| `cpus`          | Number of CPUs to request for CPU-only jobs | `12`, `24`                                         |
| `memory_cpu`    | Memory for request for CPU-only jobs        | `'24GB'`, `'32GB'`                                 |

---

**🌟 Notes:**

- Parameters set to `null` indicate optional or user-defined inputs.
- Ensure paths are updated to your environment when running ProteinDJ.
- Filtering parameters can be set to `null` to disable filtering.

---

[⬅️ Back to Main README](../README.md)
