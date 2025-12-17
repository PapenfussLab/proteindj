[🏠 ProteinDJ](../README.md) > Parameter Guide

# ProteinDJ Parameter Guide

This guide summarizes the key parameters used in ProteinDJ's Nextflow pipeline configuration. Parameters are essential for controlling the design mode, input/output locations, model choices, advanced options, filtering thresholds, and cluster resource allocation.

---

## Essential Parameters

These parameters are required for ProteinDJ and are used by every mode.

| Parameter         | Default         | Description                                                                                                                                                                                                         |
| ----------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `design_mode`     | null            | Pipeline mode. Choose from: `monomer_denovo`, `monomer_foldcond`, `monomer_motifscaff`, `monomer_partialdiff`, `bindcraft_denovo`, `binder_denovo`, `binder_foldcond`, `binder_motifscaff`, or `binder_partialdiff` |
| `num_designs`     | 8               | Number of designs to generate using RFdiffusion or Bindcraft                                                                                                                                                        |
| `seqs_per_design` | 8               | Number of sequences to generate per design                                                                                                                                                                          |
| `out_dir`         | `./pdj_results` | Output directory for results. Existing results will be overwritten                                                                                                                                                  |

---

## Mode-Specific Parameters

These parameters are used for some of the ProteinDJ modes.

| Parameter                         | Default | Description                                                                                                                                                                                                                                                                                                                       | Required for Modes                                                                                                         |
| --------------------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `design_length`                   | null    | Design length of binder or monomer. e.g. `'60'` or `'60-70'`. Pipeline will randomly sample different sizes between these values.                                                                                                                                                                                                 | `bindcraft_denovo`, `binder_denovo`, `monomer_denovo`                                                                      |
| `input_pdb`                       | null    | Path to input PDB file (e.g. target for binders, `'./target.pdb'`). Required for several modes.                                                                                                                                                                                                                                   | `binder_denovo`, `binder_foldcond`, `binder_motifscaff`, `binder_partialdiff`, `monomer_motifscaff`, `monomer_partialdiff` |
| `hotspot_residues`                | null    | Hotspot residues for binder design, e.g. `A56,A115,A123`. Optional for `binder_denovo` and `binder_foldcond`.                                                                                                                                                                                                                     |                                                                                                                            |
| `rfd_contigs`                     | null    | Contigs specification strings for residues to include from input PDB and/or design length, e.g. `[A17-145/0 50-100]`. See docs for examples. Optional for multiple modes, required by `binder_motifscaff` and `monomer_motifscaff`. If null, contigs will be auto generated from input PDB (all residues) and/or `design_length`. | `binder_motifscaff`, `monomer_motifscaff`                                                                                  |
| `rfd_scaffold_dir`                | null    | Directory containing scaffold secondary structure and block adjacency files (e.g. `'./binderscaffolds/scaffolds_assorted'`). Required for fold conditioning modes.                                                                                                                                                                | `binder_foldcond`, `monomer_foldcond`                                                                                      |
| `rfd_mask_loops`                  | true    | Whether to ignore loops during scaffold secondary structure constraints. Optional for `binder_foldcond` and `monomer_foldcond` modes.                                                                                                                                                                                             |                                                                                                                            |
| `rfd_inpaint_seq`                 | null    | Residues with masked sequence during diffusion, e.g. `[A17-19/A143-145]`. Optional for motif scaffolding modes `monomer_motifscaff`, `monomer_foldcond` .                                                                                                                                                                         |                                                                                                                            |
| `rfd_length`                      | null    | Length of diffused chain during motif scaffolding, e.g. `100-100` or `100-120`. Optional for motif scaffolding modes `binder_motifscaff`, `monomer_motifscaff`.                                                                                                                                                                   |                                                                                                                            |
| `rfd_partial_diffusion_timesteps` | null    | Number of timesteps for partial diffusion (1-50) e.g. 20. Required for partial diffusion modes.                                                                                                                                                                                                                                   | `binder_partialdiff`, `monomer_partialdiff`                                                                                |

### Parameter Requirements by Mode

`-` = ignored

| Parameter                       | monomer denovo | monomer foldcond | monomer motifscaff | monomer partialdiff | binder denovo | binder foldcond | binder motifscaff | binder partialdiff | bindcraft_denovo |
| ------------------------------- | -------------- | ---------------- | ------------------ | ------------------- | ------------- | --------------- | ----------------- | ------------------ | ---------------- |
| design_length                   | Required       | -                | -                  | -                   | Required      | -               | -                 | -                  | Required         |
| input_pdb                       | -              | -                | Required           | Required            | Required      | Required        | Required          | Required           | Required         |
| hotspot_residues                | -              | -                | -                  | -                   | _Optional_    | _Optional_      | -                 | -                  | _Optional_       |
| rfd_contigs                     | _Optional_     | -                | Required           | _Optional_          | _Optional_    | -               | Required          | _Optional_         | -                |
| rfd_scaffold_dir                | -              | Required         | -                  | -                   | -             | Required        | -                 | -                  | -                |
| rfd_mask_loops                  | -              | _Optional_       | -                  | -                   | -             | _Optional_      | -                 | -                  | -                |
| rfd_inpaint_seq                 | -              | -                | _Optional_         | -                   | -             | -               | _Optional_        | -                  | -                |
| rfd_length                      | -              | -                | _Optional_         | -                   | -             | -               | _Optional_        | -                  | -                |
| rfd_partial_diffusion_timesteps | -              | -                | -                  | Required            | -             | -               | -                 | Required           | -                |

---

## Workflow Advanced Parameters

These parameters control the workflow of ProteinDJ.

| Parameter            | Default | Description                                                                                                                                                                                                              |
| -------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `seq_method`         | 'mpnn'  | Sequence design method. Options: `'mpnn'` (ProteinMPNN Fast-Relax) or `'fampnn'` (Full-Atom MPNN).                                                                                                                       |
| `pred_method`        | 'af2'   | Structure prediction method. Options: `'af2'` (AlphaFold2 Initial-Guess) or `'boltz'` (Boltz-2).                                                                                                                         |
| `zip_pdbs`           | true    | Whether to compress output final designs in a tar.gz archive. If false, results will be output as uncompressed PDB files.                                                                                                |
| `rank_designs`       | true   | Whether to rank final designs by prediction quality metrics and output `ranked_designs.csv` and `ranked_designs` folder (or `.tar.gz`). |
| `ranking_metric`     | null   | Structure prediction metric to use for ranking designs (e.g., `'af2_pae_interaction'`, `'boltz_ptm_interface'`, `'boltz_ipSAE_min'`). Must match prediction method prefix (`af2_*` or `boltz_*`). If null, defaults are set depending on design_mode and pred_method: binder modes use `'af2_pae_interaction'` (AlphaFold2) or `'boltz_ipSAE_min'` (Boltz-2); monomer modes use `'af2_plddt_overall'` (AlphaFold2) or `'boltz_ptm'` (Boltz-2). See [Metrics Guide](metrics.md) for all available metrics. |
| `max_designs`   | null    | Maximum number of top-ranked designs to output (e.g. 100). If null, all designs are ranked and output.                                                                    |
| `max_seqs_per_fold` | null | Maximum number of sequences per fold to keep after ranking (e.g. 3). Helps increase fold diversity by limiting sequences for highly successful folds. If null, no limit is applied. |
| `run_fold_only`      | false   | Whether to run only fold design and skip sequence design, prediction, and analysis.                                                                                                                                      |
| `skip_fold`          | false   | Skip fold design, run sequence design, prediction, and analysis only. Requires valid `skip_input_dir` containing PDB and JSON files with metadata. Binder design PDBs must have binder as chain A and target as chain B. |
| `skip_fold_seq`      | false   | Skip fold design and sequence design, run structure prediction and analysis only. Requires valid `skip_input_dir` containing PDB files. Binder design PDBs must have binder as chain A and target as chain B.            |
| `skip_fold_seq_pred` | false   | Skip fold design, sequence design, and prediction, run analysis only. Requires valid `skip_input_dir` containing PDB files. Binder design PDBs must have binder as chain A and target as chain B.                        |
| `skip_input_dir`     | null    | Directory path for files when skipping stages (e.g. `'./rfd_results'`).                                                                                                                                                  |

---

## RFdiffusion Advanced Parameters

Advanced parameters to control the behaviour of RFdiffusion. Can be used with any mode.

| Parameter           | Default | Description                                                                                                                                                                                                                                                                         |
| ------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `rfd_ckpt_override` | null    | Overrides the diffusion model checkpoint. Options include `'active_site'`, `'base'`, `'base_epoch8'`, `'complex_base'`, `'complex_beta'`, `'complex_fold_base'`, `'inpaint_seq'`, `'inpaint_seq_fold'`. Not generally recommended except for `binder_denovo` with `'complex_beta'`. |
| `rfd_noise_scale`   | null    | Scale of noise applied to translations and rotations. Default is 1 for monomer modes (increases diversity), recommended 0-0.5 for binders (increases success rates).                                                                                                                |
| `rfd_extra_config`  | null    | Additional raw configuration passed to RFdiffusion not covered by Nextflow parameters e.g. `'contigmap.inpaint_str=[B165-178]'`.                                                                                                                                                    |

---

## BindCraft Advanced Parameters

Advanced parameters to control the behaviour of BindCraft.

| Parameter                   | Default     | Description                                                                                                                                                                                                                                        |
| --------------------------- | ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `bc_chains`                 | null        | Optional specification of input PDB chains. Other chains will be ignored during design. Can be one or multiple chain IDs, in a comma-separated format e.g. 'A,C'. If null, will include all chains from input PDB.                                 |
| `bc_design_protocol`        | `'default'` | Which binder design protocol to run? "default" is recommended. "betasheet" promotes the design of more beta sheeted proteins, but requires more sampling. "peptide" is optimised for helical peptide binders.                                      |
| `bc_template_protocol`      | `'default'` | What target template protocol to use? "default" allows for limited amount of flexibility. "flexible" allows for greater target flexibility on both sidechain and backbone level.                                                                   |
| `bc_omit_AAs`               | `'C'`       | Residue types to avoid during sequence design (comma separated list, one letter case-insensitive) (default: 'C') e.g. 'C,H'. These residue types can still occur if no other options are possible in the position.                                 |
| `bc_fix_interface_residues` | true        | Whether to preserve/fix the interface residues of the binder designed by BindCraft during the subsequent sequence design stage (Recommended to leave as true).                                                                                     |
| `bc_advanced_json`          | true        | Path to a custom advanced settings json file. Not recommended unless you know what you are doing. Some parameters will be ignored as they apply to BindCraft routines downstream of fold design that are not implemented here e.g. MPNN parameters |

---

## ProteinMPNN-FastRelax Advanced Parameters

Advanced parameters to control the behaviour of ProteinMPNN-FastRelax.

| Parameter                           | Default      | Description                                                                                                                                                                             |
| ----------------------------------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mpnn_omitAAs`                      | `'CX'`       | Residue types to exclude during design (one-letter code, case insensitive).                                                                                                             |
| `mpnn_temperature`                  | 0.1          | Temperature for sequence sampling; higher values increase diversity. Recommended to lower significantly for binders to improve success (e.g., 0.0001).                                  |
| `mpnn_checkpoint_type`              | `'soluble'`  | Checkpoint selection: `'vanilla'`, `'soluble'` or `'hyper'`.                                                                                                                            |
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

| Parameter                      | Default | Description                                                                                                  |
| ------------------------------ | ------- | ------------------------------------------------------------------------------------------------------------ |
| `fampnn_fix_target_sidechains` | false   | Fix target side-chain positions when performing binder sequence design (default: false)                      |
| `fampnn_psce_threshold`        | 0.3     | Will only keep sidechains below this PSCE threshold during design. Null means no filtering.                  |
| `fampnn_temperature`           | 0.1     | Temperature for sampling; higher increases sequence diversity. Recommended lower for binders (e.g., 0.0001). |
| `fampnn_exclude_cys`           | true    | Exclude cysteine residues from design.                                                                       |
| `fampnn_extra_config`          | null    | Additional raw configuration passed to Full-Atom MPNN. e.g. 'timestep_schedule.num_steps=100 seq_only=true'. |

---

## Prediction Advanced Parameters

| Parameter                  | Default | Description                                                                                                                                                                          |
| -------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `uncropped_target_pdb`     | null    | Path to uncropped target PDB for prediction (e.g., full complex).                                                                                                                    |
| `af2_initial_guess`        | true    | Use an initial guess for target chains in AlphaFold2 structure predictions.                                                                                                          |
| `af2_extra_config`         | null    | Additional raw parameters passed to AlphaFold2 Initial-Guess. e.g. '-recycle 3'                                                                                                      |
| `boltz_recycling_steps`    | 3       | Number of recycling steps in Boltz-2 predictions.                                                                                                                                    |
| `boltz_diffusion_samples`  | 1       | Number of diffusion samples in Boltz-2 predictions.                                                                                                                                  |
| `boltz_sampling_steps`     | 200     | Number of sampling steps in Boltz-2 predictions.                                                                                                                                     |
| `boltz_use_potentials`     | false   | Use physics-based potentials during inference to improve physical plausibility of predictions (also known as Boltz-2x). Increases prediction time.                                   |
| `boltz_use_templates`      | true    | Enable template-guided structure prediction. In binder modes, the target (chain B) is used as template, similar to AlphaFold2 Initial-Guess. In monomer modes, chain A is templated. |
| `boltz_template_force`     | false   | Enforce template structure with potential during prediction. Requires `boltz_use_templates = true`.                                                                                  |
| `boltz_template_threshold` | null    | Distance threshold in Angstroms for template deviation. If specified with `boltz_template_force = true`, constrains prediction near template.                                        |
| `boltz_input_msa` | null | Path to a multiple sequence alignment (.a3m format) for the input PDB. For binder modes, MSA is applied to chain B (target), and for monomer modes, MSA is applied to chain A. MSA must match entire sequence of chain or will be ignored by Boltz. If null, single-sequence mode is used (msa: empty in YAML). e.g. 'lib/pdl1_msa.a3m' |
| `boltz_extra_config`       | null    | Additional raw parameters for Boltz-2 predictions. e.g. '--msa_pairing_strategy complete'                                                                                            |

---

## Filtering Parameters <a name="params-filter"></a>

Due to the inherently stochastic nature of protein design, often we see problematic results during the pipeline. It can save computation time to discard designs mid-pipeline that fail to meet success criteria. We have implemented four filtering stages that can be used to reject poor designs:

- **Fold Filtering** - Filters designs according to the number of secondary structure elements and radius of gyration.
- **Sequence Filtering** - Filters designs according to the score of the generated sequence
- **AlphaFold2/Boltz-2 Filtering** - Filters designs according to the quality of the structure prediction
- **Analysis Filtering** - Filters designs according to detailed biophysical metrics calculated by PyRosetta and BioPython, including interface quality, energy, and sequence properties

The most powerful predictors of experimental success are structure prediction metrics, but some metrics are more effective than others. Here are some recommended filters for binder design from the literature and their corresponding parameters in ProteinDJ:

| Parameter                     | RFdiffusion paper<sup>1</sup> | BindCraft paper <sup>2</sup> | AlphaProteo whitepaper<sup>3</sup> |
| ----------------------------- | ----------------------------- | ---------------------------- | ---------------------------------- |
| `af2_max_pae_interaction`     | 10                            | 10.5                         | 7                                  |
| `af2_min_plddt_overall`       | 80                            | 80                           | 90                                 |
| `af2_max_rmsd_binder_bndaln`  | 1                             |                              | 1.5                                |
| `af2_max_rmsd_binder_tgtaln`  |                               | 6                            |                                    |
| `boltz_max_rmsd_overall`      |                               |                              | 2.5                                |
| `boltz_min_ptm_binder`        |                               |                              | 0.8                                |
| `pr_min_intface_shpcomp`      |                               | 0.6                          |                                    |
| `pr_min_intface_hbonds`       |                               | 3                            |                                    |
| `pr_max_intface_unsat_hbonds` |                               | 4                            |                                    |
| `pr_max_surfhphobics`         |                               | 35                           |                                    |

<sup> 1. Watson, J.L. et al. Nature 620, 1089–1100 (2023). https://doi.org/10.1038/s41586-023-06415-8; 2. Pacesa, M. et al. Nature 646, 483-492 (2025). https://doi.org/10.1038/s41586-025-09429-6 3. Zambaldi, V. et al. arXiv (2024). https://doi.org/10.48550/arXiv.2409.08022
</sup>

We recommend disabling other filters for small-scale and pilot experiments, and using these results to decide on values to use for filtering large-scale runs. Note that BindCraft has built-in filtering of designs and will automatically reject designs that meet any of the following criteria:

- Low confidence (pLDDT < 0.7)
- Severe clashes (clashes detected between C-alpha atoms)
- Insufficient contact between binder and target (less than three residues contacting the target)

If a design fails to meet these criteria, BindCraft will generate a new design until it finds one that passes. This can lead to long run times compared to RFdiffusion but tends to give binder designs that are more likely to succeed in the subsequent Structure Prediction stage.

### Fold Filtering Parameters

Fold Filtering Parameters that can be used to filter designs by RFdiffusion and BindCraft according to secondary structure and radius of gyration. Metrics are calculated on the binder chain only in binder design modes, otherwise all chains are used in calculations.

| Parameter          | Description                                                             |
| ------------------ | ----------------------------------------------------------------------- |
| `fold_min_helices` | Minimum number of alpha-helices required.                               |
| `fold_max_helices` | Maximum number of alpha-helices allowed.                                |
| `fold_min_strands` | Minimum number of beta-strands required.                                |
| `fold_max_strands` | Maximum number of beta-strands allowed.                                 |
| `fold_min_ss`      | Minimum number of secondary structure elements (α-helices + β-strands). |
| `fold_max_ss`      | Maximum number of secondary structure elements (α-helices + β-strands). |
| `fold_min_rog`     | Minimum radius of gyration (Å).                                         |
| `fold_max_rog`     | Maximum radius of gyration (Å).                                         |

### Sequence Filtering Parameters

Sequence Filtering Parameters for ProteinMPNN and Full-Atom MPNN. Recommended to disable unless you know what you are doing.

| Parameter         | Description                                                                                                      |
| ----------------- | ---------------------------------------------------------------------------------------------------------------- |
| `mpnn_max_score`  | Maximum MPNN score (negative log likelihood). A lower score means ProteinMPNN is more confident in the sequence. |
| `fampnn_max_psce` | Max PSCE score for designed side-chains. A lower score means FAMPNN is more confident in the sequence.           |

### AlphaFold2 Filtering Parameters

AlphaFold2 Initial-Guess Filtering Parameters.

| Parameter                    | Description                                                                                    |
| ---------------------------- | ---------------------------------------------------------------------------------------------- |
| `af2_max_pae_interaction`    | Max predicted aligned error for interactions                                                   |
| `af2_max_pae_overall`        | Max predicted aligned error for all chains                                                     |
| `af2_max_pae_binder`         | Max predicted aligned error for binder                                                         |
| `af2_max_pae_target`         | Max predicted aligned error for target                                                         |
| `af2_max_rmsd_overall`       | Max C-alpha RMSD between AF2 prediction and input design when all chains are aligned           |
| `af2_max_rmsd_binder_bndaln` | Max binder C-alpha RMSD between AF2 prediction and input design when binder chains are aligned |
| `af2_max_rmsd_binder_tgtaln` | Max binder C-alpha RMSD between AF2 prediction and input design when target chains are aligned |
| `af2_max_rmsd_target`        | Max target C-alpha RMSD between AF2 prediction and input design when target chains are aligned |
| `af2_min_plddt_overall`      | Min average pLDDT score for all chains                                                         |
| `af2_min_plddt_binder`       | Min pLDDT score required for binder                                                            |
| `af2_min_plddt_target`       | Min pLDDT score required for target                                                            |

### Boltz-2 Filtering Parameters

Boltz-2 Filtering Parameters.

| Parameter                   | Description                                                                                       |
| --------------------------- | ------------------------------------------------------------------------------------------------- |
| `boltz_max_rmsd_overall`    | Max C-alpha RMSD between all chains of Boltz-2 prediction and input design                        |
| `boltz_max_rmsd_binder`     | Max C-alpha RMSD between binder chains of Boltz-2 prediction and input design. Binder modes only. |
| `boltz_max_rmsd_target`     | Max C-alpha RMSD between target chains of Boltz-2 prediction and input design. Binder modes only. |
| `boltz_min_conf_score`      | Minimum confidence score of the prediction                                                        |
| `boltz_min_ipSAE_min`       | Minimum value allowed for the minimum interaction prediction Score from Aligned Errors (ipSAE) of target and binder chains (0 to 1).      |
| `boltz_min_LIS`             | Minimum Local Interaction Score (> 0)                                           |
| `boltz_min_pDockQ2_min`     | Minimum value allowed for the minimum predicted DockQ Score v2 of target and binder chains (0 to 1).                                       |
| `boltz_max_pae_interaction`   | Maximum predicted aligned error at interaction interfaces (0 to ~30 Å)                        |
| `boltz_min_ptm`             | Minimum predicted template modelling score of the prediction                                      |
| `boltz_min_ptm_binder`      | Minimum predicted template modelling score of the binder chain. Binder modes only.                |
| `boltz_min_ptm_target`      | Minimum predicted template modelling score of the target chain. Binder modes only.                |
| `boltz_min_ptm_interface`   | Minimum predicted template modelling score of the prediction interface                            |
| `boltz_min_plddt`           | Minimum pLDDT score of the prediction                                                             |
| `boltz_min_plddt_interface` | Minimum pLDDT score of the prediction interface                                                   |
| `boltz_max_pde`             | Maximum predicted distance error of the prediction                                                |
| `boltz_max_pde_interface`   | Maximum predicted distance error for the prediction interface                                     |

### Analysis Filtering Parameters

Analysis Filtering Parameters for the final Analysis stage using PyRosetta and BioPython. These metrics provide detailed biophysical characterization of the predicted structures, particularly useful for binder design. Note that interface metrics are only calculated for binder design modes.

| Parameter                     | Description                                                                                   |
| ----------------------------- | --------------------------------------------------------------------------------------------- |
| `pr_min_helices`              | Minimum number of alpha-helices in predicted structure                                        |
| `pr_max_helices`              | Maximum number of alpha-helices in predicted structure                                        |
| `pr_min_strands`              | Minimum number of beta-strands in predicted structure                                         |
| `pr_max_strands`              | Maximum number of beta-strands in predicted structure                                         |
| `pr_min_total_ss`             | Minimum total secondary structure elements (α-helices + β-strands) in predicted structure     |
| `pr_max_total_ss`             | Maximum total secondary structure elements (α-helices + β-strands) in predicted structure     |
| `pr_min_rog`                  | Minimum radius of gyration (Å) of predicted structure                                         |
| `pr_max_rog`                  | Maximum radius of gyration (Å) of predicted structure                                         |
| `pr_min_intface_bsa`          | Minimum buried surface area (Å²) at the binding interface                                     |
| `pr_min_intface_shpcomp`      | Minimum shape complementarity of interface (0-1 scale; 1 is optimal)                          |
| `pr_min_intface_hbonds`       | Minimum number of hydrogen bonds at the interface                                             |
| `pr_max_intface_unsat_hbonds` | Maximum number of buried, unsatisfied hydrogen bonds at the interface                         |
| `pr_max_intface_deltag`       | Maximum solvation free energy gain at interface (Rosetta Energy Units; lower is better)       |
| `pr_max_intface_deltagtobsa`  | Maximum ratio of delta-G to buried surface area                                               |
| `pr_min_intface_packstat`     | Minimum packing quality of the interface (0-1 scale; higher is better)                        |
| `pr_max_tem`                  | Maximum total energy metric score (Rosetta Energy Units; lower indicates more stable designs) |
| `pr_max_surfhphobics`         | Maximum percentage of hydrophobic residues exposed on the surface                             |
| `pr_max_sap`                  | Maximum mean residue Spatial Aggregation Propensity for monomer/binder (solubility prediction; lower is better) |
| `pr_max_sap_complex`          | Maximum mean residue Spatial Aggregation Propensity for binder in complex (solubility prediction; lower is better) |
| `seq_min_ext_coef`            | Minimum extinction coefficient at 280nm (M⁻¹cm⁻¹)                                             |
| `seq_max_ext_coef`            | Maximum extinction coefficient at 280nm (M⁻¹cm⁻¹)                                             |
| `seq_min_pi`                  | Minimum isoelectric point (pI) of the sequence                                                |
| `seq_max_pi`                  | Maximum isoelectric point (pI) of the sequence                                                |

## Cluster Parameters

The cluster parameters may need adjusting depending on your HPC setup and available hardware. You must ensure that the paths to the containers and models for RFdiffusion, AlphaFold2 and Boltz-2 are valid and contain the expected files (see ['Installation Guide'](installation.md)).

| Parameter      | Example Value                  | Description                                 |
| -------------- | ------------------------------ | ------------------------------------------- |
| `rfd_models`   | `"${projectDir}/models/rfd"`   | Path to the RFdiffusion model checkpoints.  |
| `af2_models `  | `"${projectDir}/models/af2"`   | Path to the AlphaFold2 models.              |
| `boltz_models` | `"${projectDir}/models/boltz"` | Path to the Boltz-2 models.                 |
| `gpu_model`    | `'A30'`                        | GPU model to request, e.g., 'A30'.          |
| `gpus`         | `1`, `2`, `4`, `8`             | Number of GPUs to request                   |
| `cpus_per_gpu` | `8`, `12`                      | Number of CPUs to request per GPU           |
| `memory_gpu`   | `'24GB'`, `'48GB'`             | Memory to request for GPU jobs              |
| `cpus`         | `12`, `24`                     | Number of CPUs to request for CPU-only jobs |
| `memory_cpu`   | `'24GB'`, `'32GB'`             | Memory for request for CPU-only jobs        |

---

**🌟 Notes:**

- Parameters set to `null` indicate optional or user-defined inputs.
- Ensure paths are updated to your environment when running ProteinDJ.
- Filtering parameters can be set to `null` to disable filtering.

---

[⬅️ Back to Main README](../README.md)
