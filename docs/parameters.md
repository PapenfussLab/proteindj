[🏠 ProteinDJ](../README.md) > Parameter Guide

# ProteinDJ Parameter Guide

This guide summarizes the key parameters used in ProteinDJ's Nextflow pipeline configuration. Parameters are essential for controlling the design mode, input/output locations, model choices, advanced options, filtering thresholds, and cluster resource allocation.

---

## Essential Parameters

These parameters are required for ProteinDJ and are used by every mode.

| Parameter         | Default         | Description                                                                                                                                                                                                         |
| ----------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `design_mode`     | null            | Pipeline mode. Choose from: `rfd_denovo`, `rfd_foldcond`, `rfd_motifscaff`, `rfd_partialdiff` (each automatically runs as monomer or binder design depending on `input_pdb`), `bindcraft_denovo`, `boltzgen_denovo`, or `boltzgen_motifscaff` |
| `num_designs`     | 8               | Number of designs to generate using RFdiffusion or Bindcraft                                                                                                                                                        |
| `seqs_per_design` | 8               | Number of sequences to generate per design                                                                                                                                                                          |
| `out_dir`         | `./pdj_results` | Output directory for results. Existing results will be overwritten                                                                                                                                                  |

---

## Mode-Specific Parameters

These parameters are used for some of the ProteinDJ modes.

| Parameter                         | Default | Description                                                                                                                                                                                                                                                                                                                       | Required for Modes                                                                                                         |
| --------------------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `design_length`                   | null    | Design length of binder or monomer. e.g. `'60'` or `'60-70'`. Pipeline will randomly sample different sizes between these values.                                                                                                                                                                                                 | `bindcraft_denovo`, `rfd_denovo`, `boltzgen_denovo`                                                                      |
| `input_pdb`                       | null    | Path to input PDB file (e.g. target for binders, `'./target.pdb'`). Required for all binder modes.                                                                                                                                                                                                                                   | `bindcraft_denovo`, `boltzgen_motifscaff`, `rfd_motifscaff`, `rfd_partialdiff` |
| `hotspot_residues`                | null    | Hotspot residues on target chain(s) for binder design. Comma-separated tokens, each a chain-qualified residue (`A56`), chain-qualified range (`A115-120`), or bare chain ID meaning the whole chain (`B`), e.g. `A56,A115-120,B`. A chain identifier is always required. Optional for `rfd_denovo`/`rfd_foldcond`, `bindcraft_denovo`, `boltzgen_denovo`/`boltzgen_motifscaff` binder design modes.                                                                                                                                                                                                                     |                                                                                                                            |
| `motifscaff_spec`                   | null    | Chain A architecture spec: an ordered, comma-separated list of chain-A keep tokens (`A<start>-<end>`/`A<n>`) and insert tokens (bare `<n>`/`<min>-<max>`). Used by `boltzgen_motifscaff` and `rfd_motifscaff`.                                                                                                                                                                                                        |                                                                                                    |
| `motifscaff_inpaint_seq`            | null    | Chain A residues, within those kept by `motifscaff_spec`, whose sequence is allowed to change while structure stays fixed/conditioned. Comma-separated e.g. `'A10-50,A60'`. Used by `boltzgen_motifscaff` and `rfd_motifscaff`.                                                                                                                                                                                                                                  |                                                                                                    |
| `flexible_residues`               | null    | Residues whose structure should NOT be conditioned on, useful for disordered/flexible regions e.g. loops or IDPs. Comma-separated tokens, each a chain-qualified residue e.g. `'A16'`, a chain-qualified range e.g. `'B10-13'`, or a bare chain ID e.g. `'C'` meaning the whole chain. Used by `boltzgen_denovo`, `boltzgen_motifscaff`, `rfd_motifscaff`, `rfd_partialdiff`, and binder design with `rfd_denovo`/`rfd_foldcond`.                                                                                                              |                                                                                                    |
| `rfd_foldcond_scaffold_dir`                | null    | Directory containing scaffold secondary structure and block adjacency files (e.g. `'./binderscaffolds/scaffolds_assorted'`). Required for fold conditioning modes.                                                                                                                                                                | `rfd_foldcond`                                                                                      |
| `rfd_foldcond_mask_loops`                  | true    | Whether to ignore loops during scaffold secondary structure constraints. Optional for `rfd_foldcond` mode.                                                                                                                                                                                             |                                                                                                                            |
| `rfd_motifscaff_length`                      | null    | Length of diffused chain during motif scaffolding, e.g. `100-100` or `100-120`. Optional for motif scaffolding mode `rfd_motifscaff`.                                                                                                                                                                   |                                                                                                                            |
| `rfd_partialdiff_timesteps` | null    | Number of timesteps for partial diffusion (1-50) e.g. 20. Required for partial diffusion mode.                                                                                                                                                                                                                                   | `rfd_partialdiff`                                                                                |
| `rfd_partialdiff_fixed_seq`           | null    | Design chain residues whose sequence should be kept fixed during partial diffusion, e.g. `A17-19,A143-145`. References the design chain only (chain A in binder mode, the sole chain in monomer mode). If null, the sequence is fully redesignable. Optional for `rfd_partialdiff` mode.                                       |                                                                                                   |
| `rfd_partialdiff_spec`                | null    | Chain A partial diffusion spec: an ordered, comma-separated list of chain-A keep tokens (`A<start>-<end>`/`A<n>`, left fully fixed) and fixed-count diffuse tokens (bare `<n>` only, no ranges since length cannot change), e.g. `A1-19,59`. If null, the whole chain is partially diffused.                                       |                                                                                                   |

### Parameter Requirements by Mode

`-` = ignored

#### Monomer design

| Parameter                       | rfd_denovo | rfd_foldcond | rfd_motifscaff | rfd_partialdiff | boltzgen_denovo | boltzgen_motifscaff |
| -------------------------------- | ---------- | ------------ | --------------- | ---------------- | ---------------- | ------------------- |
| design_length                    | Required   | -            | -               | -                 | Required          | -                   |
| input_pdb                        | -          | -            | Required        | Required          | -                 | Required            |
| hotspot_residues                 | -          | -            | -               | -                 | -                 | -                   |
| bg_not_binding_residues          | -          | -            | -               | -                 | -                 | -                   |
| motifscaff_spec                    | -          | -            | _Optional_      | -                 | -                 | _Optional_          |
| motifscaff_inpaint_seq             | -          | -            | _Optional_      | -                 | -                 | _Optional_          |
| flexible_residues                | -          | -            | _Optional_      | _Optional_        | -                 | _Optional_          |
| rfd_foldcond_scaffold_dir                 | -          | Required     | -               | -                 | -                 | -                   |
| rfd_foldcond_mask_loops                   | -          | _Optional_   | -               | -                 | -                 | -                   |
| rfd_motifscaff_length                       | -          | -            | _Optional_      | -                 | -                 | -                   |
| rfd_partialdiff_timesteps  | -          | -            | -               | Required          | -                 | -                   |
| rfd_partialdiff_fixed_seq            | -          | -            | -               | _Optional_        | -                 | -                   |
| rfd_partialdiff_spec                 | -          | -            | -               | _Optional_        | -                 | -                   |

#### Binder design

| Parameter                       | rfd_denovo | rfd_foldcond | rfd_motifscaff | rfd_partialdiff | bindcraft_denovo | boltzgen_denovo | boltzgen_motifscaff |
| -------------------------------- | ---------- | ------------ | --------------- | ---------------- | ----------------- | ----------------- | ------------------- |
| design_length                    | Required   | -            | -               | -                 | Required           | Required           | -                   |
| input_pdb                        | Required   | Required     | Required        | Required          | Required           | Required           | Required            |
| hotspot_residues                 | _Optional_ | _Optional_   | -               | -                 | _Optional_         | _Optional_         | _Optional_          |
| bc_chains                        | -          | -            | -               | -                 | _Optional_         | -                  | -                   |
| bc_design_protocol               | -          | -            | -               | -                 | _Optional_         | -                  | -                   |
| bc_template_protocol             | -          | -            | -               | -                 | _Optional_         | -                  | -                   |
| bc_omit_AAs                      | -          | -            | -               | -                 | _Optional_         | -                  | -                   |
| bc_fix_interface_residues        | -          | -            | -               | -                 | _Optional_         | -                  | -                   |
| bc_advanced_json                 | -          | -            | -               | -                 | _Optional_         | -                  | -                   |
| bg_not_binding_residues          | -          | -            | -               | -                 | -                  | _Optional_         | _Optional_          |
| motifscaff_spec                    | -          | -            | _Optional_      | -                 | -                  | -                  | _Optional_          |
| motifscaff_inpaint_seq             | -          | -            | _Optional_      | -                 | -                  | -                  | _Optional_          |
| flexible_residues                | _Optional_ | _Optional_   | _Optional_      | _Optional_        | -                  | _Optional_         | _Optional_          |
| rfd_foldcond_scaffold_dir                 | -          | Required     | -               | -                 | -                  | -                  | -                   |
| rfd_foldcond_mask_loops                   | -          | _Optional_   | -               | -                 | -                  | -                  | -                   |
| rfd_motifscaff_length                       | -          | -            | _Optional_      | -                 | -                  | -                  | -                   |
| rfd_partialdiff_timesteps  | -          | -            | -               | Required          | -                  | -                  | -                   |
| rfd_partialdiff_fixed_seq            | -          | -            | -               | _Optional_        | -                  | -                  | -                   |
| rfd_partialdiff_spec                 | -          | -            | -               | _Optional_        | -                  | -                  | -                   |

---

## Workflow Advanced Parameters

These parameters control the workflow of ProteinDJ.

| Parameter            | Default | Description                                                                                                                                                                                                              |
| -------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `seq_method`         | 'mpnn'  | Sequence design method. Options: `'mpnn'` (ProteinMPNN Fast-Relax) or `'fampnn'` (Full-Atom MPNN).                                                                                                                       |
| `pred_method`        | 'af2'   | Structure prediction method. Options: `'af2'` (AlphaFold2 Initial-Guess), `'boltz'` (Boltz-2), or `'af2_boltz'` (runs AF2 first, filters, then runs Boltz-2 on the survivors for a final, more accurate prediction). |
| `zip_pdbs`           | true    | Whether to compress output final designs in a tar.gz archive. If false, results will be output as uncompressed PDB files.                                                                                                |
| `rank_designs`       | true   | Whether to rank final designs by prediction quality metrics and output `ranked_designs.csv` and `ranked_designs` folder (or `.tar.gz`). |
| `ranking_metric`     | null   | Structure prediction metric to use for ranking designs (e.g., `'af2_pae_interaction'`, `'boltz_iptm'`, `'boltz_ipSAE_min'`). Must match prediction method prefix (`af2_*` or `boltz_*`; either prefix is accepted when `pred_method` is `'af2_boltz'`). If null, defaults are set depending on design_mode and pred_method: binder modes use `'af2_pae_interaction'` (AlphaFold2) or `'boltz_ipSAE_min'` (Boltz-2 or af2_boltz); monomer modes use `'af2_plddt_overall'` (AlphaFold2) or `'boltz_ptm'` (Boltz-2 or af2_boltz). See [Metrics Guide](metrics.md) for all available metrics. |
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
| `rfd_ckpt_override` | null    | Overrides the diffusion model checkpoint. Options include `'active_site'`, `'base'`, `'base_epoch8'`, `'complex_base'`, `'complex_beta'`, `'complex_fold_base'`, `'inpaint_seq'`, `'inpaint_seq_fold'`. Not generally recommended except for `rfd_denovo` binder design with `'complex_beta'`. `rfd_motifscaff` binder design already defaults to `'complex_base'` internally (unless `motifscaff_inpaint_seq`/`flexible_residues` is set) since RFdiffusion's own checkpoint auto-selection would otherwise fall back to the monomer-only checkpoint. |
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

## BoltzGen Advanced Parameters

Advanced parameters to control the behaviour of BoltzGen (see [BoltzGen Design mode guide](modes.md#boltzgendesign) for detailed usage examples).

| Parameter                  | Default | Description                                                                                                                                                                                                                                                                                                                                                                                                             | Applies to        |
| --------------------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| `bg_not_binding_residues`   | null    | Target residues the binder should NOT contact. Comma-separated tokens, each a chain-qualified residue e.g. `'A200'`, a chain-qualified range e.g. `'A210-215'`, or a bare chain ID e.g. `'B'` meaning the whole chain.                                                                                | `boltzgen_denovo`, `boltzgen_motifscaff` |

---

## ProteinMPNN-OpenMMRelax Advanced Parameters

Advanced parameters to control the behaviour of ProteinMPNN-OpenMMRelax.

| Parameter                           | Default      | Description                                                                                                                                                                             |
| ----------------------------------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `mpnn_omitAAs`                      | `'CX'`       | Residue types to exclude during design (one-letter code, case insensitive).                                                                                                             |
| `mpnn_temperature`                  | 0.1          | Temperature for sequence sampling; higher values increase diversity. Recommended to lower significantly for binders to improve success (e.g., 0.0001).                                  |
| `mpnn_checkpoint_type`              | `'soluble'`  | Checkpoint selection: `'vanilla'`, `'soluble'` or `'hyper'`.                                                                                                                            |
| `mpnn_checkpoint_model`             | `'v_48_020'` | Checkpoint model variant indicating backbone noise level used during training. e.g. 'v_48_020' noised with 0.20Å Gaussian noise. Choose from 'v_48_002','v_48_010,'v_48_020','v_48_030' |
| `mpnn_backbone_noise`               | 0            | Std dev Gaussian noise added to backbone atoms. 0 = none, 0.1-0.3 = mild perturbation.                                                                                                  |
| `mpnn_relax_max_cycles`             | 0            | Max fast relaxation cycles per sequence; 0 disables OpenMMRelax functionality.                                                                                                            |
| `mpnn_relax_seqs_per_cycle`         | 1            | Number of sequences generated between relaxation steps; best scored sequence is kept.                                                                                                   |
| `mpnn_relax_output`                 | false        | Whether to run relaxation before saving output.                                                                                                                                         |
| `mpnn_relax_convergence_rmsd`       | 0.2          | RMSD early convergence threshold for relaxation cycles. Design is considered converged if the C-alpha RMSD (Å) between cycles is <= this threshold.                                     |
| `mpnn_relax_convergence_score`      | 0.1          | Score early convergence threshold for relaxation cycles. Design is considered converged if the improvement in score between cycles is <= this threshold.                                |
| `mpnn_relax_convergence_max_cycles` | 1            | Design is considered converged if it meets both convergence criteria for n consecutive cycles.                                                                                          |
| `mpnn_extra_config`                 | null         | Additional raw configuration passed to ProteinMPNN-OpenMMRelax. e.g. `-protein_features="full"'`                                                                                          |

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
| `af2_jax_compilation_cache_dir` | null | Optional existing path for the persistent JAX compilation cache used by AlphaFold2. Relative paths are resolved against the workflow launch directory. The path is bind-mounted into the container. When null, `JAX_COMPILATION_CACHE_DIR` is not set. |
| `boltz_recycling_steps`    | 3       | Number of recycling steps in Boltz-2 predictions.                                                                                                                                    |
| `boltz_diffusion_samples`  | 1       | Number of diffusion samples in Boltz-2 predictions.                                                                                                                                  |
| `boltz_sampling_steps`     | 200     | Number of sampling steps in Boltz-2 predictions.                                                                                                                                     |
| `boltz_use_potentials`     | false   | Use physics-based potentials during inference to improve physical plausibility of predictions (also known as Boltz-2x). Increases prediction time.                                   |
| `boltz_use_templates`      | true    | Enable template-guided structure prediction. In binder modes, the target (all non-binder chains, e.g. chain B, or B,C... if there are multiple target chains) is used as template, similar to AlphaFold2 Initial-Guess. In monomer modes, chain A is templated. |
| `boltz_template_force`     | false   | Enforce template structure with potential during prediction. Requires `boltz_use_templates = true`.                                                                                  |
| `boltz_template_threshold` | null    | Distance threshold in Angstroms for template deviation. If specified with `boltz_template_force = true`, constrains prediction near template.                                        |
| `boltz_input_msa` | null | Comma-separated list of one-or-more multiple sequence alignment files (.a3m format) for the input PDB. For binder modes, applied to all non-binder (target) chain(s): a single file is applied to all target chains (e.g. an identical-sequence multi-chain target), or multiple files are matched positionally one-per-chain in chain order. For monomer modes, applied to chain A. MSA must match entire sequence of each chain it is applied to or will be ignored by Boltz for that chain. If null, single-sequence mode is used (msa: empty in YAML). e.g. 'lib/pdl1_msa.a3m' or 'lib/chainB_msa.a3m,lib/chainC_msa.a3m' |
| `boltz_extra_config`       | null    | Additional raw parameters for Boltz-2 predictions. e.g. '--msa_pairing_strategy complete'                                                                                            |
| `boltz_predict_unbound_binder` | false | Also predict the binder alone (no target) as a target-free monomer prediction with Boltz-2, to check whether the binder folds consistently to the design in the absence of the target. Binder design modes only. Runs as a pure single-sequence prediction (no template, no MSA). |

---

## Filtering Parameters <a name="params-filter"></a>

Due to the inherently stochastic nature of protein design, often we see problematic results during the pipeline. It can save computation time to discard designs mid-pipeline that fail to meet success criteria. We have implemented four filtering stages that can be used to reject poor designs:

- **Fold Filtering** - Filters designs according to the number of secondary structure elements and radius of gyration.
- **Sequence Filtering** - Filters designs according to the score of the generated sequence
- **AlphaFold2/Boltz-2 Filtering** - Filters designs according to the quality of the structure prediction
- **Analysis Filtering** - Filters designs according to detailed biophysical metrics calculated by arpeggia, PRODIGY, and BioPython on the energy-minimized structure, including interface quality, energy, and sequence properties

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
| `seq_min_ext_coef`            | Minimum extinction coefficient at 280nm (M⁻¹cm⁻¹)                                             |
| `seq_max_ext_coef`            | Maximum extinction coefficient at 280nm (M⁻¹cm⁻¹)                                             |
| `seq_min_pi`                  | Minimum isoelectric point (pI) of the sequence                                                |
| `seq_max_pi`                  | Maximum isoelectric point (pI) of the sequence                                                |

### AlphaFold2 Filtering Parameters

AlphaFold2 Initial-Guess Filtering Parameters.

| Parameter                    | Description                                                                                    |
| ---------------------------- | ---------------------------------------------------------------------------------------------- |
| `af2_min_iptm`               | Min predicted interface TM-score (ipTM) between binder and target. |
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
| `boltz_min_iptm`   | Minimum predicted template modelling score of the prediction interface                            |
| `boltz_min_plddt`           | Minimum pLDDT score of the prediction                                                             |
| `boltz_min_iplddt` | Minimum pLDDT score of the prediction interface                                                   |
| `boltz_max_pde`             | Maximum predicted distance error of the prediction                                                |
| `boltz_max_ipde`   | Maximum predicted distance error for the prediction interface                                     |
| `boltz_max_unbound_rmsd`    | Max C-alpha RMSD between the unbound (target-free) binder prediction and the design*|
| `boltz_min_unbound_conf_score` | Minimum confidence score of the unbound (target-free) binder prediction*|
| `boltz_min_unbound_ptm`     | Minimum pTM score of the unbound (target-free) binder prediction*|
| `boltz_min_unbound_plddt`   | Minimum pLDDT score of the unbound (target-free) binder prediction* |
| `boltz_max_unbound_pde`     | Maximum PDE score of the unbound (target-free) binder prediction* |

*Requires `boltz_predict_unbound_binder = true` to generate unbound binder predictions.

### Analysis Filtering Parameters

Analysis Filtering Parameters for the final Analysis stage using PDBFixer/OpenMM, arpeggia, PRODIGY, and BioPython. These metrics provide detailed biophysical characterization of the predicted structures, particularly useful for binder design. Note that interface metrics are only calculated for binder design modes.

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
| `pr_max_intface_deltag`       | Maximum predicted binding free energy at the interface (PRODIGY IC-NIS model, kcal/mol; more negative is better)       |
| `pr_max_intface_deltagtobsa`  | Maximum ratio of delta-G to buried surface area                                               |
| `pr_max_surfhphobics`         | Maximum percentage of hydrophobic residues exposed on the surface                             |
| `pr_max_sap`                  | Maximum mean residue Spatial Aggregation Propensity for monomer/binder (solubility prediction; lower is better) |
| `pr_max_sap_complex`          | Maximum mean residue Spatial Aggregation Propensity for binder in complex (solubility prediction; lower is better) |

## Cluster Parameters

The cluster parameters may need adjusting depending on your HPC setup and available hardware. You must ensure that the paths to the containers and models for RFdiffusion, AlphaFold2 and Boltz-2 are valid and contain the expected files (see ['Installation Guide'](installation.md)).

| Parameter      | Example Value                  | Description                                 |
| -------------- | ------------------------------ | ------------------------------------------- |
| `rfd_models`   | `"${projectDir}/models/rfd"`   | Path to the RFdiffusion model checkpoints.  |
| `af2_models `  | `"${projectDir}/models/af2"`   | Path to the AlphaFold2 models.              |
| `boltz_models` | `"${projectDir}/models/boltz"` | Path to the Boltz-2 models.                 |
| `mpnn_models`  | `"${projectDir}/models/mpnn"`  | Path to the ProteinMPNN model checkpoints.  |
| `bg_models`    | `"${projectDir}/models/boltzgen"` | Path to the BoltzGen model checkpoints.  |
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
