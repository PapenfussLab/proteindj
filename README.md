# ProteinDJ & BindSweeper

<img height="240" src="img/logo.png"><img height="240" alt="bindsweeper_logo" src="img/bindsweeper_logo.png" />

[Link to Publication](https://onlinelibrary.wiley.com/doi/10.1002/pro.70464)

ProteinDJ is a Nextflow pipeline for protein design that installs and utilises multiple external software packages, including:

- AlphaFold2 Initial Guess (from https://github.com/nrbennet/dl_binder_design)
- Arpeggia (from https://github.com/y1zhou/arpeggia)
- BindCraft (from https://github.com/martinpacesa/BindCraft) and FreeBindCraft (https://github.com/cytokineking/FreeBindCraft)
- BioPython (from https://biopython.org/)
- Boltz-2 (from https://github.com/jwohlwend/boltz)
- BoltzGen (from https://github.com/HannesStark/boltzgen)
- Full-Atom MPNN (from https://github.com/richardshuai/fampnn)
- HyperMPNN (from https://github.com/meilerlab/HyperMPNN)
- OpenMM (from https://openmm.org/) and PDBFixer (from https://github.com/openmm/pdbfixer)
- PRODIGY (from https://github.com/haddocking/prodigy)
- ProteinMPNN-FastRelax (from https://github.com/nrbennet/dl_binder_design)
- RFdiffusion (from https://github.com/RosettaCommons/RFdiffusion)

BindSweeper provides a convenient wrapper for ProteinDJ, enabling sweeping of different parameters for binder design.

If you find ProteinDJ or BindSweeper useful in your research, please [cite us](https://onlinelibrary.wiley.com/doi/10.1002/pro.70464) and the developers of the software listed above that make protein design pipelines like this possible. We have provided a list of citations [here](#citations).

<sup>_Logo image credit: Lyn Deng, Joshua Hardy_</sup>

## Table of Contents

- [Installation](#install)
- [Using ProteinDJ](#execution)
- [Tutorial - De novo binder design](#tutorial)
- [Advanced Parameters](#params)
- [Filtering Designs](#params-filter)
- [Metrics and Metadata](#metrics)
- [BindSweeper](#bindsweeper)
- [Appendices](#append)
  - [Known limitations](#limitations)
  - [Seqera Support](#seqera)
  - [Troubleshooting and common errors](#errors)
  - [Data used for benchmarking](#append-bench)
  - [Citations for software packages used in ProteinDJ](#citations)

## Installation <a name="install"></a>

> **Note: v3 is not backwards compatible with v2.** Containers and models have been updated, so a fresh installation is required.

ProteinDJ requires that [Apptainer](https://apptainer.org/docs/admin/main/installation.html) and [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ v24.04) are installed and accessible to your environment. For v26 of Nextflow onwards, you will need to export NXF_SYNTAX_PARSER=v1 to your shell environment, as ProteinDJ does not yet follow the strict syntax enforced in recent Nextflow versions.

First, clone the repo for ProteinDJ:

```
git clone https://github.com/PapenfussLab/proteindj
cd proteindj
```

Next, download the models for AF2, Boltz, RFdiffusion etc. (~16 GB) using the download script . This may take a while depending on your internet connection. Note that this only needs to be done once on a cluster as long as the files and containers are in a location that can be accessed by all users (see [Installation Guide](docs/installation.md) for more details):

```
bash scripts/download_models.sh
```

Apptainer will automatically fetch containers as needed during the Nextflow run and cache them locally to the location specified by the environment variable `NXF_APPTAINER_CACHEDIR`. If you would like to build containers locally, you can follow our [Installation Guide](docs/installation.md).

## Using ProteinDJ <a name="execution"></a>

The ProteinDJ consists of four stages:

1. Fold Design - Using RFdiffusion, BindCraft, or BoltzGen
2. Sequence Design - Using ProteinMPNN or Full-Atom MPNN
3. Structure Prediction - Using AlphaFold2 Initial Guess, Boltz-2, or both sequentially
4. Analysis and Reporting - Using a combination of tools including PDBFixer/OpenMM, arpeggia, PRODIGY, and BioPython

<img src="img/pipelineoverview.png" height="200">

Due to the creative nature of protein design and the complexity of RFdiffusion there are many ways you can use ProteinDJ. To help with delineating this, we have created design modes for ProteinDJ. Each mode is described in detail in our [Guide to Design Modes](docs/modes.md), but for now, here's a quick summary of each one with a simple illustration of each mode in action:

- **rfd_denovo** – diffusion of new monomers/binders from noise
- **rfd_foldcond** – diffusion of new monomers/binders with fold-conditioning on scaffolds/templates
- **rfd_motifscaff** – inpainting/extension of input monomers, or diffusion of binding motifs in input scaffolds / diffusion of binding motifs in scaffolds
- **rfd_partialdiff** – partial diffusion of an input monomer or binder
- **bindcraft_denovo** - hallucination of new binders using BindCraft
- **boltzgen_denovo** – generative design of new binders against a target using BoltzGen
- **boltzgen_motifscaff** – redesign/rediffusion of an existing binder using BoltzGen

<img src="img/modes_overview.png" height="720">

All the settings and parameters for ProteinDJ can be found in the `nextflow.config` file. It contains a lot of optional parameters, but there are 4 essential parameters to pay attention to: the protein design mode (`design_mode`), the number of designs (`num_designs`) and sequences you want to generate (`seqs_per_design`), and the output directory path (`out_dir`).

| Parameter         | Example Value      | Description                                                                                                                                                                                     |
| ----------------- | ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `design_mode`        | `'rfd_denovo'` | Pipeline mode. Choose from 'rfd_denovo', 'rfd_foldcond', 'rfd_motifscaff', 'rfd_partialdiff', 'bindcraft_denovo', 'boltzgen_denovo', or 'boltzgen_motifscaff'. Each automatically runs as monomer or binder design depending on inputs. |
| `num_designs` | `8`                | Number of designs to generate using BindCraft, BoltzGen, or RFdiffusion.                                                                                                                                                |
| `seqs_per_design` | `8`                | Number of sequences to generate per design.                                                                                                                                         |
| `out_dir`         | `'./pdj_results'`  | Output directory for results. Existing results in this directory will be overwritten.                                                                                                           |

If a parameter has a prefix, it indicates that it is specific to an underlying program e.g. 'bg_' for BoltzGen. Other parameters are shared by multiple programs `input_pdb`, `hotspot_residues` etc. You can find a detailed description of all ProteinDJ parameters [here](docs/parameters.md).

To launch a design campaign, launch the pipeline with `nextflow run main.nf` command from the root of the `proteindj` repository followed by the ProteinDJ parameters:

```
# Example monomer design run
nextflow run main.nf --design_mode boltzgen_denovo --num_designs 2 --seqs_per_design 1 --design_length 50 --out_dir pdj_monomer_example
```

This will launch the nextflow pipeline and show you progress in your terminal window. If you are running this over an ssh connection, you might want to use screen or tmux to avoid cancelling the process upon disconnect.

Note that by default, this will use the parameters in `nextflow.config` file in the installation folder, but you can specify a different config file. This can be useful if you want to queue different design campaigns with alternative parameters (just be sure to specify unique output directories for each):

`nextflow run main.nf -c CONFIGFILE`

A helpful feature of the nextflow.config file is profiles. Profiles can be used to override parameters but are easier to edit. e.g. the
profile for the rfd_denovo monomer mode looks like this. We recommend using the existing profiles as a reference for each mode.

```
rfd_denovo_monomer {
    params {
        design_mode = 'rfd_denovo'
        design_length = "60-80"
        seq_method = 'fampnn'
        pred_method = 'boltz'
    }
}
```

In this example, Nextflow will use all of the default parameter values from the params section except for `design_mode`, `design_length`, `seq_method` and `pred_method` (in this case to generate a de novo monomer 60-80 residues in length with RFdiffusion, Full-Atom MPNN and Boltz-2). You can use profiles by adding the `-profile` flag. e.g.

`nextflow run main.nf -profile rfd_denovo_monomer`

You can specify multiple profiles, for example, to combine a profile for your HPC environment (e.g. Milton at WEHI) with a design mode profile:

`nextflow run main.nf -profile milton,rfd_denovo_binder`

After running the pipeline, you can find all the results as well as intermediate files and logs in your specified `out_dir` organised as below:

```
out_dir/
├── configs/                # Config files used for run
├── inputs/                 # Input files used in run e.g. PDB files for binder design
├── run/                    # Intermediate results and log files with subfolders for each process
├── results/                # Final results and metadata
    ├── best_designs/       # Directory containing PDB files of designs that passed all filters
    ├── ranked_designs/     # Directory containing ranked PDB files of designs that passed all filters
    ├── all_designs.csv     # CSV file with metadata for all designs
    ├── best_designs.csv    # CSV file with metadata for designs that passed all filters
    └── ranked_designs.csv  # CSV file with metadata for ranked designs that passed all filters
└── nextflow.log            # Copy of Nextflow log from run
```

> Tip: If your run gets interrupted you can resume from the last completed step by using the -resume flag e.g. `nextflow run main.nf -profile rfd_denovo_monomer -resume`

### Tutorial - De novo Binder Design <a name="tutorial"></a>

New to ProteinDJ? Follow our step-by-step [De novo Binder Design Tutorial](tutorial/tutorial_binderdesign.md) to design a binder against the insulin receptor.

## Advanced Parameters <a name="params"></a>

We have aimed to provide as much functionality as possible of the underlying software packages and there are many parameters you can adjust. Here is a [Parameter Guide](docs/parameters.md) to all of the parameters that are configurable within the `nextflow.config` file.

## Filtering Designs <a name="params-filter"></a>

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

If a design fails to meet these criteria, BindCraft will generate a new design until it finds one that passes. This can lead to long run times compared to BoltzGen/RFdiffusion but tends to give binder designs that are more likely to succeed in the subsequent Structure Prediction stage.
We have prepared a [Filtering Guide](docs/parameters.md/#filtering-parameters) on all the filters available in ProteinDJ with recommended values for each.

## Metrics and metadata <a name="metrics"></a>

ProteinDJ generates and captures metadata for all designs in a CSV file '`all_designs.csv`' and the best designs if filtering is applied to '`best_designs.csv`'. See our [Metrics and Metadata Guide](docs/metrics.md) for a description of each metric.

## BindSweeper <a name="bindsweeper"></a>

BindSweeper is a python-based tool that can launch multiple instances of ProteinDJ to 'sweep' different binder design parameters e.g. hotspots, timesteps. For detailed information about installing and using BindSweeper, see the [BindSweeper User Guide](docs/bindsweeper.md).

<img src="img/bindsweeper_workflow.png" height="500">

## Appendices <a name="append"></a>

### Known limitations <a name="limitations"></a>

- Ligands / non-natural amino acids (e.g. PTMs) are not compatible with ProteinDJ

### Seqera Support <a name="seqera"></a>

We have designed ProteinDJ to be compatible with the [Seqera platform](https://seqera.io/platform/), so that jobs can be executed in the cloud or deployed on HPC. If you have access to Seqera, you can use the schema files in `schemas/`. There is a schema file for each mode with relevant parameters, defaults and input validation.

### Troubleshooting and common errors <a name="errors"></a>

`KeyError: 'P1L'` - A non-standard amino acid code (e.g. P1L) is present in your input PDB and included in contigs. RFdiffusion only takes natural amino acids (i.e. 'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL') and unknown or masked amino acids ('UNK','MAS').

`Unknown variable 'metadata_ch_fold'` - We are using topic channels for metadata, and this feature is only available in Nextflow v24.04 and above. This error occurs with earlier versions of Nextflow.

`FATAL:   While pulling image from oci registry: error fetching image to cache: unable to Download Image: error writing layer: stream error: stream ID 1; PROTOCOL_ERROR; received from peer` - Downloading large containers can sometimes fail from the GitHub Container Registry (GHCR), due to timeout or protocol errors. Try increasing `apptainer.pullTimeout` in nextflow.config or changing protocol e.g. `export GODEBUG=http2client=0`

### Data used for benchmarking <a name="append-bench"></a>

We used five structures for testing and benchmarking our pipeline.

| Protein           | PDB ID | Filename       | Domain boundaries                                                | Hotspots           |
| ----------------- | ------ | -------------- | ------------------------------------------------------ | ------------------ |
| Influenza A H1 HA | 5VLI   | 5vli_ha.pdb    | A4-53,A79-83,A110-114,A261-325,B501-568,B580-670 | B521,B545,B552 |
| IL-7Rα            | 3DI3   | 3di3_il7ra.pdb | B17-209                                            | B58,B80,B139   |
| InsR              | 4ZXB   | 4zxb_ir.pdb    | E6-150 |                                             E64,E88,E96    |
| PD-L1             | 5O45   | 5o45_pd-l1.pdb | A17-131                                            | A56,A115,A123  |
| TrkA              | 1WWW   | 1www_trka.pdb  | X282-382                                           | X294,X296,X333 |

### Citations for software packages used in ProteinDJ <a name="citations"></a>

ProteinDJ - Silke, D., Iskander, J., Pan, J., Thompson, A.P., Papenfuss, A.T., Lucet, I.S., Hardy, J.M. ProteinDJ: a high-performance and modular protein design pipeline. Prot Sci (2026). https://doi.org/10.1002/pro.70464

AlphaFold2 - Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2 · [Read online](https://webeditions.page/works/highly-accurate-protein-structure-prediction-with-alphafold/?via=placement:49677f46bcc6f384cb56)

AlphaFold2 Initial Guess and ProteinMPNN-FastRelax - Bennett, N.R., Coventry, B., Goreshnik, I. et al. Improving de novo protein binder design with deep learning. Nat Commun 14, 2625 (2023). https://doi.org/10.1038/s41467-023-38328-5

Arpeggia - Zhou, Y. Arpeggia: calculation of interatomic interactions in molecular structures. https://github.com/y1zhou/arpeggia

BindCraft - Pacesa, M., Nickel, L., Schellhaas, C. et al. One-shot design of functional protein binders with BindCraft. Nature 646, 483-492 (2025). https://doi.org/10.1038/s41586-025-09429-6

BioPython - Cock, P. J., Antao, T. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422-1423, (2009). https://doi.org/10.1093/bioinformatics/btp163

Boltz-2 - Wohlwend, J., et al. Boltz-2 Democratizing Biomolecular Interaction Modeling, bioRxiv 2024.11.19.624167 (2024). https://doi.org/10.1101/2024.11.19.624167

BoltzGen - Stark, H., Faltings, F., Choi, M. et al. BoltzGen: Toward Universal Binder Design, bioRxiv 2025.11.20.689494 (2025). https://www.biorxiv.org/content/10.1101/2025.11.20.689494v2

Full-Atom MPNN - Shuai, R.W., et al. Sidechain conditioning and modeling for full-atom protein sequence design with FAMPNN, bioRxiv 2025.02.13.637498 (2025). https://doi.org/10.1101/2025.02.13.637498

HyperMPNN - Ertelt, M., Schlegel, P., Beining, M. et al. HyperMPNN-A general strategy to design thermostable proteins learned from hyperthermophiles. bioRxiv (2024) https://doi.org/10.1101/2024.11.26.625397

iPSAE score scripts - Digital Biotechnology Lab. (2025). Overath, M. D., Rygaard, A., Jacobsen, C. P., Brasas, V., Morell, O., Sormanni, P., & Jenkins, T. P. (2025). Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders. bioRxiv. https://doi.org/10.1101/2025.08.14.670059v1

OpenMM and PDBFixer - Eastman, P., Swails, J., Chodera, J.D. et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. PLOS Comput Biol 13(7), e1005659 (2017). https://doi.org/10.1371/journal.pcbi.1005659

PRODIGY - Xue, L.C., Rodrigues, J.P., Kastritis, P.L., Bonvin, A.M.J.J., Vangone, A. PRODIGY: a web server for predicting the binding affinity of protein-protein complexes. Bioinformatics 32(23), 3676-3678 (2016). https://doi.org/10.1093/bioinformatics/btw514
ProteinMPNN - Dauparas, J., et al. Robust deep learning–based protein sequence design using ProteinMPNN. Science 378, 49-56 (2022). https://doi.org/10.1126/science.add2187

RFdiffusion - Watson, J.L., Juergens, D., Bennett, N.R. et al. De novo design of protein structure and function with RFdiffusion. Nature 620, 1089–1100 (2023). https://doi.org/10.1038/s41586-023-06415-8
