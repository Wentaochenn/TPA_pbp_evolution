# TPA PBP evolution

Analysis code for **Recurrent evolution of penicillin-binding proteins across global Treponema pallidum populations**.

## Reproduce the published figures and tables

Download and extract the companion Zenodo data archive. The input directory is `Data_for_Zenodo`; use a new output directory.

```bash
bash scripts/setup_environment.sh
conda activate pbp_figures_release
conda env create -f scripts/structure_environment.yml
bash master.sh reproduce --data /path/to/Data_for_Zenodo \
  --out /path/to/reproduced_outputs --font-dir /path/to/Arial
```

The R workflow regenerates five main figures, four supplementary figures, three supplementary tables and the data used by every figure panel. Inputs are checked against their archived SHA-256 values. Figures use Arial regular, bold, italic and bold italic; provide genuine font files with `--font-dir`, or omit that option if these four faces are installed. PDF fonts and editable SVG text are checked after export. Arial font files are supplied separately and are not distributed in this repository.

Figure 5 shows the existing common-ancestor node-date 95% HPD intervals and thin grey guides from sampled tips to their variant-state tracks. The node intervals and guide coordinates are included in the corresponding Source Data.

For a strict comparison with archived figures, use the code commit recorded in `Data_for_Zenodo/RELEASE.json` and add `--compare`. `REPRODUCTION_CHECK.json` then reports numerical and graphical agreement. The current code can also use earlier frozen analytical inputs to generate the updated figure styling.

| Directory | Contents |
|---|---|
| `Figures/` | PDF, SVG and PNG figures |
| `Figure_Source_Data/` | Data organised by figure and panel |
| `Supplementary_Tables/` | Genome metadata and amino-acid states; 40 observed substitutions with ESM-2 scores; FUBAR topology sensitivity |
| `Reproducibility_Data/` | Supporting analytical tables |

Supplementary Table 1 contains 1,861 genomes, 42 surveyed amino-acid positions across five genes and `mrcA_protein_type`. Codon states that cannot be resolved are marked `uncallable`; unresolved protein types are marked `untyped`. The workbook includes the callability criteria.

Reference-structure images are checked against the archived model, rendering script and camera settings. Compatible renderings are reused; otherwise the reference PDB is rendered using the specified PyMOL environment and fixed camera. Structure labels are editable Arial text. Add `--rerender-structure` to request a fresh rendering, or use `--pymol-python /path/to/python` for another PyMOL interpreter.

## Reanalyse the study inputs

Upstream analysis scripts are in `scripts/upstream/`, with their environment specifications in `scripts/upstream/envs/`. Additional analysis modules are in `scripts/analysis/`. Python requirements are listed in `scripts/analysis/requirements.txt`.

Prepare the deposited input files and retrieve paired reads:

```bash
bash master.sh inputs prepare --data /path/to/Data_for_Zenodo \
  --out /path/to/source_data
bash master.sh inputs download-reads --source /path/to/source_data
bash master.sh upstream prepare --source-project /path/to/source_data \
  --run-root /path/to/new_analysis
bash master.sh upstream run check --project /path/to/new_analysis
bash master.sh upstream run STAGE --project /path/to/new_analysis
```

The downloader uses the 1,856 study run accessions and verifies paired-end layout and ENA MD5 checksums. The five assembly genomes, reference FASTAs, metadata and sampling dates are supplied in `Upstream_Inputs/`. Downloads can be limited with `--runs RUN1 RUN2`.

| Stage | Analysis |
|---|---|
| `wgs`, `p0`, `p0_02` | Read processing, consensus alignment and assembly variant calls |
| `p1` | Genome-wide phylogeny |
| `p2` | Codon-aware variant extraction and callability |
| `p3` | Ancestral reconstruction, homoplasy and ER/ARD stochastic mapping |
| `selection` | FUBAR site-specific selection |
| `temporal`, `treetime` | Temporal-signal assessment |
| `beast_prepare`, `beast_run`, `pbp_temporal` | Clock inputs, BEAST2 and acquisition dates |

ER/ARD fitting includes a numerical fallback if the primary optimizer fails, followed by stochastic mapping. BEAST2 reruns use the deposited study XML files. The raw-read workflow uses Snakemake and SGE; adjust computing resources to the target system.

The `pbp_temporal` stage uses the same two-clock implementation as `analysis clock-histories`. Copy `scripts/analysis/clock_histories/config.example.json` into the analysis project as `config/clock_histories.json` and set the strict and relaxed XML, trace, posterior-tree and MCC paths before running this stage. `PBP_CLOCK_HISTORY_CONFIG` selects another configuration file.

Configure external references with `PBP_HOST_REF_HUMAN`, `PBP_HOST_REF_RABBIT`, `PBP_TRIMMOMATIC_ADAPTERS` and `PBP_STRAINSEEKER_DB`. Raw-read annotation uses the `TPA_SS14` and `TPA_Nichols` snpEff databases. Database identities, access routes and available files are recorded in `Upstream_Inputs/resource_provenance.json`. Figure reproduction uses the frozen analytical inputs and does not require raw reads or these databases.

Additional analyses use `bash master.sh analysis MODULE -- ARGUMENTS`; configuration templates accompany each module. Run `bash master.sh analysis MODULE -- --help` for its arguments.

| Module | Analysis |
|---|---|
| `clustered-dates` | Clustered date randomisation and subset assessment |
| `selection-topologies` | FUBAR across genome-wide bootstrap topologies |
| `protein-types` | MrcA sequence extraction and protein-type grouping |
| `clock-histories` | Clock-chain summaries and posterior acquisition dates |

## ESM-2 sequence-context scores

In a GPU environment with Biopython, PyTorch, transformers and pandas, run:

```bash
bash master.sh esm2 -- recompute --project /path/to/new_analysis \
  --model /path/to/esm2_t36_3B_UR50D --out /path/to/esm2_results
```

Reference proteins are translated from complete Nichols CDS features using bacterial translation table 11 and checked against the annotated GenBank proteins, including initiating methionine. The unchanged scoring implementation evaluates 19 alternate residues at each of 2,594 positions in four proteins, yielding 49,286 scores. Supplementary Table 2 contains the 40 observed substitutions.

To compare with an existing complete score table, add `--previous-scores /path/to/scores.tsv`. To regenerate figures using a new scan, add `--esm2-results /path/to/esm2_results/inference` to `master.sh reproduce`.

Model weights are retrieved separately:

```bash
bash master.sh inputs download-model --spec /path/to/source_data/esm2_model.json \
  --out /path/to/esm2_t36_3B_UR50D
```

An optional `--revision` selects a model repository revision. Model and tokenizer files must match the deposited SHA-256 values.

Run ESM-2 within a GPU allocation provided by the local scheduler. Reference translation, scoring and figure input binding use the same validated output structure. Comparing with previous scores is optional.

## Package a completed analysis

```bash
bash master.sh inputs package --project /path/to/completed_analysis \
  --data /path/to/Data_for_Zenodo --out /path/to/new_package \
  --resources /path/to/resource_provenance.json
```

This writes new code and data archives and checksums. Existing results are retained. Raw reads, database files and model weights are retrieved separately. The deposited resource metadata provide their identities and access information; these metadata are reused if `--resources` is omitted.

## Citation

Please cite the accompanying study and its final archived code and data records. Software citation metadata are supplied in `CITATION.cff`.
