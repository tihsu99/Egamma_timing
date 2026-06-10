# Egamma timing study

This repository contains the HGCAL timing comparison workflow for online and offline egamma objects, with a dedicated `v4` vs `v5` comparison pipeline.

The comparison pipeline does four things:
1. rerun the HLT side products in process `HLTX`
2. rerun offline `RAW2DIGI,L1Reco,RECO,RECOSIM,PAT`
3. build comparison ntuples for `online_v4`, `online_v5`, `offline_v4`, and `offline_v5`
4. convert the merged ROOT outputs into parquet tables and summary plots

## Layout

- `EgammaTimingTools/plugins`
  C++ ntuplizers and helpers
- `EgammaTimingTools/python`
  CMSSW configs and submission scripts
- `EgammaTimingTools/config`
  workflow and sample definitions
- `EgammaTimingTools/analysis`
  merge, parquet production, and plotting scripts

## Important comparison choices

- Both `v4` and `v5` use `ticlTrackstersMerge`
- `v5` timing is propagated back to the origin before comparison
- The workflow config that controls this is [comparison_workflow_example.yaml](/Users/tihsu/PycharmProjects/Egamma_timing/EgammaTimingTools/config/comparison_workflow_example.yaml)

## Setup

The production step should be run inside a CMSSW area.

```bash
cmsrel CMSSW_15_1_0_pre6
cd CMSSW_15_1_0_pre6/src
cmsenv
git clone https://github.com/tihsu99/Egamma_timing.git
scram b -j 8
cd Egamma_timing
```

If you already have a suitable CMSSW release, use that instead, but keep `cmsenv` active for the production step.

## Step 1: Prepare the workflow config

Start from:

```text
EgammaTimingTools/config/comparison_workflow_example.yaml
EgammaTimingTools/config/comparison_samples_example.yaml
```

What to edit:

- `cmssw-dir`
  Optional. Leave empty to auto-detect the current CMSSW area.
- `timing-dir`
  Usually keep as `..`
- `sample-config`
  Path to your sample yaml
- `samples`
  Add your datasets or local/remote ROOT paths in `comparison_samples_example.yaml`

The sample file supports both:

```yaml
dataset: /A/B/C
dbs_instance: prod/global
```

and

```yaml
path: /full/path/or/xrootd/location
```

## Step 2: Produce the comparison ntuples

This step runs the full rerun chain and writes one comparison ntuple per input file.

Example:

```bash
cd /path/to/CMSSW_X_Y_Z/src/Egamma_timing/EgammaTimingTools/python

mkdir -p $HOME/tmp
voms-proxy-init -voms cms -out $HOME/tmp/x509up
export X509_USER_PROXY=$HOME/tmp/x509up

python3 submit_condor_comparison.py \
  --config ../config/comparison_workflow_example.yaml \
  --sample-config ../config/comparison_samples_example.yaml \
  --mode produce \
  --variant v4 v5 \
  --outdir /path/to/output/comparison_ntuple \
  --particle electron \
  --region endcap
```

Useful filters: (Note that `max-events` > 60 will likely cause the job to exceed the 2 GB memory limit on Condor, normally `workday` or `tomorrow` is a safer option for large jobs)

- `--particle electron` or `--particle photon`
- `--region endcap`
- `--sample DY_wTICL`
- `--max-files 2`
- `--max-events 50`
- `--test`
  only generate the shell scripts, do not submit

Output layout:

```text
/path/to/output/comparison_ntuple/
  v4/<particle>/<region>/<sample>/comparisonNtuple_*.root
  v5/<particle>/<region>/<sample>/comparisonNtuple_*.root
```

## Step 3: Merge v4 and v5 ntuples

This step matches the `v4` and `v5` ntuples file by file and writes one merged ROOT file containing:

- `gen_oriented`
- `offline_v4`
- `offline_v5`
- `online_v4`
- `online_v5`

Example:

```bash
cd /path/to/CMSSW_X_Y_Z/src/Egamma_timing/EgammaTimingTools/python

python3 submit_condor_comparison.py \
  --config ../config/comparison_workflow_example.yaml \
  --mode merge \
  --v4_dir /path/to/output/comparison_ntuple/v4 \
  --v5_dir /path/to/output/comparison_ntuple/v5 \
  --outdir /path/to/output/comparison_merged \
  --particle electron \
  --region endcap
```

Output layout:

```text
/path/to/output/comparison_merged/
  <particle>/<region>/<sample>/comparisonNtuple_*.root
```

## Step 4: Build parquet tables for isolation and timing studies

Run the analysis step outside `cmsenv`. The repository already provides an LCG environment helper:

```bash
cd /path/to/CMSSW_X_Y_Z/src/Egamma_timing/EgammaTimingTools/analysis
source script/env.sh
```

Example:

```bash
python3 analysis_comparison_hgcal.py \
  --inputDir /path/to/output/comparison_merged \
  --outDir /path/to/output/comparison_iso \
  --particle electron \
  --region endcap \
  --signal DY_wTICL \
  --background TT_wTICL \
  --n-workers 4
```

Output layout:

```text
/path/to/output/comparison_iso/
  <particle>/<region>/online_v4/<process>_iso.parquet
  <particle>/<region>/online_v5/<process>_iso.parquet
  <particle>/<region>/offline_v4/<process>_iso.parquet
  <particle>/<region>/offline_v5/<process>_iso.parquet
```

## Step 5: Make summary plots

Example:

```bash
python3 summary_comparison_hgcal.py \
  --isoDir /path/to/output/comparison_iso \
  --mergedDir /path/to/output/comparison_merged \
  --plotDir /path/to/output/comparison_plots \
  --particle electron \
  --region endcap \
  --signal DY_wTICL \
  --background TT_wTICL \
  --n-workers 4
```

This produces ROC curves, AUC scans, and timing distributions for:

- `online_v4`
- `online_v5`
- `offline_v4`
- `offline_v5`

## Minimal end-to-end example

```bash
cd /path/to/CMSSW_X_Y_Z/src/Egamma_timing/EgammaTimingTools/python

python3 submit_condor_comparison.py \
  --config ../config/comparison_workflow_example.yaml \
  --sample-config ../config/comparison_samples_example.yaml \
  --mode produce \
  --variant v4 v5 \
  --outdir /tmp/comparison_ntuple \
  --particle electron \
  --region endcap \
  --sample DY_wTICL \
  --max-files 1 \
  --max-events 20 \
  --test

python3 submit_condor_comparison.py \
  --config ../config/comparison_workflow_example.yaml \
  --mode merge \
  --v4_dir /tmp/comparison_ntuple/v4 \
  --v5_dir /tmp/comparison_ntuple/v5 \
  --outdir /tmp/comparison_merged \
  --particle electron \
  --region endcap

cd ../analysis
source script/env.sh

python3 analysis_comparison_hgcal.py \
  --inputDir /tmp/comparison_merged \
  --outDir /tmp/comparison_iso \
  --particle electron \
  --region endcap \
  --signal DY_wTICL \
  --background TT_wTICL

python3 summary_comparison_hgcal.py \
  --isoDir /tmp/comparison_iso \
  --mergedDir /tmp/comparison_merged \
  --plotDir /tmp/comparison_plots \
  --particle electron \
  --region endcap \
  --signal DY_wTICL \
  --background TT_wTICL
```

## Notes

- `produce` mode uses Condor submission by default. Use `--test` first if you want to inspect the generated shell scripts.
- `merge` mode currently runs locally.
- If `uproot` complains about broken baskets after `hadd`, use the original per-process merged files instead of a second-level hadd output.
- The analysis scripts expect the merged ROOT format written by `merge_comparison_trees.py`.
