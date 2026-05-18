# PFM-Depression

This repository is an archived MATLAB codebase for delineating functional brain
networks entirely within individual subjects, a procedure called precision
functional mapping (PFM; see Gordon et al. 2017 Neuron).

The original tutorial remains in `PFM-Tutorial/pfm_tutorial.m`. It calls the
legacy MATLAB utilities in `PFM-Tutorial/Utilities` and can be run with the
example PFM-style dataset on OpenNeuro:
https://openneuro.org/datasets/ds005118/versions/1.0.0. The original
step-by-step instructions are in `PFM-Tutorial/pfm_tutorial_instructions.pdf`.

## Archive status

This repository is preserved for reproducibility and for users who want the
original MATLAB workflow. New PFM development has moved into the
`pfm-mefmri` release of the multi-echo fMRI pipeline:

`/home/charleslynch/MultiEchofMRI-Pipeline-Beta/release/pfm-mefmri`

The newer pipeline implements both strategies below without MATLAB. The MATLAB
versions in this repository are provided for users who want to stay close to the
original tutorial structure.

## Newer PFM strategies mirrored here

### 1. Infomap plus consensus labeling

Infomap plus consensus labeling remains the gold-standard and preferred
approach for high-quality, densely sampled PFM datasets. It is intentionally
close to the original PFM workflow: identify subject-specific Infomap
communities, assign those communities to canonical networks, review the results,
and adjust assignments when needed.

The updated workflow mainly changes the algorithmic details of assigning
community identities. Instead of labeling one Infomap density column, it labels
communities across multiple density columns, scores each community against
canonical FC and spatial priors, writes review tables, and produces consensus
outputs across thresholds. It also writes probability maps that represent
weighted network support, which can be useful because winner-take-all labels can
be ambiguous when there is evidence supporting more than one plausible network
assignment.

Manual review remains essential. The automated assignments are a starting point,
not a substitute for inspecting the maps, reviewing the community tables, and
changing labels when the subject-level organization supports a different
interpretation.

The MATLAB port is:

`PFM-Tutorial/Utilities/pfm_infomap_consensus_labeling.m`

Primary outputs:

- `<OutFile>_DensityXX.dlabel.nii`: canonical labels for each Infomap density
- `<OutFile>_DensityXX_CommunityTable.csv`: community-level scoring table
- `<OutFile>_ModeConsensus.dlabel.nii`: mode consensus across densities
- `<OutFile>_ProbabilityConsensus.dtseries.nii`: per-network consensus support

This strategy keeps subject-specific Infomap communities, but makes the final
network identity more stable by considering all requested graph densities rather
than relying on one selected column.

### 2. Ridge fusion

Ridge fusion directly assigns each grayordinate to a canonical network by
combining two sources of evidence:

- a ridge-regression fit between the vertex's functional connectivity
  fingerprint and canonical network FC priors
- a spatial-prior term that favors anatomically plausible assignments

The main advantage is practical robustness and speed. Ridge fusion may be more
appropriate for noisier datasets or applications where fully manual,
Infomap-centered PFM is not feasible. It incorporates priors into vertex-level
network assignments in a Bayesian-like manner, the weight of those priors can be
adjusted as needed, and it is much faster than Infomap while requiring far less
RAM. It can be sufficient for many applications.

For high-quality, densely sampled PFM data, ridge fusion is less suitable than
Infomap plus consensus labeling because it moves directly to vertex-level
canonical labels rather than preserving subject-specific community structure as
the primary object of review.

The MATLAB port is:

`PFM-Tutorial/Utilities/pfm_ridge_fusion.m`

Primary outputs:

- `<OutFile>.dlabel.nii`: winner-take-all ridge-fusion network labels
- `<OutFile>_ProbMaps.dtseries.nii`: per-network assignment probabilities
- `<OutFile>_R2.dtseries.nii`: approximate FC fingerprint fit quality
- `<OutFile>.L.border` and `<OutFile>.R.border`: Workbench border files

## Example use

A compact example that starts from an already concatenated and smoothed CIFTI is
provided here:

`PFM-Tutorial/pfm_mefmri_modern_example.m`

It shows both newer paths:

1. run classic Infomap across graph densities
2. label every Infomap density and write mode/probability consensus maps
3. run ridge fusion from a concatenated CIFTI, distance matrix, and priors

The example assumes the same core MATLAB dependencies as the original tutorial:
the FieldTrip-style CIFTI helpers from the MSC codebase
(`ft_read_cifti_mod`, `ft_write_cifti_mod`), Connectome Workbench, Infomap, and
the GIfTI MATLAB reader. We previously forgot to call out that the
`ft_read_cifti_mod` and `ft_write_cifti_mod` functions are needed from the
Midnight Scan Club codebase:
https://github.com/MidnightScanClub/MSCcodebase. Sorry about that.

Contact Chuck at cjl2007@med.cornell.edu if questions arise.



