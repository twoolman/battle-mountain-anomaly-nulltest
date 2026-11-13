# Battle Mountain anomaly photograph: null-hypothesis tests

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24139413.svg)](https://doi.org/10.5281/zenodo.24139413)

Analysis code and source data for a formal test of whether an anomalous group of
bright pixels in a single photograph is consistent with ordinary optical and
sensor mechanisms.

The photograph was taken on 18 July 2015 at 6:42:14 PM EDT with a Samsung
Galaxy S6, on a property in Rappahannock County, Virginia.

## What the code tests

The global null hypothesis is declarative: the anomaly group of pixels is
consistent with normally expected lighting artifacts in this photograph, namely
sensor noise, lens-flare ghosting, or the bright patches that ordinary light
produces in a shaded tree line.

Three component tests are run, and the global null is rejected only if all three
reject.

| Test | Mechanism | Method |
|---|---|---|
| 1 | Sensor noise | Monte Carlo of Gaussian noise at the measured local shade level, 100,000 replicates |
| 2 | Lens-flare ghost | Distance from the ghost position predicted by the located sun disk, Rayleigh tail |
| 3 | Scene lighting | Scan statistic over every bright patch in the tree line, with a generalized Pareto tail so the look-elsewhere effect is corrected |

Test 3 runs in three strata:

- **3a** lower tree line, glare corrected, which is the anomaly's own setting
- **3b** whole tree line including sky gaps, glare corrected, which is conservative
- **3c** deep-shade band on raw luminance, comparing only blocks whose raw median
  lies in the anomaly's own shade band

Rejecting the global null supports only the statement that the pixel group is
not explained by the tested mechanisms. It does not identify what the shapes
are, and the code makes no claim about that.

## Method notes

**Veiling glare correction.** The sun is in frame, and its glare adds a smooth
additive haze that lifts the darkest shadows from roughly 5 far from the sun to
roughly 65 near it on the 0-255 scale. The glare floor is estimated as a low
quantile of luminance within each block, smoothed, and subtracted before patches
are measured. Without this step only the darkest corner of the frame is
comparable with the anomaly.

**NIEL** is the normalized integrated excess luminance: the sum over a patch of
the glare-corrected luminance divided by the local scene level, minus one. It is
the patch's total excess brightness in units of its own local scene level, which
makes it invariant to illumination level and responsive to both brightness and
size.

**Look-elsewhere correction.** The anomaly was found by looking at the
photograph, so it is compared against the maximum of N ordinary patches rather
than against a typical one. A generalized Pareto distribution is fitted above a
quantile of the null statistic and the probability is computed for the maximum
of N draws. A bootstrap of the tail exceedances gives the uncertainty of that
p-value, which is wide when the tail holds few patches.

**Alpha.** 0.001, Bonferroni-divided across the three component tests, so each
test is evaluated at 0.000333.

**Sensitivity.** Test 3 is re-run across three detection thresholds and three
tail quantiles in each of the three strata, and the share of settings in which
the null is rejected is reported per stratum.

## Requirements

R, with:

- `jpeg`
- `evd`
- `spatstat.geom`

The script installs any that are missing.

## Input data

The source photograph is included in this repository as
`20150718_184215.jpg`, so the analysis is reproducible as distributed.

It is the unmodified file as received, renamed only so that the filename the
script expects matches the file on disk. No resampling, cropping, color
correction or recompression has been applied, because every test in this
analysis measures luminance and any of those operations would alter the
quantity under test.

### Photograph provenance

The photograph was taken by the author, on land the author owned at the time,
using the author's own phone. The original file was subsequently lost from that
handset and was recovered from a copy the author had sent to a family member,
who returned the unaltered high-resolution original by email. The file
retained in this project under the name
`20150718_184215_original_from_email.jpg` is that recovered original, and the
copy distributed here as `20150718_184215.jpg` is identical to it apart from
the filename.

This history is recorded because it explains the original filename and because
it is the whole of the chain of custody. The file has at no point passed
through an editor, an uploader that recompresses, or any processing step.

### Photograph rights

Copyright in the photograph is held by the author, Thomas A. Woolman.

The MIT license in `LICENSE` covers the source code only and does not extend
to the photograph. The photograph is made available under the Creative Commons
Attribution 4.0 International license, CC BY 4.0, which permits reuse with
attribution to the author.

## Running

```r
source("bm_null_test_v3.R")
```

The script is written to run from the folder holding the photograph, which in
this repository is the repository root. The random seed is fixed at 20150718,
so a run reproduces exactly.

## Outputs

All written to `results/`:

- `bm_null_test_<timestamp>.txt`, the full run log
- `patches_lower_treeline.csv`, `patches_whole_treeline.csv`, `patches_deep_shade_raw.csv`
- `sensitivity_test3.csv`
- `fig_null_distribution.png`
- `fig_patch_map.png`
- `sessionInfo.txt`, written on every run so that the R and package versions
  quoted in any report are read from a file the run itself produced rather than
  reconstructed from memory

`results/` is excluded from version control, since every file in it is
regenerated by any run.

## Versions

- `bm_null_test_v3.R` is current. It adds the `sessionInfo.txt` record.
- `bm_null_test_v2.R` and `bm_null_test.R` are retained for provenance.

## Status

This repository holds the analysis code and the source photograph. It contains
no manuscript text and no results.

Release v3.0.0 is archived on Zenodo at
[10.5281/zenodo.24139413](https://doi.org/10.5281/zenodo.24139413).

The manuscript reporting these results is posted as a preprint on Research
Square:

> Woolman, T. A. (2026). *Testing an Anomalous Photograph Against Ordinary
> Image-Forming Mechanisms: A Null Hypothesis Framework and a Single-Case
> Demonstration from Battle Mountain, Virginia, USA* [Preprint]. Research
> Square. https://doi.org/10.21203/rs.3.rs-11707102/v1

The preprint is not peer reviewed.

## How to cite

Cite the archived release rather than the repository, so that the version used
is identifiable:

> Woolman, T. A. (2026). *Analysis code for null-hypothesis tests of an
> anomalous pixel group in a single photograph* [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.24139413

That is the **concept DOI**, which resolves to all versions and always
redirects to the most recent one. Use it when citing the software in general.

To cite the exact code that produced a particular result, use the **version
DOI** instead, which for v3.0.0 is `10.5281/zenodo.24139414`.

`CITATION.cff` carries the same metadata in machine-readable form.

## Author

Thomas A. Woolman, PhD
On Target Technologies, Inc., Ridgeley, West Virginia, United States
ORCID: [0000-0002-4120-2157](https://orcid.org/0000-0002-4120-2157)

## License

Source code: MIT. See `LICENSE`.
Photograph: CC BY 4.0, copyright Thomas A. Woolman. See **Photograph rights**
under **Input data** above.
