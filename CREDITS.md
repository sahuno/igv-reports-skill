# Credits and Attribution

This project — `igv-reports-skill` — is a driver and verifier layer built
on top of the upstream `igv-reports` package. We did not write the HTML
rendering engine. This file documents what we reuse, who wrote it, and
how to cite it.

## Upstream project

**igv-reports** — <https://github.com/igvteam/igv-reports>

Maintained by the IGV team at the Broad Institute (Jim Robinson, Helga
Thorvaldsdóttir, and contributors). Provides the `create_report` CLI
that slices BAM/VCF/track data for each genomic region of interest and
embeds the slices into a self-contained HTML page driven by `igv.js`.

This skill pins `igv-reports >= 1.16.0` as a runtime dependency
(`pyproject.toml`); the upstream package is installed from PyPI via
`pip install igv-reports` and is **not vendored** into this repository.

### Upstream license (verbatim)

The following text is reproduced verbatim from
<https://github.com/igvteam/igv-reports/blob/master/LICENSE>:

```
The MIT License (MIT)

Copyright (c) 2018-2019 The Broad Institute and The Regents of the University of California

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:


The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.


THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

## What is reused vs. what this repo adds

| Component | Source | What it does |
|---|---|---|
| `create_report` CLI + the HTML rendering | upstream `igv-reports` (MIT, © Broad / UC Regents) | Slices BAM/VCF/track data per region; emits a single self-contained HTML driven by igv.js. **This is the engine.** |
| `igv.js` (embedded in the rendered HTML) | upstream IGV project (MIT, © Broad / UC Regents) | In-browser genome viewer that powers the rendered HTML. Transitively bundled by `igv-reports`. |
| `scripts/build_igvreports.py` | this repo (MIT, © Samuel Ahuno) | Cohort-aware driver: samplesheet mode, default-track auto-resolution per genome, SIF auto-detect, run logging, opt-in anchor capture. |
| `scripts/verify_report.py` | this repo | Per-sample structural verifier: parses the HTML's embedded `tableJson` + `sessionDictionary` and confirms region count / coordinates / track names match inputs. |
| `scripts/verify_cohort.py` | this repo | Cohort-level structural verifier: cross-sample contamination scan, sample-id consistency. |
| `scripts/verify_anchors.py` | this repo | Content verifier: re-counts embedded BAM slices with `samtools` and compares to anchors frozen at build time. |
| `scripts/generate_tracks_json.py` | this repo | YAML → tracks.json builder with ONT-methylation defaults (colorBy=basemod2, min:0/max:100, Okabe-Ito pairing). |
| `scripts/prep_track.sh` | this repo | gunzip → sort → bgzip → tabix utility for annotation tracks. |
| `SKILL.md`, `references/`, `tests/`, `examples/` | this repo | Documentation, eval harness, and reference invocations. |

Bug reports against the HTML rendering itself (igv.js viewer behaviour,
data-slicing edge cases, BED/VCF parsing inside `create_report`) belong
upstream at <https://github.com/igvteam/igv-reports/issues>. Issues with
the driver, the verifiers, the ONT methylation presets, or this repo's
documentation belong at
<https://github.com/sahuno/igv-reports-skill/issues>.

## Citation

If you use the rendered HTML reports in a publication, please cite the
igv.js paper (the HTML's embedded viewer):

> Robinson, J. T., Thorvaldsdóttir, H., Turner, D., & Mesirov, J. P.
> (2023). igv.js: an embeddable JavaScript implementation of the
> Integrative Genomics Viewer (IGV). *Bioinformatics*, 39(1), btac830.
> <https://doi.org/10.1093/bioinformatics/btac830>

BibTeX:

```bibtex
@article{robinson2023igvjs,
  author    = {Robinson, James T. and Thorvaldsd{\'o}ttir, Helga and Turner, Douglass and Mesirov, Jill P.},
  title     = {igv.js: an embeddable {JavaScript} implementation of the {Integrative Genomics Viewer} ({IGV})},
  journal   = {Bioinformatics},
  volume    = {39},
  number    = {1},
  pages     = {btac830},
  year      = {2023},
  doi       = {10.1093/bioinformatics/btac830},
  publisher = {Oxford University Press}
}
```

If you used this repo's cohort driver, verifiers, or ONT methylation
presets, please also cite this repository by URL:

> Ahuno, S. (2026). *igv-reports-skill: cohort-aware driver and
> verifiers for igv-reports.* <https://github.com/sahuno/igv-reports-skill>

(A Zenodo DOI may be assigned to a future tagged release; until then the
repository URL is the citable record.)

## Acknowledgments

Thank you to the IGV team at the Broad Institute — Jim Robinson, Helga
Thorvaldsdóttir, Douglass Turner, Jill Mesirov, and contributors — for
maintaining `igv-reports` and `igv.js` as open-source MIT-licensed
projects. This driver would not exist without their engine.

## Trademarks

"IGV" and "Integrative Genomics Viewer" are names associated with the
Broad Institute's IGV project. This repository is an independent driver
built on top of `igv-reports` and is **not** an official IGV team product.
