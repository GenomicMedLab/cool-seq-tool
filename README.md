<h1 align="center">
Cool-Seq-Tool
</h1>

[![image](https://img.shields.io/pypi/v/cool-seq-tool.svg)](https://pypi.python.org/pypi/cool-seq-tool) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14007783.svg)](https://doi.org/10.5281/zenodo.14007783) [![image](https://img.shields.io/pypi/l/cool-seq-tool.svg)](https://pypi.python.org/pypi/cool-seq-tool) [![image](https://img.shields.io/pypi/pyversions/cool-seq-tool.svg)](https://pypi.python.org/pypi/cool-seq-tool) [![Actions status](https://github.com/genomicmedlab/cool-seq-tool/actions/workflows/checks.yaml/badge.svg)](https://github.com/genomicmedlab/cool-seq-tool/actions/checks.yaml)

---

**[Documentation](https://coolseqtool.readthedocs.io/stable/)** · [Installation](https://coolseqtool.readthedocs.io/stable/install.html) · [Usage](https://coolseqtool.readthedocs.io/stable/usage.html) · [API reference](https://coolseqtool.readthedocs.io/stable/reference/index.html)

---

## Overview

<!-- description -->
The Common Operations On Lots-Of Sequences Tool, **Cool-Seq-Tool**, provides:

- A Pythonic API on top of sequence data of interest to tertiary analysis tools, including mappings between gene names and transcripts, [MANE transcript](https://www.ncbi.nlm.nih.gov/refseq/MANE/) descriptions, and transcript alignment data from the [Universal Transcript Archive](https://github.com/biocommons/uta)
- Augmented access to the [SeqRepo](https://github.com/biocommons/biocommons.seqrepo) database, including multiple additional methods and tools
- Mapping tools, including a transcript selection algorithm for selecting a representative transcript defined [here](https://coolseqtool.readthedocs.io/stable/transcript_selection.html), that combine the above to support translation between references sequences, annotation layers, and transcripts
<!-- /description -->

---

## Install

Cool-Seq-Tool is available on [PyPI](https://pypi.org/project/cool-seq-tool)

```shell
python3 -m pip install cool-seq-tool
```

See the [installation instructions](https://coolseqtool.readthedocs.io/stable/install.html) in the documentation for a description of dependency setup requirements.

---

## Usage

All Cool-Seq-Tool resources can be initialized by way of a top-level class instance:

```pycon
>>> from cool_seq_tool import CoolSeqTool
>>> from cool_seq_tool.schemas import AnnotationLayer, CoordinateType
>>> cst = CoolSeqTool()
>>> result = await cst.mane_transcript.get_mane_transcript(
...     "NP_004324.2",
...     599,
...     AnnotationLayer.PROTEIN,
...     coordinate_type=CoordinateType.INTER_RESIDUE,
... )
>>> result.gene, result.refseq, result.status
('EGFR', 'NM_005228.5', <TranscriptPriority.MANE_SELECT: 'mane_select'>)
```

---

## Feedback and contributing

We welcome bug reports, feature requests, and code contributions from users and interested collaborators. The [documentation](https://coolseqtool.readthedocs.io/stable/contributing.html) contains guidance for submitting feedback and contributing new code.
