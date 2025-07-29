# Data Repository

`misosoup`: A metabolic modeling tool for identifying minimal microbial communities reveals pervasive cross-feeding–driven niche expansion

## Overview

- `data/` any necessary input data
- `out/` stores all computed results
- `data/strains/` contains all reconstructed marine microbial models
- `data/media.yaml` defines the media compositions for all environments
- `config.yaml` specifies the configurations to be tested
- `Snakefile` defines the processing pipeline for minimal community search

## Usage

Install [`uv`](https://github.com/astral-sh/uv):

```
pipx install uv
```

Run minimal community search with `snakemake` and `misosoup`:

```
uv run snakemake --cores 1
```

Run growth analysis in isolation with `cobra`:

```
uv run scripts/cobra_growth.py
```
