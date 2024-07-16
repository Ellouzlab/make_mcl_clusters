# INSTALLATION

## Method 1:
1. Create a new environment:

```mamba create -n mobile_clusters python=3.10```


2. Navigate into package. Install via pip.

```python -m pip install -e .```

3. Install dependencies

```mamba install -c bioconda -c conda-forge diamond=2.0.0 mmseqs2=15.6 mcl=22-282```

NOTE: mamba can be substituted for conda if you prefer; however, mamba is much faster.

## Method 2

```mamba install -c conda-forge mobile_clusters```

(NOT YET AVAILABLE)