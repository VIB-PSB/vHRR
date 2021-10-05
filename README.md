# vHRR: co-expression based gene function prediction

## Requirements

The pipelines were tested with Python3.6.5, Nextflow v20.10.0 and

```
numpy==1.16.1
pandas==1.1.5
```


## Installation

1. Clone repo
2. Add path to `go.obo` file in `go_config.sh`
3. Run `bash go_config.sh`

## Run

The folders `input` and `nextflow_config` currently contain example data. Refer to this to check data formats.

1. Replace the example data `./input/` with your own.
2. Run `bash preprocess_go_data.sh` (replace `go_gene.tsv` with the filename of your input GO annotation file).
3. Replace the input files in `nextflow_config/vHRR_single.config` or `nextflow_config/vHRR_multiple.config` with your filenames and adapt the nextflow parameters to  fit your needs and system.
4. For a single atlas: run `bash vHRR_single.sh`. For multiple atlases: run `bash vHRR_multiple.sh`. 

## Citation

If you publish results generated using vHRR, please cite:

Depuydt, T. and Vandepoele, K. (2021), Multi-omics network-based functional annotation of unknown Arabidopsis genes. Plant J. Accepted Author Manuscript. https://doi.org/10.1111/tpj.15507
