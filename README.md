# Principal component (PC) projector for ancestry

PC-projector is an R script that calculates reference-projected principal components (PCs) based on genetic data. Projection of PCs is is useful when you want the PC to be independent of the target sample. The R-script does this using PLINK v1.9. First using the reference genetic data provided, a list of SNPs that are LD independent and suitable for calculating ancestry is made. Then PCs based on the reference sample are projected into the target sample. The PCs are then scaled to the reference sample (i.e. the reference mean is subtracted from each score, and divided by the reference standard deviation).

## Getting started

### Prerequisites

* R and the required packages:

```R
install.packages(c('data.table','optparse','foreach','doMC'))
```

* [PLINK 1.9 software](https://www.cog-genomics.org/plink2)

* Reference genetic data
  * Binary PLINK format (.bed/.bim/.fam)

* Target sample genetic data
  * Binary PLINK format (.bed/.bim/.fam) 
  * RSIDs should match the reference data
  * The data should not contain any duplicate RSIDs

### Parameters

| Flag                | Description                                                  | Default |
| :------------------ | ------------------------------------------------------------ | :-----: |
| --target_plink_chr | Path to per chromosome target PLINK files (.bed/.bim/.fam) | NA |
| --ref_plink_chr | Path to per chromosome reference PLINK files (.bed/.bim/.fam) | NA |
| --plink | Path PLINK software binary | plink |
| --output | Name of output files. Set it to be in an empty directory. | ./PC_projector_output/Output |
| --batch_size | Number of individuals in each batch | 5000 |
| --n_pcs | Number of PCs to extract | 20 |

### Output files

In the specified output directory, the following files will be produced:

| Name                   | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| .eigenvec | File containing projected PC scores for the target sample |
| .log | Log file. |

## Examples

These examples use FUSION 1KG reference data ([download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)).

##### When using default settings:

```shell
Rscript ./PC_projector.R \
--target_plink_chr ./FUSION/LDREF/1000G.EUR. \
--ref_plink_chr ./FUSION/LDREF/1000G.EUR. \
--plink ./Software/plink1.9 \
--output ./PC_projector_demo/PC_projections \
--n_pcs 100
```
## Help

This script was written by Dr Oliver Pain.

If you have any questions or comments please contact me at oliver.pain@kcl.ac.uk







