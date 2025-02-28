# Podarcis adaptation

Data and code repository of the individual-level metagenomic data of Podarcis muralis and Podarcis liolepis lizards throughout a captivity experiment with microbiota manipulation. The project aims to understand the role of gut microbes in the thermal adaptation to high temperatures using a hologenomic approach.


## 1. Study design

The samples correspond to the first experiment performed the summer of 2023. Samples correspond to both Podarcis muralis and Podarcis liolepis, collected in different sampling points: Wild (0), After acclimation (1), After antibiotics (2), Transplant 1 (3), Transplant 2 (4), Post FMT 1 or first sample taken after the fecal microbiota transplant (5) and Post FMT 2 or second sample taken after the fecal microbiota transplant (6). P.muralis individuals were catch in Aiako Harria (a.k.a North) and P.liolepis in Lizarra (a.k.a South). First location corresponding to the cold and wet environment and the second to the hot and dry one.

## 2. Bioinformatic procedures

Data processing to generate annotated metagenome-assembled genomes and genome count tables was conducted using the following Snakemake pipeline: [mg_assembly](https://github.com/3d-omics/mg_assembly). Data analysis procedures source from the outputs of this pipeline.

## 3. Data analysis

Samples from both species were analysed together and the code used can be found in the Rmd files stored in the root directory of this repository, while the bookdown-rendered webbook is available at:

[alberdilab.github.io/Podarcis_adaptation](https://alberdilab.github.io/Podarcis_adaptation)

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the Rmd documents. Note that the code chunks that require heavy computation have been tuned off using 'eval=FALSE'. To re-render the webbook, you can use the following code:

```
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```

