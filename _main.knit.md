---
title: "AlberdiLab | Martin-Bideguren et al. 2024"
subtitle: "Study title to be added"
author:
  - Garazi Martin-Bideguren^[University of Copenhagen, garazi.bideguren@sund.ku.dk], Ostaizka Aizpurua^[University of Copenhagen, ostaizka.aizpurua@sund.ku.dk], Carlos Cabido^[Sociedad de Ciencias Aranzadi-Departamento de Herpetología, ccabido@aranzadi.eus] and Antton Alberdi^[University of Copenhagen, antton.alberdi@sund.ku.dk]
date: "Last update: 2024-08-19"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
url: https://alberdilab.github.io/Podarcis_adaptation
description: |
  Data analysis code for the study on the recovery of metagenome‑assembled genomes and derived microbial communities from lizard faecal samples.
link-citations: yes
github-repo: alberdilab/Podarcis_adaptation
---



# Introduction

This webbook contains all the code used for data analysis in study of the individual-level metagenomic data of Podarcis muralis and Podarcis liolepis lizards from different environments during an experimental setup.

## Prepare the R environment

### Environment

To reproduce all the analyses locally, clone this repository in your computer using:

```
RStudio > New Project > Version Control > Git
```

And indicating the following git repository:

> https://github.com/alberdilab/Podarcis_adaptation.git

Once the R project has been created, follow the instructions and code chunks shown in this webbook.

### Libraries

The following R packages are required for the data analysis.


```{.r .script-source}
# Base
library(R.utils)
library(knitr)
library(tidyverse)
library(devtools)
library(tinytable)
library(janitor)

# For tree handling
library(ape)
library(phyloseq)
library(phytools)

# For plotting
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggnewscale)
library(gridExtra)
library(ggtreeExtra)
library(ggtree)
library(ggh4x)

# For statistics
library(spaa)
library(vegan)
library(Rtsne)
library(geiger)
library(hilldiv2)
library(distillR)
library(broom.mixed)
library(corrplot)
library(nlme)
library(pairwiseAdonis)
```

<!--chapter:end:index.Rmd-->

# Prepare data

## Load data

Load the original data files outputted by the bioinformatic pipeline.

### Sample metadata


```{.r .script-source}
sample_metadata <- read_tsv("data/metadata.tsv")
```

### Read counts


```{.r .script-source}
read_counts_raw <- read_tsv("data/genome.count.tsv") %>%
    rename(genome=1)

#Transformation of read_counts to combine data from both sequence rounds
merge_and_rename <- function(read_counts_raw) {
  read_counts_raw %>%
    # Gather the columns into long format
    pivot_longer(cols = -genome, names_to = "col") %>%
    # Extract prefix
    mutate(prefix = gsub("^(.*?)_.*", "\\1", col)) %>%
    # Group by prefix and genome, then summarize
    group_by(prefix, genome) %>%
    summarize(value = sum(value)) %>%
    # Spread the data back to wide format
    pivot_wider(names_from = prefix, values_from = value)
}

read_counts <- merge_and_rename(read_counts_raw)
```

### Genome base hits


```{.r .script-source}
genome_hits_raw <- read_tsv("data/genome.covered_bases.tsv") %>%
    rename(genome=1)

#Transformation of genome_hits to combine data from both sequence rounds
merge_and_rename <- function(genome_hits_raw) {
  genome_hits_raw %>%
    # Gather the columns into long format
    pivot_longer(cols = -genome, names_to = "col") %>%
    # Extract prefix
    mutate(prefix = gsub("^(.*?)_.*", "\\1", col)) %>%
    # Group by prefix and genome, then summarize
    group_by(prefix, genome) %>%
    summarize(value = sum(value)) %>%
    # Spread the data back to wide format
    pivot_wider(names_from = prefix, values_from = value)
}

genome_hits <- merge_and_rename(genome_hits_raw)
```

### Genome taxonomy


```{.r .script-source}
genome_taxonomy <- read_tsv("data/gtdbtk.summary.tsv") %>%
  select(mag_id = user_genome, classification) %>%
  separate(
    classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";") %>%
    rename(genome=1)
```

### Genome quality


```{.r .script-source}
genome_quality <- read_tsv("data/quality_report.tsv") %>%
  select(
    genome = Name, 
    completeness = Completeness, 
    contamination = Contamination,
    length = Genome_Size, 
    gc = GC_Content
  )
genome_quality<-genome_quality %>% 
  mutate (genome=str_remove_all(genome,".fa"))

#Filter MAGs with over 70% completeness and less than 10% contamination
genome_quality <- genome_quality %>%
  filter(completeness > 70 & contamination < 10)
```

### Genome tree


```{.r .script-source}
genome_tree <- read_tree("data/gtdbtk.backbone.bac120.classify.tree")
genome_tree$tip.label <- str_replace_all(genome_tree$tip.label,"'", "") #remove single quotes in MAG names

#Filter genome_taxonomy to keep MAGs with over 70% completeness and less than 10% contamination
genome_taxonomy <- genome_taxonomy %>%
  semi_join(genome_quality, by = "genome")
genome_tree <- keep.tip(genome_tree, tip=genome_taxonomy$genome) # keep only MAG tips
```

### Genome annotations


```{.r .script-source}
genome_annotations <- read_tsv("data/annotations.tsv.xz") %>%
    rename(gene=1, genome=2, contig=3)

#Filter only the MAGs with over 70% completeness and less than 10% contamination
genome_annotations <- genome_annotations %>%
  semi_join(genome_quality, by = "genome")
```

## Create working objects

Transform the original data files into working objects for downstream analyses.

### Merge genome taxonomy and quality


```{.r .script-source}
genome_metadata <- genome_taxonomy %>%
  inner_join(genome_quality,by=join_by(genome==genome)) #join quality
```

### Calculate genome coverage


```{.r .script-source}
#Filter genome_hits for the MAGs with over 70% completeness and less than 10% contamination
genome_hits <- genome_hits %>%
  semi_join(genome_quality, by = "genome")

genome_coverage <- genome_hits %>%
  mutate(across(where(is.numeric), ~ ./genome_metadata$length))
```

## Filtering

4 samples are removed from the analysis due to their low sequencing depth.


```{.r .script-source}
#Counts_raw
columns_to_exclude <- c("AD16","AD23", "AD25", "AD91") # Columns to exclude
read_counts <- read_counts %>%
  select(-columns_to_exclude)

#Metadata
sample_metadata <- sample_metadata %>%
  filter(Tube_code != "AD16") %>%
  filter(Tube_code != "AD23") %>%
  filter(Tube_code != "AD25") %>%
  filter(Tube_code != "AD91")

#Coverage_table
columns_to_exclude <- c("AD16", "AD23", "AD25", "AD91")  # Columns to exclude
genome_coverage <- genome_coverage %>%
  select(-columns_to_exclude)
```


```{.r .script-source}
##Removal of samples of individual LI1_2nd_6
#Counts_raw
columns_to_exclude <- c("AFV13","AD06", "AD31", "AE03", "AF12") # Columns to exclude
read_counts <- read_counts %>%
  select(-columns_to_exclude)

#Metadata
sample_metadata <- sample_metadata %>%
  filter(Tube_code != "AFV13") %>%
  filter(Tube_code != "AD06") %>%
  filter(Tube_code != "AD31") %>%
  filter(Tube_code != "AE03") %>%
  filter(Tube_code != "AF12") 

#Coverage_table
columns_to_exclude <- c("AFV13","AD06", "AD31", "AE03", "AF12")  # Columns to exclude
genome_coverage <- genome_coverage %>%
  select(-columns_to_exclude)
```


### Filter reads by coverage


```{.r .script-source}
#Filter read_counts for the MAGs with over 70% completeness and less than 10% contamination
read_counts <- read_counts %>%
  semi_join(genome_quality, by = "genome")

min_coverage=0.3
read_counts_filt <- genome_coverage %>%
  mutate(across(where(is.numeric), ~ ifelse(. > min_coverage, 1, 0))) %>%
  mutate(across(-1, ~ . * read_counts[[cur_column()]])) 
```

### Transform reads into genome counts


```{.r .script-source}
readlength=150
genome_counts <- read_counts %>%
  mutate(across(where(is.numeric), ~ . / (genome_metadata$length / readlength) ))
```


```{.r .script-source}
readlength=150
genome_counts_filt <- read_counts_filt %>%
  mutate(across(where(is.numeric), ~ . / (genome_metadata$length / readlength) ))
```

### Distill annotations into GIFTs 


```{.r .script-source}
genome_gifts <- distill(genome_annotations,GIFT_db,genomecol=2,annotcol=c(9,10,19))
```

## Load data statistics

### Raw reads


```{.r .script-source}
raw_reads <-
  "data/multiqc_general_stats_2.txt" %>%
  read_tsv() %>%
  select(
    sample = Sample,
    raw_reads = `total_sequences`
  ) %>%
  mutate(sample_prefix = str_extract(sample, "^[^_]+")) %>%
  group_by(sample_prefix) %>%
  summarize(raw_reads = sum(raw_reads, na.rm = TRUE)) %>%
  rename(sample = sample_prefix)  %>%
  mutate(sample = str_remove(sample, "^fastp \\|\\s*"))
```

### Quality-filtered reads


```{.r .script-source}
fastp_reads <-
  "data/multiqc_general_stats.txt" %>%
  read_tsv() %>%
  select(
    sample = Sample,
    trimmed_reads = `total_sequences`
  ) %>%
  mutate(sample_prefix = str_extract(sample, "^[^_]+")) %>%
  group_by(sample_prefix) %>%
  summarize(trimmed_reads = sum(trimmed_reads, na.rm = TRUE)) %>%
  rename(sample = sample_prefix) %>% 
  filter(!str_detect(sample, "nonlizard \\|")) %>%
  filter(!str_detect(sample, "lizard \\|")) %>%
  filter(!str_detect(sample, "refseq500 \\|"))  %>%
  mutate(sample = str_remove(sample, "^fastp \\|\\s*"))
```

### Host-mapped reads


```{.r .script-source}
host_mapped <-
  "data/multiqc_general_stats.txt" %>%
  read_tsv() %>%
  filter(str_detect(Sample, "lizard")) %>%
  select(
    sample = Sample,
    host_mapped = `reads_mapped`,
    mapping_total = `raw_total_sequences`
  ) %>%
  mutate(
    host_unmapped = mapping_total - host_mapped
  ) %>%
  filter(!is.na(host_mapped)) %>%
  separate(
    col = sample,
    into = c("host_name", "sample"),
    sep = " \\| "
  ) %>%
  rename(mapped = host_mapped, unmapped = host_unmapped) %>%
  select(-mapping_total) %>%
  pivot_longer(-host_name:-sample) %>%
  mutate(
    name = str_glue("{name}_{host_name}")
  ) %>%
  select(-host_name) %>%
  pivot_wider() %>%
  mutate(sample = str_extract(sample, "^[^_]+")) %>%
  group_by(sample) %>%
  summarize(
    mapped_lizard = sum(mapped_lizard),
    unmapped_lizard = sum(unmapped_lizard)
  ) 
```

### Prokaryotic fraction


```{.r .script-source}
singlem <-
  "data/singleM.tsv" %>%
  read_tsv() %>%
  distinct() %>%
  mutate(
    sample = sample %>% str_remove_all("_1$"),
   read_fraction = read_fraction %>% str_remove("%") %>% as.numeric(),
   read_fraction = read_fraction / 100
  ) %>%
  select(
   sample,
   singlem_prokaryotic_bases = bacterial_archaeal_bases,
   singlem_metagenome_size = metagenome_size,
   singlem_read_fraction = read_fraction,
  ) %>%
mutate(sample_prefix = str_extract(sample, "^[^_]+")) %>%
  group_by(sample_prefix) %>%
  summarize(
    singlem_prokaryotic_bases = sum(singlem_prokaryotic_bases),
    singlem_metagenome_size = sum(singlem_metagenome_size),
    singlem_read_fraction = mean(singlem_read_fraction)
  ) %>%
  rename(sample = sample_prefix)
```

### MAG-mapped reads


```{.r .script-source}
mag_mapping <-
  read_tsv("data/contig.count.tsv.xz") %>%
  pivot_longer(-sequence_id) %>%
  summarise(value = sum(value), .by = "name") %>%
  rename(sample = name, mapped_mags = value) %>%
  mutate(sample_prefix = str_extract(sample, "^[^_]+")) %>%
    group_by(sample_prefix) %>%
    summarize(
      mapped_mags = sum(mapped_mags)
    ) %>%
    rename(sample = sample_prefix)
```

### Wrap data statistics


```{.r .script-source}
data_stats <- raw_reads %>%
  left_join(fastp_reads) %>%
  left_join(host_mapped) %>%
  left_join(singlem) %>%
  left_join(mag_mapping)

data_stats<- data_stats %>%
  filter(!str_detect(sample, "nonlizard \\|")) %>%
  filter(!str_detect(sample, "lizard \\|")) %>%
  filter(!str_detect(sample, "refseq500 \\|"))
```

## Prepare color scheme

[AlberdiLab](www.alberdilab.dk) projects use unified color schemes developed for the [Earth Hologenome Initiative](www.earthhologenome.org), to facilitate figure interpretation.


```{.r .script-source}
phylum_colors <- read_tsv("https://raw.githubusercontent.com/earthhologenome/EHI_taxonomy_colour/main/ehi_phylum_colors.tsv") %>%
    right_join(genome_metadata, by=join_by(phylum == phylum)) %>%
    arrange(match(genome, genome_tree$tip.label)) %>%
    select(phylum, colors) %>% 
    unique() %>%
    arrange(phylum) %>%
    pull(colors, name=phylum)
```

## Wrap working objects

All working objects are wrapped into a single Rdata object to facilitate downstream usage.


```{.r .script-source}
save(sample_metadata, 
     genome_metadata, 
     read_counts, 
     genome_counts, 
     genome_counts_filt, 
     genome_tree,
     genome_gifts, 
     phylum_colors,
     data_stats,
     file = "data/data.Rdata")
```

<!--chapter:end:01_data_preparation.Rmd-->

# Data statistics


```{.r .script-source}
load("data/data.Rdata")
```

## Sequencing reads statistics


```{.r .script-source}
data_stats$raw_reads %>% sum()
```

```{.script-output}
[1] 4995505290
```

```{.r .script-source}
data_stats$raw_reads %>% mean()
```

```{.script-output}
[1] 28875753
```

```{.r .script-source}
data_stats$raw_reads %>% sd()
```

```{.script-output}
[1] 13447526
```

## DNA fractions

```{.r .script-source}
#Overall
data_stats %>%
    mutate(mapped_perc=mapped_mags/trimmed_reads) %>%
    summarise(mean=mean(mapped_perc),sd=sd(mapped_perc))
```

```{.script-output}
# A tibble: 1 × 2
   mean    sd
  <dbl> <dbl>
1 0.438 0.207
```



```{.r .script-source}
data_stats %>%
  mutate(
    low_quality = raw_reads - trimmed_reads,
    unmapped_reads = trimmed_reads - mapped_lizard - mapped_mags
  ) %>%
  select(sample, low_quality, mapped_lizard, mapped_mags, unmapped_reads) %>%
  pivot_longer(-sample) %>%
  mutate(name=factor(name,levels=c("low_quality","mapped_lizard","unmapped_reads","mapped_mags"))) %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(!is.na(time_point)) %>%
  ggplot(aes(x = sample, y = value, fill = name)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(name="Sequence type",
                    breaks=c("low_quality","mapped_lizard","unmapped_reads","mapped_mags"),
                    labels=c("Low quality","Mapped to host","Unmapped","Mapped to MAGs"),
                    values=c("#CCCCCC", "#bcdee1", "#d8b8a3","#93655c"))+
      facet_grid(~time_point, scales = "free", labeller = labeller(time_point=c("0_Wild"="Wild", "1_Acclimation"="Acclimation", "2_Antibiotics"="Antibiotics",
                                                                                "3_Transplant1"="Transplant1", "4_Transplant2"="Transplant2", "5_Post-FMT1"="Post1",
                                                                                "6_Post-FMT2"="Post2"))) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=0)) +
      labs(y="DNA sequence fraction",x="Samples")
```

<img src="_main_files/figure-html/data_fractions_plot-1.png" width="960" />

## Recovered microbial fraction


```{.r .script-source}
data_stats %>%
  mutate(
    unmapped_reads = trimmed_reads - mapped_lizard - mapped_mags,
    mag_proportion = mapped_mags / (mapped_mags + unmapped_reads),
    singlem_read_fraction = singlem_read_fraction
  ) %>%
  select(sample, mag_proportion, singlem_read_fraction) %>%
  mutate(
    mag_proportion = if_else(singlem_read_fraction == 0, 0, mag_proportion),
    singlem_read_fraction = if_else(singlem_read_fraction == 0, NA, singlem_read_fraction),
    singlem_read_fraction = if_else(singlem_read_fraction < mag_proportion, NA, singlem_read_fraction),
    singlem_read_fraction = if_else(singlem_read_fraction > 1, 1, singlem_read_fraction)
  ) %>%
  pivot_longer(-sample, names_to = "proportion", values_to = "value") %>%
  mutate(
    proportion = factor(
      proportion,
      levels = c("mag_proportion", "singlem_read_fraction")
    )
  ) %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(!is.na(time_point)) %>%
  ggplot(aes(x = sample, y = value, color = proportion)) +
      geom_line(aes(group = sample), color = "#f8a538") +
      geom_point() +
      scale_color_manual(name="Proportion",
                    breaks=c("mag_proportion","singlem_read_fraction"),
                    labels=c("Recovered","Estimated"),
                    values=c("#52e1e8", "#876b53"))+
      theme_minimal() +
      facet_grid(~time_point, scales = "free", space="free", labeller = labeller(time_point=c("0_Wild"="Wild", "1_Acclimation"="Acclimation", "2_Antibiotics"="Antibiotics",
                                                                                "3_Transplant1"="Transplant1", "4_Transplant2"="Transplant2", "5_Post-FMT1"="Post1",
                                                                                "6_Post-FMT2"="Post2")))+
      labs(y = "Samples", x = "Prokaryotic fraction") +
      scale_y_continuous(limits = c(0, 1)) +
      theme(
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1, size = 0),
        legend.position = "right"
      )
```

<img src="_main_files/figure-html/data_estimations_plot-1.png" width="960" />

### Domain-adjusted mapping rate (DAMR)


```{.r .script-source}
data_stats %>%
  mutate(
    unmapped_reads = trimmed_reads - mapped_lizard - mapped_mags,
    mag_proportion = mapped_mags / (mapped_mags + unmapped_reads),
    singlem_read_fraction = singlem_read_fraction
  ) %>%
  mutate(damr=pmin(1, mag_proportion/singlem_read_fraction)) %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(!is.na(time_point)) %>%
  select(sample,damr, time_point, species, type) %>%
  #group_by(type) %>%
  #summarise(mean=mean(damr),sd=sd(damr)) %>%
  tt()
```

```{=html}
<!DOCTYPE html> 
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tinytable_nlvf7q7dt4hrysram3pb</title>
    <style>
.table td.tinytable_css_gbm849g017i44dxx1b9u, .table th.tinytable_css_gbm849g017i44dxx1b9u {    border-bottom: solid 0.1em #d3d8dc; }
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']]
      },
      svg: {
        fontCache: 'global'
      }
    };
    </script>
  </head>

  <body>
    <div class="container">
      <table class="table table-borderless" id="tinytable_nlvf7q7dt4hrysram3pb" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <thead>
        
              <tr>
                <th scope="col">sample</th>
                <th scope="col">damr</th>
                <th scope="col">time_point</th>
                <th scope="col">species</th>
                <th scope="col">type</th>
              </tr>
        </thead>
        
        <tbody>
                <tr>
                  <td>AC79</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC80</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC81</td>
                  <td>0.9813013</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC82</td>
                  <td>0.7490771</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC83</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC84</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC85</td>
                  <td>0.8363777</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC86</td>
                  <td>0.8800275</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC87</td>
                  <td>0.9672854</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC88</td>
                  <td>0.3142108</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC89</td>
                  <td>0.9338962</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC90</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC91</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC92</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AC93</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC94</td>
                  <td>0.9159348</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC95</td>
                  <td>0.9239904</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC96</td>
                  <td>0.9425962</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AC97</td>
                  <td>0.9394886</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AC98</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AC99</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD01</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD02</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD03</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD04</td>
                  <td>0.8896057</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD05</td>
                  <td>1.0000000</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD07</td>
                  <td>0.9515379</td>
                  <td>1_Acclimation</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD08</td>
                  <td>0.7045003</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD09</td>
                  <td>1.0000000</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD10</td>
                  <td>0.9598720</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD11</td>
                  <td>0.6230176</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD12</td>
                  <td>0.7487463</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD13</td>
                  <td>1.0000000</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD14</td>
                  <td>0.4536529</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD15</td>
                  <td>1.0000000</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD17</td>
                  <td>0.8774871</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD18</td>
                  <td>0.4991585</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD19</td>
                  <td>0.4060523</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD20</td>
                  <td>0.5227482</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD21</td>
                  <td>0.6551882</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD22</td>
                  <td>0.9106725</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD24</td>
                  <td>0.5972915</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD26</td>
                  <td>0.9204825</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD27</td>
                  <td>0.5215059</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD28</td>
                  <td>0.8180607</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD29</td>
                  <td>0.7175978</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD30</td>
                  <td>0.9937574</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD32</td>
                  <td>0.6710817</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD33</td>
                  <td>0.6778368</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD34</td>
                  <td>0.4642326</td>
                  <td>2_Antibiotics</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD36</td>
                  <td>0.8288365</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD37</td>
                  <td>0.9530848</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD38</td>
                  <td>0.9247925</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD39</td>
                  <td>1.0000000</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD40</td>
                  <td>0.9991552</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD41</td>
                  <td>0.9550964</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD42</td>
                  <td>1.0000000</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD43</td>
                  <td>1.0000000</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD44</td>
                  <td>0.9108807</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD46</td>
                  <td>0.8708468</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD47</td>
                  <td>0.8328653</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD49</td>
                  <td>0.8940294</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD50</td>
                  <td>0.8757711</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD51</td>
                  <td>0.8479193</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD52</td>
                  <td>0.8596870</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD53</td>
                  <td>0.8528974</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD54</td>
                  <td>1.0000000</td>
                  <td>3_Transplant1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD55</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD56</td>
                  <td>0.9327538</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD57</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD58</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD59</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD60</td>
                  <td>0.7988844</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD61</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD62</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD63</td>
                  <td>0.8992048</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD65</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD66</td>
                  <td>0.8821027</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD68</td>
                  <td>0.8820274</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD70</td>
                  <td>0.9590793</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD71</td>
                  <td>0.8360540</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD72</td>
                  <td>0.8835253</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD73</td>
                  <td>1.0000000</td>
                  <td>4_Transplant2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD74</td>
                  <td>0.9017271</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD75</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD76</td>
                  <td>0.9862581</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD77</td>
                  <td>0.8963376</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD78</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD79</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD80</td>
                  <td>0.8808155</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD81</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD82</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD83</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD84</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD85</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD86</td>
                  <td>0.9889084</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD87</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD88</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AD89</td>
                  <td>0.9562489</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD90</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AD93</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD94</td>
                  <td>0.7947750</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD95</td>
                  <td>0.9506977</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD96</td>
                  <td>0.9171413</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD97</td>
                  <td>0.9257685</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD98</td>
                  <td>0.9048349</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AD99</td>
                  <td>0.9966152</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AE01</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AE02</td>
                  <td>1.0000000</td>
                  <td>5_Post-FMT1</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AE04</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AE05</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AE06</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE07</td>
                  <td>0.7797794</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE08</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE09</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE91</td>
                  <td>0.9259146</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AE92</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE93</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE94</td>
                  <td>0.9152960</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AE95</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AE96</td>
                  <td>0.7786112</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AE97</td>
                  <td>0.9365037</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE98</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AE99</td>
                  <td>0.8513354</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AF01</td>
                  <td>0.9046369</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AF02</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AF03</td>
                  <td>0.9030087</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AF04</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF05</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF06</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF07</td>
                  <td>0.9526443</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF08</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF09</td>
                  <td>1.0000000</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF10</td>
                  <td>0.9716162</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF11</td>
                  <td>0.9654152</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AF13</td>
                  <td>0.8818742</td>
                  <td>6_Post-FMT2</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFU87</td>
                  <td>0.8810513</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFU88</td>
                  <td>0.8317800</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFU91</td>
                  <td>0.8923699</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFU92</td>
                  <td>0.8094776</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFU93</td>
                  <td>0.8517368</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFU94</td>
                  <td>0.8325385</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFU95</td>
                  <td>0.8419270</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFU96</td>
                  <td>0.8326820</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFU97</td>
                  <td>0.8107271</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFU98</td>
                  <td>0.7506522</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFU99</td>
                  <td>0.8582371</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFV01</td>
                  <td>0.9331539</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFV02</td>
                  <td>0.8316460</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Treatment</td>
                </tr>
                <tr>
                  <td>AFV03</td>
                  <td>0.8752591</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFV04</td>
                  <td>0.9180285</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFV05</td>
                  <td>1.0000000</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFV06</td>
                  <td>1.0000000</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFV07</td>
                  <td>0.8460805</td>
                  <td>0_Wild</td>
                  <td>Podarcis_muralis</td>
                  <td>Control</td>
                </tr>
                <tr>
                  <td>AFV08</td>
                  <td>0.7497043</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV09</td>
                  <td>0.5412999</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV10</td>
                  <td>0.8002499</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV11</td>
                  <td>0.8225298</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV12</td>
                  <td>0.7925988</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV14</td>
                  <td>0.8106269</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV15</td>
                  <td>0.9691106</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV16</td>
                  <td>0.8218990</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
                <tr>
                  <td>AFV17</td>
                  <td>0.8091152</td>
                  <td>0_Wild</td>
                  <td>Podarcis_liolepis</td>
                  <td>Hot_control</td>
                </tr>
        </tbody>
      </table>
    </div>

    <script>
      function styleCell_tinytable_ei74qqt4fqqo1xvoydd2(i, j, css_id) {
        var table = document.getElementById("tinytable_nlvf7q7dt4hrysram3pb");
        table.rows[i].cells[j].classList.add(css_id);
      }
      function insertSpanRow(i, colspan, content) {
        var table = document.getElementById('tinytable_nlvf7q7dt4hrysram3pb');
        var newRow = table.insertRow(i);
        var newCell = newRow.insertCell(0);
        newCell.setAttribute("colspan", colspan);
        // newCell.innerText = content;
        // this may be unsafe, but innerText does not interpret <br>
        newCell.innerHTML = content;
      }
      function spanCell_tinytable_ei74qqt4fqqo1xvoydd2(i, j, rowspan, colspan) {
        var table = document.getElementById("tinytable_nlvf7q7dt4hrysram3pb");
        const targetRow = table.rows[i];
        const targetCell = targetRow.cells[j];
        for (let r = 0; r < rowspan; r++) {
          // Only start deleting cells to the right for the first row (r == 0)
          if (r === 0) {
            // Delete cells to the right of the target cell in the first row
            for (let c = colspan - 1; c > 0; c--) {
              if (table.rows[i + r].cells[j + c]) {
                table.rows[i + r].deleteCell(j + c);
              }
            }
          }
          // For rows below the first, delete starting from the target column
          if (r > 0) {
            for (let c = colspan - 1; c >= 0; c--) {
              if (table.rows[i + r] && table.rows[i + r].cells[j]) {
                table.rows[i + r].deleteCell(j);
              }
            }
          }
        }
        // Set rowspan and colspan of the target cell
        targetCell.rowSpan = rowspan;
        targetCell.colSpan = colspan;
      }

window.addEventListener('load', function () { styleCell_tinytable_ei74qqt4fqqo1xvoydd2(0, 0, 'tinytable_css_gbm849g017i44dxx1b9u') })
window.addEventListener('load', function () { styleCell_tinytable_ei74qqt4fqqo1xvoydd2(0, 1, 'tinytable_css_gbm849g017i44dxx1b9u') })
window.addEventListener('load', function () { styleCell_tinytable_ei74qqt4fqqo1xvoydd2(0, 2, 'tinytable_css_gbm849g017i44dxx1b9u') })
window.addEventListener('load', function () { styleCell_tinytable_ei74qqt4fqqo1xvoydd2(0, 3, 'tinytable_css_gbm849g017i44dxx1b9u') })
window.addEventListener('load', function () { styleCell_tinytable_ei74qqt4fqqo1xvoydd2(0, 4, 'tinytable_css_gbm849g017i44dxx1b9u') })
    </script>

  </body>

</html>
```

<!--chapter:end:02_data_statistics.Rmd-->

# MAG catalogue


```{.r .script-source}
load("data/data.Rdata")
```

## Genome phylogeny


```{.r .script-source}
# Generate the phylum color heatmap
phylum_heatmap <- read_tsv("https://raw.githubusercontent.com/earthhologenome/EHI_taxonomy_colour/main/ehi_phylum_colors.tsv") %>%
    right_join(genome_metadata, by=join_by(phylum == phylum)) %>%
    arrange(match(genome, genome_tree$tip.label)) %>%
    select(genome,phylum) %>%
    mutate(phylum = factor(phylum, levels = unique(phylum))) %>%
    column_to_rownames(var = "genome")

# Generate  basal tree
circular_tree <- force.ultrametric(genome_tree, method="extend") %>% # extend to ultrametric for the sake of visualisation
    ggtree(., layout="fan", open.angle=10, size=0.5)
```

```{.script-output}
***************************************************************
*                          Note:                              *
*    force.ultrametric does not include a formal method to    *
*    ultrametricize a tree & should only be used to coerce    *
*   a phylogeny that fails is.ultrametric due to rounding --  *
*    not as a substitute for formal rate-smoothing methods.   *
***************************************************************
```

```{.r .script-source}
# Add phylum ring
circular_tree <- gheatmap(circular_tree, phylum_heatmap, offset=0.55, width=0.1, colnames=FALSE) +
        scale_fill_manual(values=phylum_colors) +
        geom_tiplab2(size=1, hjust=-0.1) +
        theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), panel.margin = margin(0, 0, 0, 0))

# Flush color scale to enable a new color scheme in the next ring
circular_tree <- circular_tree + new_scale_fill()

# Add completeness ring
circular_tree <- circular_tree +
        new_scale_fill() +
        scale_fill_gradient(low = "#d1f4ba", high = "#f4baba") +
        geom_fruit(
                data=genome_metadata,
                geom=geom_bar,
                mapping = aes(x=completeness, y=genome, fill=contamination),
                offset = 0.55,
                orientation="y",
              stat="identity")

# Add genome-size ring
circular_tree <-  circular_tree +
        new_scale_fill() +
        scale_fill_manual(values = "#cccccc") +
        geom_fruit(
             data=genome_metadata,
             geom=geom_bar,
             mapping = aes(x=length, y=genome),
                 offset = 0.05,
                 orientation="y",
         stat="identity")

# Add text
circular_tree <-  circular_tree +
        annotate('text', x=3.4, y=0, label='              Phylum', family='ArialMT', size=3.5) +
        annotate('text', x=4.2, y=0, label='                         Genome quality', family='ArialMT', size=3.5) +
        annotate('text', x=4.7, y=0, label='                      Genome size', family='ArialMT', size=3.5)

#Plot circular tree
circular_tree %>% open_tree(30) %>% rotate_tree(90)
```

<img src="_main_files/figure-html/genome_phylogeny-1.png" width="960" />

## Genome quality


```{.r .script-source}
genome_metadata$completeness %>% mean()
```

```{.script-output}
[1] 92.98658
```

```{.r .script-source}
genome_metadata$completeness %>% sd()
```

```{.script-output}
[1] 7.170467
```

```{.r .script-source}
genome_metadata$contamination %>% mean()
```

```{.script-output}
[1] 2.154888
```

```{.r .script-source}
genome_metadata$contamination %>% sd()
```

```{.script-output}
[1] 2.194372
```


```{.r .script-source}
#Generate quality biplot
genome_biplot <- genome_metadata %>%
  select(c(genome,domain,phylum,completeness,contamination,length)) %>%
  arrange(match(genome, rev(genome_tree$tip.label))) %>% #sort MAGs according to phylogenetic tree
  ggplot(aes(x=completeness,y=contamination,size=length,color=phylum)) +
              geom_point(alpha=0.7) +
                    ylim(c(10,0)) +
                    scale_color_manual(values=phylum_colors) +
                    labs(y= "Contamination", x = "Completeness") +
                    theme_classic() +
                    theme(legend.position = "none")

#Generate contamination boxplot
genome_contamination <- genome_metadata %>%
            ggplot(aes(y=contamination)) +
                    ylim(c(10,0)) +
                    geom_boxplot(colour = "#999999", fill="#cccccc") +
                    theme_void() +
                    theme(legend.position = "none",
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        plot.margin = unit(c(0, 0, 0.40, 0),"inches")) #add bottom-margin (top, right, bottom, left)

#Generate completeness boxplot
genome_completeness <- genome_metadata %>%
        ggplot(aes(x=completeness)) +
                xlim(c(50,100)) +
                geom_boxplot(colour = "#999999", fill="#cccccc") +
                theme_void() +
                theme(legend.position = "none",
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    plot.margin = unit(c(0, 0, 0, 0.50),"inches")) #add left-margin (top, right, bottom, left)

#Render composite figure
grid.arrange(grobs = list(genome_completeness,genome_biplot,genome_contamination),
        layout_matrix = rbind(c(1,1,1,1,1,1,1,1,1,1,1,4),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3),
                              c(2,2,2,2,2,2,2,2,2,2,2,3)))
```

<img src="_main_files/figure-html/genome_quality_plot-1.png" width="960" />


## Functional overview


```{.r .script-source}
# Aggregate basal GIFT into elements
function_table <- genome_gifts %>%
    to.elements(., GIFT_db)

# Generate  basal tree
function_tree <- force.ultrametric(genome_tree, method="extend") %>%
                ggtree(., size = 0.3) 
```

```{.script-output}
***************************************************************
*                          Note:                              *
*    force.ultrametric does not include a formal method to    *
*    ultrametricize a tree & should only be used to coerce    *
*   a phylogeny that fails is.ultrametric due to rounding --  *
*    not as a substitute for formal rate-smoothing methods.   *
***************************************************************
```

```{.r .script-source}
#Add phylum colors next to the tree tips
function_tree <- gheatmap(function_tree, phylum_heatmap, offset=0, width=0.1, colnames=FALSE) +
            scale_fill_manual(values=phylum_colors) +
            labs(fill="Phylum")

#Reset fill scale to use a different colour profile in the heatmap
function_tree <- function_tree + new_scale_fill()

#Add functions heatmap
function_tree <- gheatmap(function_tree, function_table, offset=0.5, width=3.5, colnames=FALSE) +
            vexpand(.08) +
            coord_cartesian(clip = "off") +
            scale_fill_gradient(low = "#f4f4f4", high = "steelblue", na.value="white") +
            labs(fill="GIFT")

#Reset fill scale to use a different colour profile in the heatmap
function_tree <- function_tree + new_scale_fill()

# Add completeness barplots
function_tree <- function_tree +
            geom_fruit(data=genome_metadata,
            geom=geom_bar,
            grid.params=list(axis="x", text.size=2, nbreak = 1),
            axis.params=list(vline=TRUE),
            mapping = aes(x=length, y=genome, fill=completeness),
                 offset = 3.8,
                 orientation="y",
                 stat="identity") +
            scale_fill_gradient(low = "#cf8888", high = "#a2cc87") +
            labs(fill="Genome\ncompleteness")

function_tree
```

<img src="_main_files/figure-html/function_heatmap-1.png" width="960" />

## Functional ordination


```{.r .script-source}
# Generate the tSNE ordination
tSNE_function <- Rtsne(X=function_table, dims = 2, check_duplicates = FALSE)

# Plot the ordination
function_ordination <- tSNE_function$Y %>%
                as.data.frame() %>%
                mutate(genome=rownames(function_table)) %>%
                inner_join(genome_metadata, by="genome") %>%
                rename(tSNE1="V1", tSNE2="V2") %>%
                select(genome,phylum,tSNE1,tSNE2, length) %>%
                ggplot(aes(x = tSNE1, y = tSNE2, color = phylum, size=length))+
                            geom_point(shape=16, alpha=0.7) +
                            scale_color_manual(values=phylum_colors) +
                            theme_minimal() +
                labs(color="Phylum", size="Genome size") +
                guides(color = guide_legend(override.aes = list(size = 5))) # enlarge Phylum dots in legend

function_ordination
```

<img src="_main_files/figure-html/function_ordination-1.png" width="960" />

<!--chapter:end:03_mag_catalogue.Rmd-->

# Community composition


```{.r .script-source}
load("data/data.Rdata")
```

## Taxonomy overview 

### Stacked barplot


```{.r .script-source}
# Merge data frames based on sample
transplants_metadata<-sample_metadata%>%
  mutate(Tube_code=str_remove_all(Tube_code, "_a"))
transplants_metadata$newID <- paste(transplants_metadata$Tube_code, "_", transplants_metadata$individual)

merged_data<-genome_counts_filt %>%
  mutate_at(vars(-genome),~./sum(.)) %>% #apply TSS normalisation
  pivot_longer(-genome, names_to = "sample", values_to = "count") %>% #reduce to minimum number of columns
  left_join(., genome_metadata, by = join_by(genome == genome)) %>% #append genome metadata
  left_join(., transplants_metadata, by = join_by(sample == Tube_code)) %>% #append sample metadata
  filter(count > 0) #filter 0 counts

ggplot(merged_data, aes(x=newID,y=count, fill=phylum, group=phylum)) + #grouping enables keeping the same sorting of taxonomic units
    geom_bar(stat="identity", colour="white", linewidth=0.1) + #plot stacked bars with white borders
    scale_fill_manual(values=phylum_colors) +
    facet_nested(. ~ time_point + type ,  scales="free") + #facet per day and treatment
    guides(fill = guide_legend(ncol = 1)) +
    labs(fill="Phylum",y = "Relative abundance",x="Sample")+
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=0))
```

<img src="_main_files/figure-html/taxonomy_barplot-1.png" width="1728" />

#### Wild samples


```{.r .script-source}
merged_data  %>%
  filter(time_point=="0_Wild")  %>%
  ggplot(aes(x=newID,y=count, fill=phylum, group=phylum)) + #grouping enables keeping the same sorting of taxonomic units
    geom_bar(stat="identity", colour="white", linewidth=0.1) + #plot stacked bars with white borders
    scale_fill_manual(values=phylum_colors) +
    facet_nested(. ~ Population,  scales="free") + #facet per day and treatment
    guides(fill = guide_legend(ncol = 1)) +
    labs(fill="Phylum",y = "Relative abundance",x="Sample")+
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=0),
    strip.text.x = element_text(size = 12)
    )
```

<img src="_main_files/figure-html/taxonomy_barplot_wild-1.png" width="1728" />


### Phylum relative abundances

```{.r .script-source}
phylum_summary <- genome_counts_filt %>%
  mutate_at(vars(-genome),~./sum(.)) %>% #apply TSS normalisation
  pivot_longer(-genome, names_to = "sample", values_to = "count") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  left_join(genome_metadata, by = join_by(genome == genome)) %>%
  group_by(sample,phylum) %>%
  summarise(relabun=sum(count))

phylum_summary %>%
    group_by(phylum) %>%
    summarise(total_mean=mean(relabun*100, na.rm=T),
              total_sd=sd(relabun*100, na.rm=T))  %>%
    mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
    arrange(-total_mean) %>% 
    dplyr::select(phylum,total) %>% 
    tt()
```

```{=html}
<!DOCTYPE html> 
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tinytable_2o6z8s9lxcwy24wc3jct</title>
    <style>
.table td.tinytable_css_neflugve7g9s51slg47e, .table th.tinytable_css_neflugve7g9s51slg47e {    border-bottom: solid 0.1em #d3d8dc; }
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']]
      },
      svg: {
        fontCache: 'global'
      }
    };
    </script>
  </head>

  <body>
    <div class="container">
      <table class="table table-borderless" id="tinytable_2o6z8s9lxcwy24wc3jct" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <thead>
        
              <tr>
                <th scope="col">phylum</th>
                <th scope="col">total</th>
              </tr>
        </thead>
        
        <tbody>
                <tr>
                  <td>p__Bacteroidota</td>
                  <td>38.75±19.94</td>
                </tr>
                <tr>
                  <td>p__Bacillota_A</td>
                  <td>24.77±15.73</td>
                </tr>
                <tr>
                  <td>p__Bacillota</td>
                  <td>11.96±14.79</td>
                </tr>
                <tr>
                  <td>p__Pseudomonadota</td>
                  <td>9.41±15.86</td>
                </tr>
                <tr>
                  <td>p__Campylobacterota</td>
                  <td>5.49±9.42</td>
                </tr>
                <tr>
                  <td>p__Verrucomicrobiota</td>
                  <td>2.78±6.7</td>
                </tr>
                <tr>
                  <td>p__Desulfobacterota</td>
                  <td>2.37±3.68</td>
                </tr>
                <tr>
                  <td>p__Chlamydiota</td>
                  <td>1.1±6.08</td>
                </tr>
                <tr>
                  <td>p__Fusobacteriota</td>
                  <td>1.06±2.86</td>
                </tr>
                <tr>
                  <td>p__Cyanobacteriota</td>
                  <td>0.93±1.66</td>
                </tr>
                <tr>
                  <td>p__Bacillota_C</td>
                  <td>0.48±0.67</td>
                </tr>
                <tr>
                  <td>p__Spirochaetota</td>
                  <td>0.41±1.25</td>
                </tr>
                <tr>
                  <td>p__Bacillota_B</td>
                  <td>0.26±0.5</td>
                </tr>
                <tr>
                  <td>p__Actinomycetota</td>
                  <td>0.13±0.65</td>
                </tr>
                <tr>
                  <td>p__Elusimicrobiota</td>
                  <td>0.11±0.63</td>
                </tr>
        </tbody>
      </table>
    </div>

    <script>
      function styleCell_tinytable_drz0devqgjq432zgtuxf(i, j, css_id) {
        var table = document.getElementById("tinytable_2o6z8s9lxcwy24wc3jct");
        table.rows[i].cells[j].classList.add(css_id);
      }
      function insertSpanRow(i, colspan, content) {
        var table = document.getElementById('tinytable_2o6z8s9lxcwy24wc3jct');
        var newRow = table.insertRow(i);
        var newCell = newRow.insertCell(0);
        newCell.setAttribute("colspan", colspan);
        // newCell.innerText = content;
        // this may be unsafe, but innerText does not interpret <br>
        newCell.innerHTML = content;
      }
      function spanCell_tinytable_drz0devqgjq432zgtuxf(i, j, rowspan, colspan) {
        var table = document.getElementById("tinytable_2o6z8s9lxcwy24wc3jct");
        const targetRow = table.rows[i];
        const targetCell = targetRow.cells[j];
        for (let r = 0; r < rowspan; r++) {
          // Only start deleting cells to the right for the first row (r == 0)
          if (r === 0) {
            // Delete cells to the right of the target cell in the first row
            for (let c = colspan - 1; c > 0; c--) {
              if (table.rows[i + r].cells[j + c]) {
                table.rows[i + r].deleteCell(j + c);
              }
            }
          }
          // For rows below the first, delete starting from the target column
          if (r > 0) {
            for (let c = colspan - 1; c >= 0; c--) {
              if (table.rows[i + r] && table.rows[i + r].cells[j]) {
                table.rows[i + r].deleteCell(j);
              }
            }
          }
        }
        // Set rowspan and colspan of the target cell
        targetCell.rowSpan = rowspan;
        targetCell.colSpan = colspan;
      }

window.addEventListener('load', function () { styleCell_tinytable_drz0devqgjq432zgtuxf(0, 0, 'tinytable_css_neflugve7g9s51slg47e') })
window.addEventListener('load', function () { styleCell_tinytable_drz0devqgjq432zgtuxf(0, 1, 'tinytable_css_neflugve7g9s51slg47e') })
    </script>

  </body>

</html>
```

#### Wild samples


```{.r .script-source}
#Cold and wet population
phylum_summary %>%
    left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
    group_by(phylum, Population) %>%
    filter(Population=="Cold_wet") %>%
    filter(time_point=="0_Wild") %>%
    summarise(total_mean=mean(relabun*100, na.rm=T),
              total_sd=sd(relabun*100, na.rm=T))  %>%
    mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
    arrange(-total_mean) %>% 
    dplyr::select(phylum,total, Population) %>% 
    tt()
```

```{=html}
<!DOCTYPE html> 
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tinytable_5nostl892bjkhrlx7n9c</title>
    <style>
.table td.tinytable_css_fh4sjam9veg5mwduqeyf, .table th.tinytable_css_fh4sjam9veg5mwduqeyf {    border-bottom: solid 0.1em #d3d8dc; }
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']]
      },
      svg: {
        fontCache: 'global'
      }
    };
    </script>
  </head>

  <body>
    <div class="container">
      <table class="table table-borderless" id="tinytable_5nostl892bjkhrlx7n9c" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <thead>
        
              <tr>
                <th scope="col">phylum</th>
                <th scope="col">total</th>
                <th scope="col">Population</th>
              </tr>
        </thead>
        
        <tbody>
                <tr>
                  <td>p__Bacteroidota</td>
                  <td>41.93±15.4</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Bacillota_A</td>
                  <td>33.96±16.7</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Bacillota</td>
                  <td>5.73±4.73</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Campylobacterota</td>
                  <td>5.39±4.64</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Verrucomicrobiota</td>
                  <td>4.36±3.06</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Pseudomonadota</td>
                  <td>3.75±3.76</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Desulfobacterota</td>
                  <td>2.31±1.6</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Bacillota_C</td>
                  <td>0.85±0.82</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Cyanobacteriota</td>
                  <td>0.64±0.83</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Bacillota_B</td>
                  <td>0.5±0.77</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Fusobacteriota</td>
                  <td>0.37±1.28</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Spirochaetota</td>
                  <td>0.12±0.35</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Actinomycetota</td>
                  <td>0.09±0.22</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Chlamydiota</td>
                  <td>0±0</td>
                  <td>Cold_wet</td>
                </tr>
                <tr>
                  <td>p__Elusimicrobiota</td>
                  <td>0±0</td>
                  <td>Cold_wet</td>
                </tr>
        </tbody>
      </table>
    </div>

    <script>
      function styleCell_tinytable_srcqauw0euffq4ykj92j(i, j, css_id) {
        var table = document.getElementById("tinytable_5nostl892bjkhrlx7n9c");
        table.rows[i].cells[j].classList.add(css_id);
      }
      function insertSpanRow(i, colspan, content) {
        var table = document.getElementById('tinytable_5nostl892bjkhrlx7n9c');
        var newRow = table.insertRow(i);
        var newCell = newRow.insertCell(0);
        newCell.setAttribute("colspan", colspan);
        // newCell.innerText = content;
        // this may be unsafe, but innerText does not interpret <br>
        newCell.innerHTML = content;
      }
      function spanCell_tinytable_srcqauw0euffq4ykj92j(i, j, rowspan, colspan) {
        var table = document.getElementById("tinytable_5nostl892bjkhrlx7n9c");
        const targetRow = table.rows[i];
        const targetCell = targetRow.cells[j];
        for (let r = 0; r < rowspan; r++) {
          // Only start deleting cells to the right for the first row (r == 0)
          if (r === 0) {
            // Delete cells to the right of the target cell in the first row
            for (let c = colspan - 1; c > 0; c--) {
              if (table.rows[i + r].cells[j + c]) {
                table.rows[i + r].deleteCell(j + c);
              }
            }
          }
          // For rows below the first, delete starting from the target column
          if (r > 0) {
            for (let c = colspan - 1; c >= 0; c--) {
              if (table.rows[i + r] && table.rows[i + r].cells[j]) {
                table.rows[i + r].deleteCell(j);
              }
            }
          }
        }
        // Set rowspan and colspan of the target cell
        targetCell.rowSpan = rowspan;
        targetCell.colSpan = colspan;
      }

window.addEventListener('load', function () { styleCell_tinytable_srcqauw0euffq4ykj92j(0, 0, 'tinytable_css_fh4sjam9veg5mwduqeyf') })
window.addEventListener('load', function () { styleCell_tinytable_srcqauw0euffq4ykj92j(0, 1, 'tinytable_css_fh4sjam9veg5mwduqeyf') })
window.addEventListener('load', function () { styleCell_tinytable_srcqauw0euffq4ykj92j(0, 2, 'tinytable_css_fh4sjam9veg5mwduqeyf') })
    </script>

  </body>

</html>
```

```{.r .script-source}
#Hot and dry population
phylum_summary %>%
    left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
    group_by(phylum, Population) %>%
    filter(Population=="Hot_dry")%>%
    filter(time_point=="0_Wild") %>%
    summarise(total_mean=mean(relabun*100, na.rm=T),
              total_sd=sd(relabun*100, na.rm=T))  %>%
    mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
    arrange(-total_mean) %>% 
    dplyr::select(phylum,total, Population) %>% 
    tt()
```

```{=html}
<!DOCTYPE html> 
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tinytable_anxvlabbwh6foy6y8g8b</title>
    <style>
.table td.tinytable_css_h7tlaz5d1dqswu6pieph, .table th.tinytable_css_h7tlaz5d1dqswu6pieph {    border-bottom: solid 0.1em #d3d8dc; }
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']]
      },
      svg: {
        fontCache: 'global'
      }
    };
    </script>
  </head>

  <body>
    <div class="container">
      <table class="table table-borderless" id="tinytable_anxvlabbwh6foy6y8g8b" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <thead>
        
              <tr>
                <th scope="col">phylum</th>
                <th scope="col">total</th>
                <th scope="col">Population</th>
              </tr>
        </thead>
        
        <tbody>
                <tr>
                  <td>p__Bacillota_A</td>
                  <td>33.9±16.17</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Bacteroidota</td>
                  <td>27.63±17.64</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Campylobacterota</td>
                  <td>13.53±20.98</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Pseudomonadota</td>
                  <td>10.92±9.94</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Bacillota</td>
                  <td>6.09±7.66</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Spirochaetota</td>
                  <td>2.66±2.78</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Desulfobacterota</td>
                  <td>1.69±1.85</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Fusobacteriota</td>
                  <td>1.34±1.73</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Bacillota_B</td>
                  <td>0.76±0.73</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Bacillota_C</td>
                  <td>0.73±0.5</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Cyanobacteriota</td>
                  <td>0.52±0.65</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Chlamydiota</td>
                  <td>0.12±0.18</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Verrucomicrobiota</td>
                  <td>0.06±0.1</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Elusimicrobiota</td>
                  <td>0.05±0.14</td>
                  <td>Hot_dry</td>
                </tr>
                <tr>
                  <td>p__Actinomycetota</td>
                  <td>0±0</td>
                  <td>Hot_dry</td>
                </tr>
        </tbody>
      </table>
    </div>

    <script>
      function styleCell_tinytable_0zmyr6g34u16anh12ul6(i, j, css_id) {
        var table = document.getElementById("tinytable_anxvlabbwh6foy6y8g8b");
        table.rows[i].cells[j].classList.add(css_id);
      }
      function insertSpanRow(i, colspan, content) {
        var table = document.getElementById('tinytable_anxvlabbwh6foy6y8g8b');
        var newRow = table.insertRow(i);
        var newCell = newRow.insertCell(0);
        newCell.setAttribute("colspan", colspan);
        // newCell.innerText = content;
        // this may be unsafe, but innerText does not interpret <br>
        newCell.innerHTML = content;
      }
      function spanCell_tinytable_0zmyr6g34u16anh12ul6(i, j, rowspan, colspan) {
        var table = document.getElementById("tinytable_anxvlabbwh6foy6y8g8b");
        const targetRow = table.rows[i];
        const targetCell = targetRow.cells[j];
        for (let r = 0; r < rowspan; r++) {
          // Only start deleting cells to the right for the first row (r == 0)
          if (r === 0) {
            // Delete cells to the right of the target cell in the first row
            for (let c = colspan - 1; c > 0; c--) {
              if (table.rows[i + r].cells[j + c]) {
                table.rows[i + r].deleteCell(j + c);
              }
            }
          }
          // For rows below the first, delete starting from the target column
          if (r > 0) {
            for (let c = colspan - 1; c >= 0; c--) {
              if (table.rows[i + r] && table.rows[i + r].cells[j]) {
                table.rows[i + r].deleteCell(j);
              }
            }
          }
        }
        // Set rowspan and colspan of the target cell
        targetCell.rowSpan = rowspan;
        targetCell.colSpan = colspan;
      }

window.addEventListener('load', function () { styleCell_tinytable_0zmyr6g34u16anh12ul6(0, 0, 'tinytable_css_h7tlaz5d1dqswu6pieph') })
window.addEventListener('load', function () { styleCell_tinytable_0zmyr6g34u16anh12ul6(0, 1, 'tinytable_css_h7tlaz5d1dqswu6pieph') })
window.addEventListener('load', function () { styleCell_tinytable_0zmyr6g34u16anh12ul6(0, 2, 'tinytable_css_h7tlaz5d1dqswu6pieph') })
    </script>

  </body>

</html>
```


```{.r .script-source}
phylum_arrange <- phylum_summary %>%
    group_by(phylum) %>%
    summarise(mean=mean(relabun)) %>%
    arrange(-mean) %>%
    select(phylum) %>%
    pull()

phylum_summary %>%
    filter(phylum %in% phylum_arrange) %>%
    mutate(phylum=factor(phylum,levels=rev(phylum_arrange))) %>%
    ggplot(aes(x=relabun, y=phylum, group=phylum, color=phylum)) +
        scale_color_manual(values=phylum_colors[rev(phylum_arrange)]) +
        geom_jitter(alpha=0.5) + 
        theme_minimal() + 
        theme(legend.position="none") +
        labs(y="Phylum",x="Relative abundance")
```

<img src="_main_files/figure-html/taxonomy_boxplot_phylum-1.png" width="960" />

## Taxonomy boxplot

### Family

```{.r .script-source}
family_summary <- genome_counts_filt %>%
  mutate_at(vars(-genome),~./sum(.)) %>% #apply TSS normalisation
  pivot_longer(-genome, names_to = "sample", values_to = "count") %>% #reduce to minimum number of columns
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>% #append sample metadata
  left_join(., genome_metadata, by = join_by(genome == genome)) %>% #append genome metadata
  group_by(sample,family) %>%
  summarise(relabun=sum(count))

family_summary %>%
    group_by(family) %>%
    summarise(mean=mean(relabun, na.rm=TRUE),sd=sd(relabun, na.rm=TRUE)) %>%
    arrange(-mean) %>%
    tt()
```

```{=html}
<!DOCTYPE html> 
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tinytable_xtfkxaazluwug0z7bqhd</title>
    <style>
.table td.tinytable_css_722l3c3vc4uq6hyg1l9a, .table th.tinytable_css_722l3c3vc4uq6hyg1l9a {    border-bottom: solid 0.1em #d3d8dc; }
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']]
      },
      svg: {
        fontCache: 'global'
      }
    };
    </script>
  </head>

  <body>
    <div class="container">
      <table class="table table-borderless" id="tinytable_xtfkxaazluwug0z7bqhd" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <thead>
        
              <tr>
                <th scope="col">family</th>
                <th scope="col">mean</th>
                <th scope="col">sd</th>
              </tr>
        </thead>
        
        <tbody>
                <tr>
                  <td>f__Bacteroidaceae</td>
                  <td>2.260146e-01</td>
                  <td>0.1384706235</td>
                </tr>
                <tr>
                  <td>f__Lachnospiraceae</td>
                  <td>1.410833e-01</td>
                  <td>0.1062953893</td>
                </tr>
                <tr>
                  <td>f__Tannerellaceae</td>
                  <td>1.045659e-01</td>
                  <td>0.0799894745</td>
                </tr>
                <tr>
                  <td>f__Helicobacteraceae</td>
                  <td>5.448546e-02</td>
                  <td>0.0937279764</td>
                </tr>
                <tr>
                  <td>f__Mycoplasmoidaceae</td>
                  <td>3.756572e-02</td>
                  <td>0.0767893776</td>
                </tr>
                <tr>
                  <td>f__Erysipelotrichaceae</td>
                  <td>3.536287e-02</td>
                  <td>0.0452140267</td>
                </tr>
                <tr>
                  <td>f__UBA3700</td>
                  <td>3.456595e-02</td>
                  <td>0.0565495010</td>
                </tr>
                <tr>
                  <td>f__Marinifilaceae</td>
                  <td>2.794365e-02</td>
                  <td>0.0272474083</td>
                </tr>
                <tr>
                  <td>f__Rikenellaceae</td>
                  <td>2.725202e-02</td>
                  <td>0.0471735191</td>
                </tr>
                <tr>
                  <td>f__Enterobacteriaceae</td>
                  <td>2.687327e-02</td>
                  <td>0.0929254305</td>
                </tr>
                <tr>
                  <td>f__Coprobacillaceae</td>
                  <td>2.627823e-02</td>
                  <td>0.0907456387</td>
                </tr>
                <tr>
                  <td>f__</td>
                  <td>2.465083e-02</td>
                  <td>0.0781013121</td>
                </tr>
                <tr>
                  <td>f__Desulfovibrionaceae</td>
                  <td>2.373771e-02</td>
                  <td>0.0367718116</td>
                </tr>
                <tr>
                  <td>f__DTU072</td>
                  <td>2.183007e-02</td>
                  <td>0.0377159876</td>
                </tr>
                <tr>
                  <td>f__Ruminococcaceae</td>
                  <td>1.832093e-02</td>
                  <td>0.0428115194</td>
                </tr>
                <tr>
                  <td>f__Rhizobiaceae</td>
                  <td>1.579679e-02</td>
                  <td>0.0779688169</td>
                </tr>
                <tr>
                  <td>f__LL51</td>
                  <td>1.556592e-02</td>
                  <td>0.0616955422</td>
                </tr>
                <tr>
                  <td>f__UBA3830</td>
                  <td>1.512118e-02</td>
                  <td>0.0441855701</td>
                </tr>
                <tr>
                  <td>f__Akkermansiaceae</td>
                  <td>1.224165e-02</td>
                  <td>0.0317653885</td>
                </tr>
                <tr>
                  <td>f__Chlamydiaceae</td>
                  <td>1.096188e-02</td>
                  <td>0.0607502761</td>
                </tr>
                <tr>
                  <td>f__Fusobacteriaceae</td>
                  <td>1.055726e-02</td>
                  <td>0.0286383947</td>
                </tr>
                <tr>
                  <td>f__CAG-239</td>
                  <td>9.138683e-03</td>
                  <td>0.0152490206</td>
                </tr>
                <tr>
                  <td>f__Enterococcaceae</td>
                  <td>8.437601e-03</td>
                  <td>0.0473906561</td>
                </tr>
                <tr>
                  <td>f__Gastranaerophilaceae</td>
                  <td>7.848357e-03</td>
                  <td>0.0146292952</td>
                </tr>
                <tr>
                  <td>f__Oscillospiraceae</td>
                  <td>6.624721e-03</td>
                  <td>0.0075288565</td>
                </tr>
                <tr>
                  <td>f__UBA1997</td>
                  <td>6.613196e-03</td>
                  <td>0.0315103296</td>
                </tr>
                <tr>
                  <td>f__Streptococcaceae</td>
                  <td>6.600789e-03</td>
                  <td>0.0348230465</td>
                </tr>
                <tr>
                  <td>f__UBA1242</td>
                  <td>4.266475e-03</td>
                  <td>0.0147768596</td>
                </tr>
                <tr>
                  <td>f__Brevinemataceae</td>
                  <td>4.098862e-03</td>
                  <td>0.0125062564</td>
                </tr>
                <tr>
                  <td>f__Acutalibacteraceae</td>
                  <td>3.498766e-03</td>
                  <td>0.0111374416</td>
                </tr>
                <tr>
                  <td>f__RUG11792</td>
                  <td>2.921450e-03</td>
                  <td>0.0255676374</td>
                </tr>
                <tr>
                  <td>f__Clostridiaceae</td>
                  <td>2.855351e-03</td>
                  <td>0.0174153876</td>
                </tr>
                <tr>
                  <td>f__UBA660</td>
                  <td>2.600647e-03</td>
                  <td>0.0118148140</td>
                </tr>
                <tr>
                  <td>f__Peptococcaceae</td>
                  <td>2.556218e-03</td>
                  <td>0.0049947786</td>
                </tr>
                <tr>
                  <td>f__Acidaminococcaceae</td>
                  <td>1.980431e-03</td>
                  <td>0.0051045211</td>
                </tr>
                <tr>
                  <td>f__CAG-508</td>
                  <td>1.874529e-03</td>
                  <td>0.0065256902</td>
                </tr>
                <tr>
                  <td>f__MGBC116941</td>
                  <td>1.783700e-03</td>
                  <td>0.0077801927</td>
                </tr>
                <tr>
                  <td>f__Moraxellaceae</td>
                  <td>1.540093e-03</td>
                  <td>0.0099192011</td>
                </tr>
                <tr>
                  <td>f__RUG14156</td>
                  <td>1.428152e-03</td>
                  <td>0.0045670616</td>
                </tr>
                <tr>
                  <td>f__Staphylococcaceae</td>
                  <td>1.411833e-03</td>
                  <td>0.0051727635</td>
                </tr>
                <tr>
                  <td>f__Anaerovoracaceae</td>
                  <td>1.410739e-03</td>
                  <td>0.0027876311</td>
                </tr>
                <tr>
                  <td>f__Elusimicrobiaceae</td>
                  <td>1.088209e-03</td>
                  <td>0.0062780187</td>
                </tr>
                <tr>
                  <td>f__CAG-288</td>
                  <td>9.840222e-04</td>
                  <td>0.0061275213</td>
                </tr>
                <tr>
                  <td>f__Anaerotignaceae</td>
                  <td>9.320656e-04</td>
                  <td>0.0041174457</td>
                </tr>
                <tr>
                  <td>f__CALVMC01</td>
                  <td>7.793540e-04</td>
                  <td>0.0044385554</td>
                </tr>
                <tr>
                  <td>f__Eggerthellaceae</td>
                  <td>6.643755e-04</td>
                  <td>0.0021620275</td>
                </tr>
                <tr>
                  <td>f__Massilibacillaceae</td>
                  <td>6.322621e-04</td>
                  <td>0.0016561037</td>
                </tr>
                <tr>
                  <td>f__Mycobacteriaceae</td>
                  <td>6.168531e-04</td>
                  <td>0.0061497354</td>
                </tr>
                <tr>
                  <td>f__UBA1820</td>
                  <td>4.705627e-04</td>
                  <td>0.0013078764</td>
                </tr>
                <tr>
                  <td>f__CAG-274</td>
                  <td>4.686117e-04</td>
                  <td>0.0022415212</td>
                </tr>
                <tr>
                  <td>f__Arcobacteraceae</td>
                  <td>3.928587e-04</td>
                  <td>0.0050156837</td>
                </tr>
                <tr>
                  <td>f__Burkholderiaceae_C</td>
                  <td>3.835606e-04</td>
                  <td>0.0048969735</td>
                </tr>
                <tr>
                  <td>f__Muribaculaceae</td>
                  <td>3.508548e-04</td>
                  <td>0.0009525792</td>
                </tr>
                <tr>
                  <td>f__UBA932</td>
                  <td>3.295199e-04</td>
                  <td>0.0011408058</td>
                </tr>
                <tr>
                  <td>f__Hepatoplasmataceae</td>
                  <td>3.099135e-04</td>
                  <td>0.0039567109</td>
                </tr>
                <tr>
                  <td>f__Rhodobacteraceae</td>
                  <td>3.068016e-04</td>
                  <td>0.0039169801</td>
                </tr>
                <tr>
                  <td>f__Weeksellaceae</td>
                  <td>2.873650e-04</td>
                  <td>0.0032049404</td>
                </tr>
                <tr>
                  <td>f__Eubacteriaceae</td>
                  <td>1.707442e-04</td>
                  <td>0.0006844943</td>
                </tr>
                <tr>
                  <td>f__Sphingobacteriaceae</td>
                  <td>1.561202e-04</td>
                  <td>0.0012685229</td>
                </tr>
                <tr>
                  <td>f__Devosiaceae</td>
                  <td>1.544841e-04</td>
                  <td>0.0015368528</td>
                </tr>
                <tr>
                  <td>f__Pumilibacteraceae</td>
                  <td>1.324439e-04</td>
                  <td>0.0007783049</td>
                </tr>
                <tr>
                  <td>f__WRAU01</td>
                  <td>9.956857e-05</td>
                  <td>0.0012712064</td>
                </tr>
                <tr>
                  <td>f__Peptostreptococcaceae</td>
                  <td>2.371535e-05</td>
                  <td>0.0003027773</td>
                </tr>
        </tbody>
      </table>
    </div>

    <script>
      function styleCell_tinytable_0xta6sf17wbc41ex1nur(i, j, css_id) {
        var table = document.getElementById("tinytable_xtfkxaazluwug0z7bqhd");
        table.rows[i].cells[j].classList.add(css_id);
      }
      function insertSpanRow(i, colspan, content) {
        var table = document.getElementById('tinytable_xtfkxaazluwug0z7bqhd');
        var newRow = table.insertRow(i);
        var newCell = newRow.insertCell(0);
        newCell.setAttribute("colspan", colspan);
        // newCell.innerText = content;
        // this may be unsafe, but innerText does not interpret <br>
        newCell.innerHTML = content;
      }
      function spanCell_tinytable_0xta6sf17wbc41ex1nur(i, j, rowspan, colspan) {
        var table = document.getElementById("tinytable_xtfkxaazluwug0z7bqhd");
        const targetRow = table.rows[i];
        const targetCell = targetRow.cells[j];
        for (let r = 0; r < rowspan; r++) {
          // Only start deleting cells to the right for the first row (r == 0)
          if (r === 0) {
            // Delete cells to the right of the target cell in the first row
            for (let c = colspan - 1; c > 0; c--) {
              if (table.rows[i + r].cells[j + c]) {
                table.rows[i + r].deleteCell(j + c);
              }
            }
          }
          // For rows below the first, delete starting from the target column
          if (r > 0) {
            for (let c = colspan - 1; c >= 0; c--) {
              if (table.rows[i + r] && table.rows[i + r].cells[j]) {
                table.rows[i + r].deleteCell(j);
              }
            }
          }
        }
        // Set rowspan and colspan of the target cell
        targetCell.rowSpan = rowspan;
        targetCell.colSpan = colspan;
      }

window.addEventListener('load', function () { styleCell_tinytable_0xta6sf17wbc41ex1nur(0, 0, 'tinytable_css_722l3c3vc4uq6hyg1l9a') })
window.addEventListener('load', function () { styleCell_tinytable_0xta6sf17wbc41ex1nur(0, 1, 'tinytable_css_722l3c3vc4uq6hyg1l9a') })
window.addEventListener('load', function () { styleCell_tinytable_0xta6sf17wbc41ex1nur(0, 2, 'tinytable_css_722l3c3vc4uq6hyg1l9a') })
    </script>

  </body>

</html>
```



```{.r .script-source}
family_arrange <- family_summary %>%
    group_by(family) %>%
    summarise(mean=sum(relabun)) %>%
    arrange(-mean) %>%
    select(family) %>%
    pull()

family_summary %>%
    left_join(genome_metadata %>% select(family,phylum) %>% unique(),by=join_by(family==family)) %>%
    left_join(sample_metadata,by=join_by(sample==Tube_code)) %>%
    filter(family %in% family_arrange[1:20]) %>%
    mutate(family=factor(family,levels=rev(family_arrange[1:20]))) %>%
    filter(relabun > 0) %>%
    ggplot(aes(x=relabun, y=family, group=family, color=phylum)) +
        scale_color_manual(values=phylum_colors[-8]) +
        geom_jitter(alpha=0.5) + 
        facet_grid(.~type)+
        theme_minimal() + 
        labs(y="Family", x="Relative abundance", color="Phylum")
```

<img src="_main_files/figure-html/taxonomy_jitterplot_family-1.png" width="960" />

### Genus

```{.r .script-source}
genus_summary <- genome_counts_filt %>%
  mutate_at(vars(-genome),~./sum(.)) %>% #apply TSS nornalisation
  pivot_longer(-genome, names_to = "sample", values_to = "count") %>% #reduce to minimum number of columns
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>% #append sample metadata
  left_join(genome_metadata, by = join_by(genome == genome)) %>% #append genome metadata
  group_by(sample,genus) %>%
  summarise(relabun=sum(count)) %>%
  filter(genus != "g__")

genus_summary %>%
    group_by(genus) %>%
    summarise(mean=mean(relabun, na.rm=TRUE),sd=sd(relabun, na.rm=TRUE)) %>%
    arrange(-mean) %>%
    tt()
```

```{=html}
<!DOCTYPE html> 
<html lang="en">
  <head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tinytable_j75duhm0x3w6p4k8bpm7</title>
    <style>
.table td.tinytable_css_ufn0rcy9i83e42rqq5qj, .table th.tinytable_css_ufn0rcy9i83e42rqq5qj {    border-bottom: solid 0.1em #d3d8dc; }
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']]
      },
      svg: {
        fontCache: 'global'
      }
    };
    </script>
  </head>

  <body>
    <div class="container">
      <table class="table table-borderless" id="tinytable_j75duhm0x3w6p4k8bpm7" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <thead>
        
              <tr>
                <th scope="col">genus</th>
                <th scope="col">mean</th>
                <th scope="col">sd</th>
              </tr>
        </thead>
        
        <tbody>
                <tr>
                  <td>g__Bacteroides</td>
                  <td>1.374441e-01</td>
                  <td>0.0923562611</td>
                </tr>
                <tr>
                  <td>g__Parabacteroides</td>
                  <td>9.843371e-02</td>
                  <td>0.0803813676</td>
                </tr>
                <tr>
                  <td>g__Phocaeicola</td>
                  <td>7.109518e-02</td>
                  <td>0.0799972428</td>
                </tr>
                <tr>
                  <td>g__Helicobacter_J</td>
                  <td>3.115738e-02</td>
                  <td>0.0603033642</td>
                </tr>
                <tr>
                  <td>g__Mycoplasmoides</td>
                  <td>3.115302e-02</td>
                  <td>0.0765633496</td>
                </tr>
                <tr>
                  <td>g__Odoribacter</td>
                  <td>2.605667e-02</td>
                  <td>0.0268723184</td>
                </tr>
                <tr>
                  <td>g__Roseburia</td>
                  <td>2.425344e-02</td>
                  <td>0.0567044550</td>
                </tr>
                <tr>
                  <td>g__NHYM01</td>
                  <td>2.332808e-02</td>
                  <td>0.0810376677</td>
                </tr>
                <tr>
                  <td>g__Alistipes</td>
                  <td>2.224698e-02</td>
                  <td>0.0287419149</td>
                </tr>
                <tr>
                  <td>g__Coprobacillus</td>
                  <td>2.070109e-02</td>
                  <td>0.0894282233</td>
                </tr>
                <tr>
                  <td>g__Agrobacterium</td>
                  <td>1.579679e-02</td>
                  <td>0.0779688169</td>
                </tr>
                <tr>
                  <td>g__Akkermansia</td>
                  <td>1.224165e-02</td>
                  <td>0.0317653885</td>
                </tr>
                <tr>
                  <td>g__Fusobacterium_A</td>
                  <td>1.046073e-02</td>
                  <td>0.0286441117</td>
                </tr>
                <tr>
                  <td>g__Kineothrix</td>
                  <td>9.107908e-03</td>
                  <td>0.0416218616</td>
                </tr>
                <tr>
                  <td>g__Proteus</td>
                  <td>8.976683e-03</td>
                  <td>0.0694135711</td>
                </tr>
                <tr>
                  <td>g__Dielma</td>
                  <td>8.687357e-03</td>
                  <td>0.0090713197</td>
                </tr>
                <tr>
                  <td>g__CAG-95</td>
                  <td>8.238073e-03</td>
                  <td>0.0207930753</td>
                </tr>
                <tr>
                  <td>g__JAAYNV01</td>
                  <td>7.265789e-03</td>
                  <td>0.0179169564</td>
                </tr>
                <tr>
                  <td>g__Desulfovibrio</td>
                  <td>7.219276e-03</td>
                  <td>0.0214990147</td>
                </tr>
                <tr>
                  <td>g__UBA866</td>
                  <td>7.016767e-03</td>
                  <td>0.0295145125</td>
                </tr>
                <tr>
                  <td>g__Enterococcus</td>
                  <td>6.966137e-03</td>
                  <td>0.0463712943</td>
                </tr>
                <tr>
                  <td>g__Lactococcus</td>
                  <td>6.600789e-03</td>
                  <td>0.0348230465</td>
                </tr>
                <tr>
                  <td>g__Ureaplasma</td>
                  <td>6.412700e-03</td>
                  <td>0.0139267552</td>
                </tr>
                <tr>
                  <td>g__Parabacteroides_B</td>
                  <td>6.132159e-03</td>
                  <td>0.0101543965</td>
                </tr>
                <tr>
                  <td>g__Lacrimispora</td>
                  <td>6.028411e-03</td>
                  <td>0.0098068179</td>
                </tr>
                <tr>
                  <td>g__CALXRO01</td>
                  <td>5.977964e-03</td>
                  <td>0.0313982647</td>
                </tr>
                <tr>
                  <td>g__Citrobacter</td>
                  <td>5.896711e-03</td>
                  <td>0.0340533686</td>
                </tr>
                <tr>
                  <td>g__NSJ-61</td>
                  <td>5.745781e-03</td>
                  <td>0.0202234473</td>
                </tr>
                <tr>
                  <td>g__Breznakia</td>
                  <td>5.530147e-03</td>
                  <td>0.0240721461</td>
                </tr>
                <tr>
                  <td>g__Clostridium_AQ</td>
                  <td>5.522246e-03</td>
                  <td>0.0123487838</td>
                </tr>
                <tr>
                  <td>g__Bilophila</td>
                  <td>5.044501e-03</td>
                  <td>0.0089558435</td>
                </tr>
                <tr>
                  <td>g__Hungatella_A</td>
                  <td>4.964136e-03</td>
                  <td>0.0096921078</td>
                </tr>
                <tr>
                  <td>g__Escherichia</td>
                  <td>4.342538e-03</td>
                  <td>0.0270859242</td>
                </tr>
                <tr>
                  <td>g__Salmonella</td>
                  <td>4.319018e-03</td>
                  <td>0.0148769561</td>
                </tr>
                <tr>
                  <td>g__UMGS1251</td>
                  <td>4.312965e-03</td>
                  <td>0.0073601071</td>
                </tr>
                <tr>
                  <td>g__MGBC136627</td>
                  <td>4.305492e-03</td>
                  <td>0.0164533523</td>
                </tr>
                <tr>
                  <td>g__Hungatella</td>
                  <td>4.150386e-03</td>
                  <td>0.0194068227</td>
                </tr>
                <tr>
                  <td>g__Clostridium_Q</td>
                  <td>4.146767e-03</td>
                  <td>0.0052575243</td>
                </tr>
                <tr>
                  <td>g__Brevinema</td>
                  <td>4.098862e-03</td>
                  <td>0.0125062564</td>
                </tr>
                <tr>
                  <td>g__Thomasclavelia</td>
                  <td>4.046233e-03</td>
                  <td>0.0110779301</td>
                </tr>
                <tr>
                  <td>g__Scatousia</td>
                  <td>3.752075e-03</td>
                  <td>0.0104403539</td>
                </tr>
                <tr>
                  <td>g__Mailhella</td>
                  <td>3.745039e-03</td>
                  <td>0.0104110785</td>
                </tr>
                <tr>
                  <td>g__Copromonas</td>
                  <td>3.643508e-03</td>
                  <td>0.0050495456</td>
                </tr>
                <tr>
                  <td>g__Enterocloster</td>
                  <td>3.613702e-03</td>
                  <td>0.0047492729</td>
                </tr>
                <tr>
                  <td>g__Ventrimonas</td>
                  <td>3.566172e-03</td>
                  <td>0.0071931788</td>
                </tr>
                <tr>
                  <td>g__Fournierella</td>
                  <td>3.313097e-03</td>
                  <td>0.0063192740</td>
                </tr>
                <tr>
                  <td>g__Limenecus</td>
                  <td>3.230504e-03</td>
                  <td>0.0066725343</td>
                </tr>
                <tr>
                  <td>g__Mucinivorans</td>
                  <td>3.006847e-03</td>
                  <td>0.0379999623</td>
                </tr>
                <tr>
                  <td>g__Lawsonia</td>
                  <td>2.916613e-03</td>
                  <td>0.0103686789</td>
                </tr>
                <tr>
                  <td>g__MGBC133411</td>
                  <td>2.902785e-03</td>
                  <td>0.0074333461</td>
                </tr>
                <tr>
                  <td>g__Caccovivens</td>
                  <td>2.887473e-03</td>
                  <td>0.0112659902</td>
                </tr>
                <tr>
                  <td>g__Sarcina</td>
                  <td>2.855351e-03</td>
                  <td>0.0174153876</td>
                </tr>
                <tr>
                  <td>g__Eisenbergiella</td>
                  <td>2.796704e-03</td>
                  <td>0.0069384489</td>
                </tr>
                <tr>
                  <td>g__Bacteroides_G</td>
                  <td>2.781473e-03</td>
                  <td>0.0352463088</td>
                </tr>
                <tr>
                  <td>g__CAJLXD01</td>
                  <td>2.730769e-03</td>
                  <td>0.0088951735</td>
                </tr>
                <tr>
                  <td>g__Acetatifactor</td>
                  <td>2.654208e-03</td>
                  <td>0.0055286194</td>
                </tr>
                <tr>
                  <td>g__Blautia</td>
                  <td>2.598789e-03</td>
                  <td>0.0062369300</td>
                </tr>
                <tr>
                  <td>g__Velocimicrobium</td>
                  <td>2.235984e-03</td>
                  <td>0.0067748392</td>
                </tr>
                <tr>
                  <td>g__C-19</td>
                  <td>2.235603e-03</td>
                  <td>0.0048296119</td>
                </tr>
                <tr>
                  <td>g__CAZU01</td>
                  <td>2.189719e-03</td>
                  <td>0.0066369837</td>
                </tr>
                <tr>
                  <td>g__Negativibacillus</td>
                  <td>2.145239e-03</td>
                  <td>0.0056002700</td>
                </tr>
                <tr>
                  <td>g__Intestinimonas</td>
                  <td>2.003816e-03</td>
                  <td>0.0035552824</td>
                </tr>
                <tr>
                  <td>g__Rikenella</td>
                  <td>1.998193e-03</td>
                  <td>0.0037323264</td>
                </tr>
                <tr>
                  <td>g__Phascolarctobacterium</td>
                  <td>1.980431e-03</td>
                  <td>0.0051045211</td>
                </tr>
                <tr>
                  <td>g__Butyricimonas</td>
                  <td>1.886974e-03</td>
                  <td>0.0042483569</td>
                </tr>
                <tr>
                  <td>g__RGIG6463</td>
                  <td>1.855727e-03</td>
                  <td>0.0040258495</td>
                </tr>
                <tr>
                  <td>g__MGBC116941</td>
                  <td>1.783700e-03</td>
                  <td>0.0077801927</td>
                </tr>
                <tr>
                  <td>g__JALFVM01</td>
                  <td>1.712574e-03</td>
                  <td>0.0038669765</td>
                </tr>
                <tr>
                  <td>g__Oscillibacter</td>
                  <td>1.546231e-03</td>
                  <td>0.0025273862</td>
                </tr>
                <tr>
                  <td>g__Acinetobacter</td>
                  <td>1.540093e-03</td>
                  <td>0.0099192011</td>
                </tr>
                <tr>
                  <td>g__Pseudoflavonifractor</td>
                  <td>1.489286e-03</td>
                  <td>0.0027026675</td>
                </tr>
                <tr>
                  <td>g__Citrobacter_A</td>
                  <td>1.444197e-03</td>
                  <td>0.0061395639</td>
                </tr>
                <tr>
                  <td>g__Staphylococcus</td>
                  <td>1.411833e-03</td>
                  <td>0.0051727635</td>
                </tr>
                <tr>
                  <td>g__14-2</td>
                  <td>1.228249e-03</td>
                  <td>0.0098038984</td>
                </tr>
                <tr>
                  <td>g__RGIG4733</td>
                  <td>1.225656e-03</td>
                  <td>0.0038271024</td>
                </tr>
                <tr>
                  <td>g__Beduini</td>
                  <td>1.217163e-03</td>
                  <td>0.0025500013</td>
                </tr>
                <tr>
                  <td>g__Scatocola</td>
                  <td>1.162463e-03</td>
                  <td>0.0045748144</td>
                </tr>
                <tr>
                  <td>g__Enterococcus_A</td>
                  <td>1.123893e-03</td>
                  <td>0.0100947603</td>
                </tr>
                <tr>
                  <td>g__UBA1436</td>
                  <td>1.088209e-03</td>
                  <td>0.0062780187</td>
                </tr>
                <tr>
                  <td>g__Faecisoma</td>
                  <td>1.057958e-03</td>
                  <td>0.0056043491</td>
                </tr>
                <tr>
                  <td>g__RGIG9287</td>
                  <td>9.993630e-04</td>
                  <td>0.0094920364</td>
                </tr>
                <tr>
                  <td>g__CAG-345</td>
                  <td>9.840222e-04</td>
                  <td>0.0061275213</td>
                </tr>
                <tr>
                  <td>g__Lachnotalea</td>
                  <td>9.593025e-04</td>
                  <td>0.0033990676</td>
                </tr>
                <tr>
                  <td>g__Blautia_A</td>
                  <td>9.546974e-04</td>
                  <td>0.0029560949</td>
                </tr>
                <tr>
                  <td>g__Ruthenibacterium</td>
                  <td>8.602962e-04</td>
                  <td>0.0024818365</td>
                </tr>
                <tr>
                  <td>g__CAG-269</td>
                  <td>8.276280e-04</td>
                  <td>0.0048017072</td>
                </tr>
                <tr>
                  <td>g__Marseille-P3106</td>
                  <td>8.230475e-04</td>
                  <td>0.0017874580</td>
                </tr>
                <tr>
                  <td>g__WRHT01</td>
                  <td>6.666234e-04</td>
                  <td>0.0027445999</td>
                </tr>
                <tr>
                  <td>g__Eggerthella</td>
                  <td>6.643755e-04</td>
                  <td>0.0021620275</td>
                </tr>
                <tr>
                  <td>g__CHH4-2</td>
                  <td>6.371240e-04</td>
                  <td>0.0020328940</td>
                </tr>
                <tr>
                  <td>g__Corynebacterium</td>
                  <td>6.168531e-04</td>
                  <td>0.0061497354</td>
                </tr>
                <tr>
                  <td>g__Serratia_A</td>
                  <td>6.076344e-04</td>
                  <td>0.0077577570</td>
                </tr>
                <tr>
                  <td>g__Anaerotruncus</td>
                  <td>6.058602e-04</td>
                  <td>0.0016447558</td>
                </tr>
                <tr>
                  <td>g__RUG14156</td>
                  <td>5.735678e-04</td>
                  <td>0.0021869659</td>
                </tr>
                <tr>
                  <td>g__RGIG1896</td>
                  <td>5.683407e-04</td>
                  <td>0.0051791669</td>
                </tr>
                <tr>
                  <td>g__IOR16</td>
                  <td>5.574841e-04</td>
                  <td>0.0016418264</td>
                </tr>
                <tr>
                  <td>g__Faecimonas</td>
                  <td>5.146607e-04</td>
                  <td>0.0054508437</td>
                </tr>
                <tr>
                  <td>g__CAG-56</td>
                  <td>5.096368e-04</td>
                  <td>0.0016613952</td>
                </tr>
                <tr>
                  <td>g__MGBC140009</td>
                  <td>4.851579e-04</td>
                  <td>0.0024491799</td>
                </tr>
                <tr>
                  <td>g__CALURL01</td>
                  <td>4.805911e-04</td>
                  <td>0.0017020401</td>
                </tr>
                <tr>
                  <td>g__Merdimorpha</td>
                  <td>4.705627e-04</td>
                  <td>0.0013078764</td>
                </tr>
                <tr>
                  <td>g__RGIG8482</td>
                  <td>4.560993e-04</td>
                  <td>0.0030287706</td>
                </tr>
                <tr>
                  <td>g__Enterobacter</td>
                  <td>4.223379e-04</td>
                  <td>0.0042068345</td>
                </tr>
                <tr>
                  <td>g__Klebsiella</td>
                  <td>4.203682e-04</td>
                  <td>0.0049802041</td>
                </tr>
                <tr>
                  <td>g__Caccenecus</td>
                  <td>4.086273e-04</td>
                  <td>0.0018112589</td>
                </tr>
                <tr>
                  <td>g__Aliarcobacter</td>
                  <td>3.928587e-04</td>
                  <td>0.0050156837</td>
                </tr>
                <tr>
                  <td>g__Scatenecus</td>
                  <td>3.851876e-04</td>
                  <td>0.0018282510</td>
                </tr>
                <tr>
                  <td>g__Alcaligenes</td>
                  <td>3.835606e-04</td>
                  <td>0.0048969735</td>
                </tr>
                <tr>
                  <td>g__Plesiomonas</td>
                  <td>3.766988e-04</td>
                  <td>0.0027593254</td>
                </tr>
                <tr>
                  <td>g__JAHHSE01</td>
                  <td>3.529590e-04</td>
                  <td>0.0014998851</td>
                </tr>
                <tr>
                  <td>g__HGM05232</td>
                  <td>3.508548e-04</td>
                  <td>0.0009525792</td>
                </tr>
                <tr>
                  <td>g__Enterococcus_B</td>
                  <td>3.475714e-04</td>
                  <td>0.0022665993</td>
                </tr>
                <tr>
                  <td>g__Egerieousia</td>
                  <td>3.295199e-04</td>
                  <td>0.0011408058</td>
                </tr>
                <tr>
                  <td>g__Stoquefichus</td>
                  <td>3.137462e-04</td>
                  <td>0.0020871798</td>
                </tr>
                <tr>
                  <td>g__Hepatoplasma</td>
                  <td>3.099135e-04</td>
                  <td>0.0039567109</td>
                </tr>
                <tr>
                  <td>g__Paracoccus</td>
                  <td>3.068016e-04</td>
                  <td>0.0039169801</td>
                </tr>
                <tr>
                  <td>g__Moheibacter</td>
                  <td>2.873650e-04</td>
                  <td>0.0032049404</td>
                </tr>
                <tr>
                  <td>g__Scatomorpha</td>
                  <td>2.738230e-04</td>
                  <td>0.0010358302</td>
                </tr>
                <tr>
                  <td>g__Emergencia</td>
                  <td>2.601331e-04</td>
                  <td>0.0013298673</td>
                </tr>
                <tr>
                  <td>g__UBA7185</td>
                  <td>2.523935e-04</td>
                  <td>0.0014817660</td>
                </tr>
                <tr>
                  <td>g__Eubacterium</td>
                  <td>1.707442e-04</td>
                  <td>0.0006844943</td>
                </tr>
                <tr>
                  <td>g__Sphingobacterium</td>
                  <td>1.561202e-04</td>
                  <td>0.0012685229</td>
                </tr>
                <tr>
                  <td>g__Devosia</td>
                  <td>1.544841e-04</td>
                  <td>0.0015368528</td>
                </tr>
                <tr>
                  <td>g__Anaerosporobacter</td>
                  <td>1.507638e-04</td>
                  <td>0.0012978048</td>
                </tr>
                <tr>
                  <td>g__Caccomorpha</td>
                  <td>1.434035e-04</td>
                  <td>0.0010730603</td>
                </tr>
                <tr>
                  <td>g__UBA2658</td>
                  <td>1.355578e-04</td>
                  <td>0.0007332702</td>
                </tr>
                <tr>
                  <td>g__Protoclostridium</td>
                  <td>1.324439e-04</td>
                  <td>0.0007783049</td>
                </tr>
                <tr>
                  <td>g__Angelakisella</td>
                  <td>1.315171e-04</td>
                  <td>0.0009387427</td>
                </tr>
                <tr>
                  <td>g__Cetobacterium_A</td>
                  <td>9.652924e-05</td>
                  <td>0.0008876688</td>
                </tr>
                <tr>
                  <td>g__Rahnella</td>
                  <td>6.708891e-05</td>
                  <td>0.0008565338</td>
                </tr>
                <tr>
                  <td>g__Peptostreptococcus</td>
                  <td>2.371535e-05</td>
                  <td>0.0003027773</td>
                </tr>
        </tbody>
      </table>
    </div>

    <script>
      function styleCell_tinytable_m9fjtv88qo6tgqya90gq(i, j, css_id) {
        var table = document.getElementById("tinytable_j75duhm0x3w6p4k8bpm7");
        table.rows[i].cells[j].classList.add(css_id);
      }
      function insertSpanRow(i, colspan, content) {
        var table = document.getElementById('tinytable_j75duhm0x3w6p4k8bpm7');
        var newRow = table.insertRow(i);
        var newCell = newRow.insertCell(0);
        newCell.setAttribute("colspan", colspan);
        // newCell.innerText = content;
        // this may be unsafe, but innerText does not interpret <br>
        newCell.innerHTML = content;
      }
      function spanCell_tinytable_m9fjtv88qo6tgqya90gq(i, j, rowspan, colspan) {
        var table = document.getElementById("tinytable_j75duhm0x3w6p4k8bpm7");
        const targetRow = table.rows[i];
        const targetCell = targetRow.cells[j];
        for (let r = 0; r < rowspan; r++) {
          // Only start deleting cells to the right for the first row (r == 0)
          if (r === 0) {
            // Delete cells to the right of the target cell in the first row
            for (let c = colspan - 1; c > 0; c--) {
              if (table.rows[i + r].cells[j + c]) {
                table.rows[i + r].deleteCell(j + c);
              }
            }
          }
          // For rows below the first, delete starting from the target column
          if (r > 0) {
            for (let c = colspan - 1; c >= 0; c--) {
              if (table.rows[i + r] && table.rows[i + r].cells[j]) {
                table.rows[i + r].deleteCell(j);
              }
            }
          }
        }
        // Set rowspan and colspan of the target cell
        targetCell.rowSpan = rowspan;
        targetCell.colSpan = colspan;
      }

window.addEventListener('load', function () { styleCell_tinytable_m9fjtv88qo6tgqya90gq(0, 0, 'tinytable_css_ufn0rcy9i83e42rqq5qj') })
window.addEventListener('load', function () { styleCell_tinytable_m9fjtv88qo6tgqya90gq(0, 1, 'tinytable_css_ufn0rcy9i83e42rqq5qj') })
window.addEventListener('load', function () { styleCell_tinytable_m9fjtv88qo6tgqya90gq(0, 2, 'tinytable_css_ufn0rcy9i83e42rqq5qj') })
    </script>

  </body>

</html>
```


```{.r .script-source}
genus_arrange <- genus_summary %>%
    group_by(genus) %>%
    summarise(mean=sum(relabun)) %>%
    filter(genus != "g__")%>%
    arrange(-mean) %>%
    select(genus) %>%
    mutate(genus= sub("^g__", "", genus)) %>%
    pull()

genus_summary %>%
    left_join(genome_metadata %>% select(genus,phylum) %>% unique(),by=join_by(genus==genus)) %>%
    left_join(sample_metadata,by=join_by(sample==Tube_code)) %>%
    mutate(genus= sub("^g__", "", genus)) %>%
    filter(genus %in% genus_arrange[1:20]) %>%
    mutate(genus=factor(genus,levels=rev(genus_arrange[1:20]))) %>%
    filter(relabun > 0) %>%
    ggplot(aes(x=relabun, y=genus, group=genus, color=phylum)) +
        scale_color_manual(values=phylum_colors[-c(3,4,6,8)]) +
        geom_jitter(alpha=0.5) + 
        facet_grid(.~type)+
        theme_minimal() + 
        labs(y="Family", x="Relative abundance", color="Phylum")
```

<img src="_main_files/figure-html/taxonomy_jitterplot_genus-1.png" width="960" />

# Diversity analysis

## Alpha diversity


```{.r .script-source}
# Calculate Hill numbers
richness <- genome_counts_filt %>%
  column_to_rownames(var = "genome") %>%
  dplyr::select(where(~ !all(. == 0))) %>%
  hilldiv(., q = 0) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(richness = 1) %>%
  rownames_to_column(var = "sample")

neutral <- genome_counts_filt %>%
  column_to_rownames(var = "genome") %>%
  dplyr::select(where(~ !all(. == 0))) %>%
  hilldiv(., q = 1) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(neutral = 1) %>%
  rownames_to_column(var = "sample")

phylogenetic <- genome_counts_filt %>%
  column_to_rownames(var = "genome") %>%
  dplyr::select(where(~ !all(. == 0))) %>%
  hilldiv(., q = 1, tree = genome_tree) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(phylogenetic = 1) %>%
  rownames_to_column(var = "sample")

# Aggregate basal GIFT into elements
dist <- genome_gifts %>%
  to.elements(., GIFT_db) %>%
  traits2dist(., method = "gower")

functional <- genome_counts_filt %>%
  column_to_rownames(var = "genome") %>%
  dplyr::select(where(~ !all(. == 0))) %>%
  hilldiv(., q = 1, dist = dist) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(functional = 1) %>%
  rownames_to_column(var = "sample") %>%
  mutate(functional = if_else(is.nan(functional), 1, functional))

# Merge all metrics
alpha_div <- richness %>%
  full_join(neutral, by = join_by(sample == sample)) %>%
  full_join(phylogenetic, by = join_by(sample == sample)) %>%
  full_join(functional, by = join_by(sample == sample))
```

### Wild samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="0_Wild") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = Population, group=Population, color=Population, fill=Population)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Population",
          breaks=c("Cold_wet","Hot_dry"),
          labels=c("Cold","Hot"),
          values=c('#008080', "#d57d2c")) +
      scale_fill_manual(name="Population",
          breaks=c("Cold_wet","Hot_dry"),
          labels=c("Cold","Hot"),
          values=c('#00808050', "#d57d2c50")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.58) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_0-1.png" width="960" />

### Acclimation samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="1_Acclimation") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = Population, group=Population, color=Population, fill=Population)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Population",
          breaks=c("Cold_wet","Hot_dry"),
          labels=c("Cold","Hot"),
          values=c('#008080', "#d57d2c")) +
      scale_fill_manual(name="Population",
          breaks=c("Cold_wet","Hot_dry"),
          labels=c("Cold","Hot"),
          values=c('#00808050', "#d57d2c50")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.58) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_1-1.png" width="960" />

### Antibiotics samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="2_Antibiotics") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = Population, group=Population, color=Population, fill=Population)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Population",
          breaks=c("Cold_wet","Hot_dry"),
          labels=c("Cold","Hot"),
          values=c('#008080', "#d57d2c")) +
      scale_fill_manual(name="Population",
          breaks=c("Cold_wet","Hot_dry"),
          labels=c("Cold","Hot"),
          values=c('#00808050', "#d57d2c50")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.58) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_2-1.png" width="960" />

### Transplant_1 samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="3_Transplant1") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = type, group=type, color=type, fill=type)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c", "#76b183")) +
      scale_fill_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.7) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_3-1.png" width="960" />

### Transplant_2 samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="4_Transplant2") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = type, group=type, color=type, fill=type)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c", "#76b183")) +
      scale_fill_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.7) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_4-1.png" width="960" />

### Post-Transplant_1 samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="5_Post-FMT1") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = type, group=type, color=type, fill=type)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c", "#76b183")) +
      scale_fill_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.7) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_5-1.png" width="960" />

### Post-Transplant_2 samples


```{.r .script-source}
alpha_div %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(., sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="6_Post-FMT2") %>%
  mutate(metric=factor(metric,levels=c("richness","neutral","phylogenetic","functional"))) %>%
      ggplot(aes(y = value, x = type, group=type, color=type, fill=type)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(alpha=0.5) +
      scale_color_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c", "#76b183")) +
      scale_fill_manual(name="Type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
      facet_wrap(. ~ metric,scales = "free") +
      coord_cartesian(xlim = c(1, NA)) +
      stat_compare_means(size=3, label.x=.7) +
      theme_classic() +
        theme(
    strip.background = element_blank(),
    panel.grid.minor.x = element_line(size = .1, color = "grey"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # Increase plot size
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
      ) +
  ylab("Alpha diversity")
```

<img src="_main_files/figure-html/alpha_div_plot_6-1.png" width="960" />

## Beta diversity


```{.r .script-source}
beta_q0n <- genome_counts_filt %>%
  column_to_rownames(., "genome") %>%
  hillpair(., q = 0)

beta_q1n <- genome_counts_filt %>%
  column_to_rownames(., "genome") %>%
  hillpair(., q = 1)

beta_q1p <- genome_counts_filt %>%
  column_to_rownames(., "genome") %>%
  hillpair(., q = 1, tree = genome_tree)

beta_q1f <- genome_counts_filt %>%
  column_to_rownames(., "genome") %>%
  hillpair(., q = 1, dist = dist)
```

<!--chapter:end:04_community_composition.Rmd-->

## Permanovas

#### Load required data


```{.r .script-source}
meta <- column_to_rownames(sample_metadata, "Tube_code")
```

### 1. Are the wild populations similar?

#### Wild: P.muralis vs P.liolepis


```{.r .script-source}
wild <- meta %>%
  filter(time_point == "0_Wild")

# Create a temporary modified version of genome_counts_filt
temp_genome_counts <- transform(genome_counts_filt, row.names = genome_counts_filt$genome)
temp_genome_counts$genome <- NULL

wild.counts <- temp_genome_counts[, which(colnames(temp_genome_counts) %in% rownames(wild))]
identical(sort(colnames(wild.counts)), sort(as.character(rownames(wild))))

wild_nmds <- sample_metadata %>%
  filter(time_point == "0_Wild")
```

#### Number of samples used


```{.script-output}
[1] 27
```


```{.r .script-source}
beta_div_richness_wild<-hillpair(data=wild.counts, q=0)
beta_div_neutral_wild<-hillpair(data=wild.counts, q=1)
beta_div_phylo_wild<-hillpair(data=wild.counts, q=1, tree=genome_tree)
beta_div_func_wild<-hillpair(data=wild.counts, q=1, dist=dist)
```

#### Richness


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.000012 0.000012 0.0012    999  0.975
Residuals 25 0.257281 0.010291                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
                  Podarcis_liolepis Podarcis_muralis
Podarcis_liolepis                              0.979
Podarcis_muralis            0.97302                 
```


|         | Df| SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|--------:|---------:|--------:|------:|
|species  |  1| 1.542719| 0.2095041| 6.625717|  0.001|
|Residual | 25| 5.820951| 0.7904959|       NA|     NA|
|Total    | 26| 7.363669| 1.0000000|       NA|     NA|


#### Neutral


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     1 0.000048 0.0000476 0.0044    999  0.948
Residuals 25 0.270114 0.0108046                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
                  Podarcis_liolepis Podarcis_muralis
Podarcis_liolepis                              0.946
Podarcis_muralis            0.94763                 
```


|         | Df| SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|--------:|---------:|--------:|------:|
|species  |  1| 1.918266| 0.2608511| 8.822682|  0.001|
|Residual | 25| 5.435610| 0.7391489|       NA|     NA|
|Total    | 26| 7.353876| 1.0000000|       NA|     NA|

#### Phylogenetic


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.03585 0.035847 2.4912    999  0.135
Residuals 25 0.35973 0.014389                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
                  Podarcis_liolepis Podarcis_muralis
Podarcis_liolepis                              0.133
Podarcis_muralis            0.12705                 
```


|         | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|---------:|---------:|--------:|------:|
|species  |  1| 0.3218613| 0.2162815| 6.899207|  0.001|
|Residual | 25| 1.1662981| 0.7837185|       NA|     NA|
|Total    | 26| 1.4881594| 1.0000000|       NA|     NA|

#### Functional


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.018367 0.018367 1.5597    999  0.203
Residuals 25 0.294402 0.011776                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
                  Podarcis_liolepis Podarcis_muralis
Podarcis_liolepis                              0.209
Podarcis_muralis            0.22328                 
```


|         | Df|  SumOfSqs|       R2|        F| Pr(>F)|
|:--------|--:|---------:|--------:|--------:|------:|
|species  |  1| 0.0858578| 0.172879| 5.225323|  0.053|
|Residual | 25| 0.4107775| 0.827121|       NA|     NA|
|Total    | 26| 0.4966352| 1.000000|       NA|     NA|


```{.r .script-source}
beta_q0n_nmds_wild <- beta_div_richness_wild$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace=FALSE) %>%
				vegan::scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(wild_nmds, by = join_by(sample == Tube_code))

beta_q1n_nmds_wild <- beta_div_neutral_wild$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace=FALSE) %>%
				vegan::scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(wild_nmds, by = join_by(sample == Tube_code))

beta_q1p_nmds_wild <- beta_div_phylo_wild$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				vegan::scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(wild_nmds, by = join_by(sample == Tube_code))

beta_q1f_nmds_wild <- beta_div_func_wild$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				vegan::scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(wild_nmds, by = join_by(sample == Tube_code))
```



<img src="_main_files/figure-html/beta_neutral_nmds_plot_final-1.png" width="960" />

### 2. Effect of acclimation


```{.r .script-source}
accli <- meta %>%
  filter(time_point == "1_Acclimation")

# Create a temporary modified version of genome_counts_filt
temp_genome_counts <- transform(genome_counts_filt, row.names = genome_counts_filt$genome)
temp_genome_counts$genome <- NULL

accli.counts <- temp_genome_counts[, which(colnames(temp_genome_counts) %in% rownames(accli))]
identical(sort(colnames(accli.counts)), sort(as.character(rownames(accli))))

accli_nmds <- sample_metadata %>%
  filter(time_point == "1_Acclimation")
```

#### Number of samples used


```{.script-output}
[1] 27
```


```{.r .script-source}
beta_div_richness_accli<-hillpair(data=accli.counts, q=0)
beta_div_neutral_accli<-hillpair(data=accli.counts, q=1)
beta_div_phylo_accli<-hillpair(data=accli.counts, q=1, tree=genome_tree)
beta_div_func_accli<-hillpair(data=accli.counts, q=1, dist=dist)
```

#### Richness


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
Groups     1 0.11796 0.117959 12.963    999  0.003 **
Residuals 25 0.22748 0.009099                        
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
          Cold_wet Hot_dry
Cold_wet             0.002
Hot_dry  0.0013711        
```


|           | Df| SumOfSqs|       R2|        F| Pr(>F)|
|:----------|--:|--------:|--------:|--------:|------:|
|Population |  1| 1.639807| 0.179834| 5.481634|  0.001|
|Residual   | 25| 7.478640| 0.820166|       NA|     NA|
|Total      | 26| 9.118447| 1.000000|       NA|     NA|


#### Neutral


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
Groups     1 0.07844 0.078443 5.2384    999  0.032 *
Residuals 25 0.37437 0.014975                       
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
         Cold_wet Hot_dry
Cold_wet            0.029
Hot_dry  0.030815        
```


|           | Df| SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|--------:|---------:|--------:|------:|
|Population |  1| 1.947003| 0.2306127| 7.493387|  0.001|
|Residual   | 25| 6.495736| 0.7693873|       NA|     NA|
|Total      | 26| 8.442739| 1.0000000|       NA|     NA|

#### Phylogenetic


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.06739 0.067395 2.9532    999  0.104
Residuals 25 0.57052 0.022821                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
         Cold_wet Hot_dry
Cold_wet            0.095
Hot_dry  0.098068        
```


|           | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|---------:|---------:|--------:|------:|
|Population |  1| 0.2441653| 0.1224638| 3.488854|  0.032|
|Residual   | 25| 1.7496100| 0.8775362|       NA|     NA|
|Total      | 26| 1.9937754| 1.0000000|       NA|     NA|

#### Functional


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.02496 0.024955 0.6729    999  0.454
Residuals 25 0.92714 0.037085                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
         Cold_wet Hot_dry
Cold_wet            0.448
Hot_dry   0.41979        
```


|           | Df|  SumOfSqs|        R2|         F| Pr(>F)|
|:----------|--:|---------:|---------:|---------:|------:|
|Population |  1| 0.0279454| 0.0248037| 0.6358634|  0.466|
|Residual   | 25| 1.0987171| 0.9751963|        NA|     NA|
|Total      | 26| 1.1266624| 1.0000000|        NA|     NA|



```{.r .script-source}
beta_q0n_nmds_accli <- beta_div_richness_accli$S %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace=FALSE) %>%
  vegan::scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(accli_nmds, by = join_by(sample == Tube_code))

beta_q1n_nmds_accli <- beta_div_neutral_accli$S %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE, trace=FALSE) %>%
  vegan::scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(accli_nmds, by = join_by(sample == Tube_code))

beta_q1p_nmds_accli <- beta_div_phylo_accli$S %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
  vegan::scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(accli_nmds, by = join_by(sample == Tube_code))

beta_q1f_nmds_accli <- beta_div_func_accli$S %>%
  metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
  vegan::scores() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(accli_nmds, by = join_by(sample == Tube_code))
```



<img src="_main_files/figure-html/beta_neutral_nmds_plot_final_accli-1.png" width="960" />

### 3. Comparison between Wild and Acclimation


```{.r .script-source}
accli1 <- meta  %>%
  filter(time_point == "0_Wild" | time_point == "1_Acclimation")

temp_genome_counts <- transform(genome_counts_filt, row.names = genome_counts_filt$genome)
temp_genome_counts$genome <- NULL

accli1.counts <- temp_genome_counts[,which(colnames(temp_genome_counts) %in% rownames(accli1))]
identical(sort(colnames(accli1.counts)),sort(as.character(rownames(accli1))))

accli1_nmds <- sample_metadata %>%
  filter(time_point == "0_Wild" | time_point == "1_Acclimation")
```

#### Number of samples used


```{.script-output}
[1] 54
```


```{.r .script-source}
beta_div_richness_accli1<-hillpair(data=accli1.counts, q=0)
beta_div_neutral_accli1<-hillpair(data=accli1.counts, q=1)
beta_div_phylo_accli1<-hillpair(data=accli1.counts, q=1, tree=genome_tree)
beta_div_func_accli1<-hillpair(data=accli1.counts, q=1, dist=dist)
```

##### Richness


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
Groups     1 0.05014 0.050145 6.2252    999  0.011 *
Residuals 52 0.41886 0.008055                       
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
                0_Wild 1_Acclimation
0_Wild                         0.014
1_Acclimation 0.015808              
```


|           | Df|   SumOfSqs|        R2|         F| Pr(>F)|
|:----------|--:|----------:|---------:|---------:|------:|
|time_point |  1|  0.6172653| 0.0360987|  3.933397|  0.001|
|species    |  1|  2.8279677| 0.1653842| 18.020647|  0.001|
|individual | 25|  9.5739861| 0.5599025|  2.440331|  0.001|
|Residual   | 26|  4.0801621| 0.2386146|        NA|     NA|
|Total      | 53| 17.0993812| 1.0000000|        NA|     NA|

##### Neutral


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     1 0.0199 0.0199035 2.1213    999  0.139
Residuals 52 0.4879 0.0093827                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
               0_Wild 1_Acclimation
0_Wild                        0.142
1_Acclimation 0.15128              
```


|           | Df|   SumOfSqs|        R2|         F| Pr(>F)|
|:----------|--:|----------:|---------:|---------:|------:|
|time_point |  1|  0.9050519| 0.0541893|  6.651487|  0.001|
|species    |  1|  3.3236300| 0.1989999| 24.426315|  0.001|
|individual | 25|  8.9352276| 0.5349902|  2.626702|  0.001|
|Residual   | 26|  3.5377576| 0.2118206|        NA|     NA|
|Total      | 53| 16.7016671| 1.0000000|        NA|     NA|


##### Phylogenetic


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.01334 0.013340 0.6524    999  0.435
Residuals 52 1.06332 0.020449                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
               0_Wild 1_Acclimation
0_Wild                        0.437
1_Acclimation 0.42294              
```


|           | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|---------:|---------:|--------:|------:|
|time_point |  1| 0.2890434| 0.0766494| 7.532050|  0.001|
|species    |  1| 0.3508889| 0.0930498| 9.143655|  0.001|
|individual | 25| 2.1332925| 0.5657133| 2.223620|  0.001|
|Residual   | 26| 0.9977533| 0.2645874|       NA|     NA|
|Total      | 53| 3.7709782| 1.0000000|       NA|     NA|

##### Functional


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq  Mean Sq      F N.Perm Pr(>F)
Groups     1 0.0123 0.012300 0.4817    999  0.482
Residuals 52 1.3277 0.025533                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
               0_Wild 1_Acclimation
0_Wild                        0.487
1_Acclimation 0.49073              
```


|           | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|---------:|---------:|--------:|------:|
|time_point |  1| 0.0448774| 0.0269021| 2.355512|  0.163|
|species    |  1| 0.0973005| 0.0583275| 5.107077|  0.034|
|individual | 25| 1.0306426| 0.6178264| 2.163841|  0.055|
|Residual   | 26| 0.4953546| 0.2969440|       NA|     NA|
|Total      | 53| 1.6681751| 1.0000000|       NA|     NA|



```{.r .script-source}
beta_richness_nmds_accli1 <- beta_div_richness_accli1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(accli1_nmds, by = c("sample" = "Tube_code"))

beta_neutral_nmds_accli1 <- beta_div_neutral_accli1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(accli1_nmds, by = c("sample" = "Tube_code"))

beta_phylo_nmds_accli1 <- beta_div_phylo_accli1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(accli1_nmds, by = join_by(sample == Tube_code))

beta_func_nmds_accli1 <- beta_div_func_accli1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(accli1_nmds, by = join_by(sample == Tube_code))
```



<img src="_main_files/figure-html/beta_neutral_nmds_plot_accli1_final-1.png" width="960" />

### 4. Do the antibiotics work?

#### Acclimation vs antibiotics


```{.r .script-source}
treat <- meta  %>%
  filter(time_point == "1_Acclimation" | time_point == "2_Antibiotics")

temp_genome_counts <- transform(genome_counts_filt, row.names = genome_counts_filt$genome)
temp_genome_counts$genome <- NULL

treat.counts <- temp_genome_counts[,which(colnames(temp_genome_counts) %in% rownames(treat))]
identical(sort(colnames(treat.counts)),sort(as.character(rownames(treat))))

treat_nmds <- sample_metadata %>%
  filter(time_point == "1_Acclimation" | time_point == "2_Antibiotics")
```

#### Number of samples used


```{.script-output}
[1] 50
```


```{.r .script-source}
beta_div_richness_treat<-hillpair(data=treat.counts, q=0)
beta_div_neutral_treat<-hillpair(data=treat.counts, q=1)
beta_div_phylo_treat<-hillpair(data=treat.counts, q=1, tree=genome_tree)
beta_div_func_treat<-hillpair(data=treat.counts, q=1, dist=dist)
```

##### Richness


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)  
Groups     1 0.025318 0.0253178 6.021    999  0.023 *
Residuals 48 0.201837 0.0042049                      
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
              1_Acclimation 2_Antibiotics
1_Acclimation                       0.027
2_Antibiotics      0.017817              
```


|           | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|---------:|---------:|--------:|------:|
|time_point |  1|  1.888584| 0.0949462| 6.361098|  0.001|
|species    |  1|  2.117109| 0.1064350| 7.130814|  0.001|
|individual | 25|  9.353701| 0.4702455| 1.260199|  0.005|
|Residual   | 22|  6.531709| 0.3283734|       NA|     NA|
|Total      | 49| 19.891103| 1.0000000|       NA|     NA|

##### Neutral


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)  
Groups     1 0.039587 0.039587 6.8387    999  0.013 *
Residuals 48 0.277854 0.005789                       
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
              1_Acclimation 2_Antibiotics
1_Acclimation                       0.015
2_Antibiotics      0.011886              
```


|           | Df|  SumOfSqs|        R2|         F| Pr(>F)|
|:----------|--:|---------:|---------:|---------:|------:|
|time_point |  1|  2.024181| 0.1063620|  9.051981|  0.001|
|species    |  1|  2.853103| 0.1499183| 12.758858|  0.001|
|individual | 25|  9.234189| 0.4852168|  1.651783|  0.001|
|Residual   | 22|  4.919584| 0.2585029|        NA|     NA|
|Total      | 49| 19.031057| 1.0000000|        NA|     NA|


##### Phylogenetic


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
Groups     1 0.58372 0.58372 35.413    999  0.001 ***
Residuals 48 0.79119 0.01648                         
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
              1_Acclimation 2_Antibiotics
1_Acclimation                       0.001
2_Antibiotics    2.9795e-07              
```


|           | Df|  SumOfSqs|        R2|         F| Pr(>F)|
|:----------|--:|---------:|---------:|---------:|------:|
|time_point |  1| 1.8065206| 0.2113909| 18.636551|  0.001|
|species    |  1| 0.7903334| 0.0924813|  8.153292|  0.001|
|individual | 25| 3.8164689| 0.4465860|  1.574869|  0.005|
|Residual   | 22| 2.1325541| 0.2495419|        NA|     NA|
|Total      | 49| 8.5458771| 1.0000000|        NA|     NA|

##### Functional


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
Groups     1 0.18591 0.185914 5.0679    999  0.023 *
Residuals 48 1.76088 0.036685                       
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
              1_Acclimation 2_Antibiotics
1_Acclimation                       0.032
2_Antibiotics      0.028989              
```


|           | Df|  SumOfSqs|        R2|          F| Pr(>F)|
|:----------|--:|---------:|---------:|----------:|------:|
|time_point |  1| 1.8020952| 0.3750193| 33.6195614|  0.001|
|species    |  1| 0.0031247| 0.0006503|  0.0582938|  0.848|
|individual | 25| 1.8208629| 0.3789249|  1.3587875|  0.222|
|Residual   | 22| 1.1792568| 0.2454055|         NA|     NA|
|Total      | 49| 4.8053396| 1.0000000|         NA|     NA|



```{.r .script-source}
beta_richness_nmds_treat <- beta_div_richness_treat$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(treat_nmds, by = c("sample" = "Tube_code"))

beta_neutral_nmds_treat <- beta_div_neutral_treat$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(treat_nmds, by = c("sample" = "Tube_code"))

beta_phylo_nmds_treat <- beta_div_phylo_treat$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(treat_nmds, by = join_by(sample == Tube_code))

beta_func_nmds_treat <- beta_div_func_treat$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(treat_nmds, by = join_by(sample == Tube_code))
```



<img src="_main_files/figure-html/beta_neutral_nmds_plot_treat_final-1.png" width="960" />


### 5. Does the FMT work?

#### Comparison between FMT2 vs Post-FMT2


```{.r .script-source}
#Create newID to identify duplicated samples
transplants_metadata<-sample_metadata%>%
  mutate(Tube_code=str_remove_all(Tube_code, "_a"))
transplants_metadata$newID <- paste(transplants_metadata$Tube_code, "_", transplants_metadata$individual)

transplant3<-transplants_metadata%>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2")%>%
  column_to_rownames("newID")

transplant3_nmds <- transplants_metadata %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2")

full_counts<-temp_genome_counts %>%
    t()%>%
    as.data.frame()%>%
    rownames_to_column("Tube_code")%>%
    full_join(transplants_metadata,by = join_by(Tube_code == Tube_code))

transplant3_counts<-full_counts %>%
  filter(time_point == "4_Transplant2" | time_point == "6_Post-FMT2") %>%
  subset(select=-c(315:324)) %>%
  column_to_rownames("newID")%>%
  subset(select=-c(1))%>%
  t() %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)

identical(sort(colnames(transplant3_counts)),sort(as.character(rownames(transplant3))))
```

#### Number of samples used


```{.script-output}
[1] 49
```



```{.r .script-source}
beta_div_richness_transplant3<-hillpair(data=transplant3_counts, q=0)
beta_div_neutral_transplant3<-hillpair(data=transplant3_counts, q=1)
beta_div_phylo_transplant3<-hillpair(data=transplant3_counts, q=1, tree=genome_tree)
beta_div_func_transplant3<-hillpair(data=transplant3_counts, q=1, dist=dist)

#Arrange of metadata dataframe
transplant3_arrange<-transplant3[labels(beta_div_neutral_transplant3$S),]
```

##### Richness


|           | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|---------:|---------:|--------:|------:|
|species    |  1|  1.180473| 0.0855095| 6.984555|  0.001|
|time_point |  1|  0.860906| 0.0623612| 5.093759|  0.001|
|type       |  1|  1.459433| 0.1057165| 8.635089|  0.001|
|individual | 24|  6.755100| 0.4893170| 1.665341|  0.001|
|Residual   | 21|  3.549250| 0.2570959|       NA|     NA|
|Total      | 48| 13.805162| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
1     Control vs Treatment  1 1.4169018 5.739828 0.15622903   0.001      0.003   *
2   Control vs Hot_control  1 2.0940966 8.509112 0.21005427   0.001      0.003   *
3 Treatment vs Hot_control  1 0.3004618 1.265034 0.04179854   0.159      0.477    
```

##### Neutral


|           | Df|   SumOfSqs|        R2|         F| Pr(>F)|
|:----------|--:|----------:|---------:|---------:|------:|
|species    |  1|  1.2800927| 0.0939787|  8.796453|  0.001|
|time_point |  1|  0.9350566| 0.0686477|  6.425458|  0.001|
|type       |  1|  1.9135997| 0.1404879| 13.149743|  0.001|
|individual | 24|  6.4363516| 0.4725281|  1.842870|  0.001|
|Residual   | 21|  3.0559984| 0.2243577|        NA|     NA|
|Total      | 48| 13.6210990| 1.0000000|        NA|     NA|


```{.script-output}
                     pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
1     Control vs Treatment  1 1.8758788  8.282671 0.21084796   0.001      0.003   *
2   Control vs Hot_control  1 2.4396317 10.635546 0.24945256   0.001      0.003   *
3 Treatment vs Hot_control  1 0.3158428  1.394345 0.04587515   0.126      0.378    
```

##### Phylogenetic


|           | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:----------|--:|---------:|---------:|--------:|------:|
|species    |  1| 0.1400466| 0.0952654| 6.956436|  0.002|
|time_point |  1| 0.1138047| 0.0774145| 5.652939|  0.001|
|type       |  1| 0.1432667| 0.0974558| 7.116383|  0.001|
|individual | 24| 0.6501795| 0.4422784| 1.345663|  0.041|
|Residual   | 21| 0.4227709| 0.2875859|       NA|     NA|
|Total      | 48| 1.4700683| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df  SumsOfSqs  F.Model         R2 p.value p.adjusted sig
1     Control vs Treatment  1 0.14387705 5.735321 0.15612552   0.001      0.003   *
2   Control vs Hot_control  1 0.22715701 9.044894 0.22036587   0.001      0.003   *
3 Treatment vs Hot_control  1 0.04648319 1.704277 0.05550617   0.129      0.387    
```

##### Functional


|           | Df|   SumOfSqs|         R2|          F| Pr(>F)|
|:----------|--:|----------:|----------:|----------:|------:|
|species    |  1|  0.0092808|  0.0077189|  0.4182529|  0.493|
|time_point |  1| -0.0061674| -0.0051295| -0.2779456|  0.895|
|type       |  1|  0.0831052|  0.0691191|  3.7452726|  0.093|
|individual | 24|  0.6501528|  0.5407359|  1.2208414|  0.345|
|Residual   | 21|  0.4659767|  0.3875556|         NA|     NA|
|Total      | 48|  1.2023481|  1.0000000|         NA|     NA|


```{.script-output}
                     pairs Df    SumsOfSqs     F.Model           R2 p.value p.adjusted sig
1     Control vs Treatment  1  0.078539743  4.59293783  0.129040706   0.069      0.207    
2   Control vs Hot_control  1  0.052468954  2.13675422  0.062593948   0.165      0.495    
3 Treatment vs Hot_control  1 -0.002340352 -0.07432315 -0.002569452   0.887      1.000    
```


```{.r .script-source}
beta_richness_nmds_transplant3 <- beta_div_richness_transplant3$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(transplant3_nmds, by = join_by(sample == newID))

beta_neutral_nmds_transplant3 <- beta_div_neutral_transplant3$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(transplant3_nmds, by = join_by(sample == newID))

beta_phylo_nmds_transplant3 <- beta_div_phylo_transplant3$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(transplant3_nmds, by = join_by(sample == newID))

beta_func_nmds_transplant3 <- beta_div_func_transplant3$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(transplant3_nmds, by = join_by(sample == newID))
```


```{.r .script-source}
p0<-beta_richness_nmds_transplant3 %>%
			group_by(individual) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
        scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y = "NMDS2", x="NMDS1 \n Richness beta diversity") +
				theme_classic() +
				theme(legend.position="none")

p1<-beta_neutral_nmds_transplant3 %>%
			group_by(individual) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
        scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
				theme_classic() +
				theme(legend.position="none")
  
p2<-beta_phylo_nmds_transplant3 %>%
			group_by(individual) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
				scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y= element_blank (), x="NMDS1 \n Phylogenetic beta diversity") +
				theme_classic() +
				theme(legend.position="none")

p3<-beta_func_nmds_transplant3 %>%
			group_by(individual) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type, shape=time_point)) +
				scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y= element_blank (), x="NMDS1 \n Functional beta diversity") +
				theme_classic()+
				theme(legend.position="none")
```

<img src="_main_files/figure-html/beta_neutral_nmds_plot_transplant3_final-1.png" width="960" />



#### Comparison between the different experimental time points
##### Comparison against Wild samples


```
The estimated time for calculating the 2850 pairwise combinations is 14 seconds.
```






```{.r .script-source}
ggarrange( p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="right")
```

<img src="_main_files/figure-html/diss_plot1_1, -1.png" width="768" />

##### Comparison against Acclimation samples


```
The estimated time for calculating the 2850 pairwise combinations is 12 seconds.
```







```{.r .script-source}
ggarrange( p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="right")
```

<img src="_main_files/figure-html/diss_plot2_!, -1.png" width="768" />

##### Comparison against Acclimation samples with combined Transplant samples




```
The estimated time for calculating the 5151 pairwise combinations is 22 seconds.
```




```{.r .script-source}
ggarrange( p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="right")
```

<img src="_main_files/figure-html/diss_plot3, -1.png" width="768" />


### 6. Are there differences between the control and the treatment group?

#### After 1 week --> Post-FMT1


```{.r .script-source}
post1 <- meta %>%
  filter(time_point == "5_Post-FMT1")

post1.counts <- temp_genome_counts[,which(colnames(temp_genome_counts) %in% rownames(post1))]
identical(sort(colnames(post1.counts)),sort(as.character(rownames(post1))))

post1_nmds <- sample_metadata %>%
  filter(time_point == "5_Post-FMT1")
```

#### Number of samples used


```{.script-output}
[1] 26
```


```{.r .script-source}
beta_div_richness_post1<-hillpair(data=post1.counts, q=0)
beta_div_neutral_post1<-hillpair(data=post1.counts, q=1)
beta_div_phylo_post1<-hillpair(data=post1.counts, q=1, tree=genome_tree)
beta_div_func_post1<-hillpair(data=post1.counts, q=1, dist=dist)

#Arrange of metadata dataframe
post1_arrange<-post1[labels(beta_div_neutral_post1$S),]
```

##### Richness


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
Groups     2 0.017675 0.0088373 2.3825    999  0.099 .
Residuals 23 0.085312 0.0037092                       
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
              Control Hot_control Treatment
Control                 0.0060000     0.661
Hot_control 0.0068795                 0.212
Treatment   0.6248469   0.2084296          
```


|         | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|---------:|---------:|--------:|------:|
|species  |  1| 0.6340254| 0.0768024| 2.065607|  0.004|
|type     |  1| 0.5615418| 0.0680222| 1.829461|  0.010|
|Residual | 23| 7.0597099| 0.8551754|       NA|     NA|
|Total    | 25| 8.2552771| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
1     Control vs Treatment  1 0.5615418 1.729004 0.1033537   0.016      0.048   .
2   Control vs Hot_control  1 0.8438429 2.793772 0.1486541   0.002      0.006   *
3 Treatment vs Hot_control  1 0.3734921 1.268929 0.0779971   0.098      0.294    
```

##### Neutral


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.011001 0.0055005 0.6303    999  0.564
Residuals 23 0.200714 0.0087267                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.21500     0.957
Hot_control 0.21166                 0.467
Treatment   0.95468     0.43604          
```


|         | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|---------:|---------:|--------:|------:|
|species  |  1| 0.7907904| 0.1076445| 3.056657|  0.001|
|type     |  1| 0.6051778| 0.0823784| 2.339205|  0.009|
|Residual | 23| 5.9503501| 0.8099772|       NA|     NA|
|Total    | 25| 7.3463184| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
1     Control vs Treatment  1 0.6051778 2.250849 0.13047758   0.009      0.027   .
2   Control vs Hot_control  1 1.0528902 4.143637 0.20570451   0.001      0.003   *
3 Treatment vs Hot_control  1 0.4150076 1.637268 0.09840968   0.046      0.138    
```

##### Phylogenetic


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.00440 0.0021994 0.1369    999   0.91
Residuals 23 0.36941 0.0160614                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.92400     0.694
Hot_control 0.91505                 0.803
Treatment   0.63312     0.73046          
```


|         | Df|  SumOfSqs|        R2|         F| Pr(>F)|
|:--------|--:|---------:|---------:|---------:|------:|
|species  |  1| 0.0560850| 0.0531376| 1.3149967|  0.280|
|type     |  1| 0.0184254| 0.0174571| 0.4320099|  0.791|
|Residual | 23| 0.9809570| 0.9294053|        NA|     NA|
|Total    | 25| 1.0554673| 1.0000000|        NA|     NA|


```{.script-output}
                     pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
1     Control vs Treatment  1 0.01842535 0.4144162 0.02688498   0.775      1.000    
2   Control vs Hot_control  1 0.05987967 1.7387847 0.09802164   0.117      0.351    
3 Treatment vs Hot_control  1 0.03212966 0.6477782 0.04139746   0.692      1.000    
```

##### Functional


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq   Mean Sq     F N.Perm Pr(>F)
Groups     2 0.00400 0.0020014 0.145    999  0.864
Residuals 23 0.31753 0.0138057                    

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.59900     0.750
Hot_control 0.59817                 0.849
Treatment   0.75141     0.83718          
```


|         | Df|  SumOfSqs|        R2|         F| Pr(>F)|
|:--------|--:|---------:|---------:|---------:|------:|
|species  |  1| 0.0024979| 0.0033024| 0.0900845|  0.638|
|type     |  1| 0.1161466| 0.1535542| 4.1887855|  0.063|
|Residual | 23| 0.6377435| 0.8431434|        NA|     NA|
|Total    | 25| 0.7563879| 1.0000000|        NA|     NA|


```{.script-output}
                     pairs Df  SumsOfSqs  F.Model         R2 p.value p.adjusted sig
1     Control vs Treatment  1 0.11614656 4.724791 0.23953568   0.069      0.207    
2   Control vs Hot_control  1 0.05000930 1.704826 0.09629160   0.224      0.672    
3 Treatment vs Hot_control  1 0.01235859 0.423812 0.02747777   0.503      1.000    
```



```{.r .script-source}
beta_richness_nmds_post1 <- beta_div_richness_post1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post1_nmds, by = join_by(sample == Tube_code))

beta_neutral_nmds_post1 <- beta_div_neutral_post1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post1_nmds, by = join_by(sample == Tube_code))

beta_phylogenetic_nmds_post1 <- beta_div_phylo_post1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post1_nmds, by = join_by(sample == Tube_code))

beta_functional_nmds_post1 <- beta_div_func_post1$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post1_nmds, by = join_by(sample == Tube_code))
```


```{.r .script-source}
p0<-beta_richness_nmds_post1 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
        scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y = "NMDS2", x="NMDS1 \n Richness beta diversity") +
				theme_classic() +
				theme(legend.position="none")

p1<-beta_neutral_nmds_post1 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
        scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
				theme_classic() +
				theme(legend.position="none")
  
p2<-beta_phylogenetic_nmds_post1 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
				scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y= element_blank (), x="NMDS1 \n Phylogenetic beta diversity") +
				theme_classic() +
				theme(legend.position="none")

p3<-beta_functional_nmds_post1 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
				scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y= element_blank (), x="NMDS1 \n Functional beta diversity") +
				theme_classic()+
				theme(legend.position="none")
```


```{.r .script-source}
ggarrange(p0, p1, p2, p3, ncol=2, nrow=2, common.legend = TRUE, legend="right")
```

<img src="_main_files/figure-html/beta_neutral_nmds_plot_post1_final-1.png" width="768" />


#### After 2 weeks -->Post-FMT2


```{.r .script-source}
post2 <- meta %>%
  filter(time_point == "6_Post-FMT2")

post2.counts <- temp_genome_counts[,which(colnames(temp_genome_counts) %in% rownames(post2))]
identical(sort(colnames(post2.counts)),sort(as.character(rownames(post2))))

post2_nmds <- sample_metadata %>%
  filter(time_point == "6_Post-FMT2")
```

#### Number of samples used


```{.script-output}
[1] 27
```


```{.r .script-source}
beta_div_richness_post2<-hillpair(data=post2.counts, q=0)
beta_div_neutral_post2<-hillpair(data=post2.counts, q=1)
beta_div_phylo_post2<-hillpair(data=post2.counts, q=1, tree=genome_tree)
beta_div_func_post2<-hillpair(data=post2.counts, q=1, dist=dist)

#Arrange of metadata dataframe
post2_arrange<-post2[labels(beta_div_neutral_post2$S),]
```

##### Richness


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.002011 0.0010056 0.1982    999   0.83
Residuals 24 0.121775 0.0050740                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.70300     0.786
Hot_control 0.67789                 0.624
Treatment   0.79246     0.59820          
```


|         | Df| SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|--------:|---------:|--------:|------:|
|type     |  2| 1.504341| 0.1967776| 2.939822|  0.001|
|Residual | 24| 6.140538| 0.8032224|       NA|     NA|
|Total    | 26| 7.644879| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
1     Treatment vs Control  1 0.6463814 2.560441 0.1379515   0.001      0.003   *
2 Treatment vs Hot_control  1 0.4796256 1.916520 0.1069694   0.001      0.003   *
3   Control vs Hot_control  1 1.1305044 4.268317 0.2105906   0.001      0.003   *
```

##### Neutral


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.008262 0.0041311 0.8024    999  0.468
Residuals 24 0.123559 0.0051483                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.46700     0.694
Hot_control 0.44675                 0.241
Treatment   0.65989     0.25095          
```


|         | Df| SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|--------:|---------:|--------:|------:|
|type     |  2| 1.923807| 0.2603795| 4.224537|  0.001|
|Residual | 24| 5.464666| 0.7396205|       NA|     NA|
|Total    | 26| 7.388473| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
1     Treatment vs Control  1 1.0227481 4.648335 0.2251191   0.001      0.003   *
2 Treatment vs Hot_control  1 0.5010202 2.206532 0.1211945   0.001      0.003   *
3   Control vs Hot_control  1 1.3619424 5.771031 0.2650785   0.001      0.003   *
```

##### Phylogenetic


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.000407 0.0002034 0.0487    999  0.959
Residuals 24 0.100305 0.0041794                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.92800     0.874
Hot_control 0.93765                 0.769
Treatment   0.83933     0.76015          
```



|         | Df|  SumOfSqs|        R2|        F| Pr(>F)|
|:--------|--:|---------:|---------:|--------:|------:|
|type     |  2| 0.1594363| 0.2042241| 3.079623|  0.001|
|Residual | 24| 0.6212564| 0.7957759|       NA|     NA|
|Total    | 26| 0.7806927| 1.0000000|       NA|     NA|


```{.script-output}
                     pairs Df  SumsOfSqs  F.Model        R2 p.value p.adjusted sig
1     Treatment vs Control  1 0.05927454 2.382025 0.1295845   0.021      0.063    
2 Treatment vs Hot_control  1 0.06906280 2.722460 0.1454115   0.005      0.015   .
3   Control vs Hot_control  1 0.11081709 4.043656 0.2017424   0.001      0.003   *
```

##### Functional


```{.script-output}

Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
Groups     2 0.01126 0.0056302 0.2861    999  0.773
Residuals 24 0.47233 0.0196806                     

Pairwise comparisons:
(Observed p-value below diagonal, permuted p-value above diagonal)
            Control Hot_control Treatment
Control                 0.54800     0.634
Hot_control 0.48255                 0.790
Treatment   0.60116     0.75643          
```


|         | Df|   SumOfSqs|         R2|          F| Pr(>F)|
|:--------|--:|----------:|----------:|----------:|------:|
|type     |  2| -0.0038724| -0.0056213| -0.0670788|  0.936|
|Residual | 24|  0.6927468|  1.0056213|         NA|     NA|
|Total    | 26|  0.6888744|  1.0000000|         NA|     NA|


```{.script-output}
                     pairs Df    SumsOfSqs     F.Model           R2 p.value p.adjusted sig
1     Treatment vs Control  1 -0.008527330 -0.46290555 -0.029793572   0.853          1    
2 Treatment vs Hot_control  1 -0.001648721 -0.04717131 -0.002956924   0.910          1    
3   Control vs Hot_control  1  0.004367477  0.13147026  0.008149924   0.690          1    
```



```{.r .script-source}
beta_richness_nmds_post2 <- beta_div_richness_post2$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post2_nmds, by = join_by(sample == Tube_code))

beta_neutral_nmds_post2 <- beta_div_neutral_post2$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post2_nmds, by = join_by(sample == Tube_code))

beta_phylogenetic_nmds_post2 <- beta_div_phylo_post2$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post2_nmds, by = join_by(sample == Tube_code))

beta_functional_nmds_post2 <- beta_div_func_post2$S %>%
				metaMDS(.,trymax = 500, k=2, verbosity=FALSE) %>%
				scores() %>%
				as_tibble(., rownames = "sample") %>%
				left_join(post2_nmds, by = join_by(sample == Tube_code))
```



```{.r .script-source}
p0<-beta_richness_nmds_post2 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
        scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y = "NMDS2", x="NMDS1 \n Richness beta diversity") +
				theme_classic() +
				theme(legend.position="none")

p1<-beta_neutral_nmds_post2 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
        scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y = "NMDS2", x="NMDS1 \n Neutral beta diversity") +
				theme_classic() +
				theme(legend.position="none")
  
p2<-beta_phylogenetic_nmds_post2 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
				scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y= element_blank (), x="NMDS1 \n Phylogenetic beta diversity") +
				theme_classic() +
				theme(legend.position="none")

p3<-beta_functional_nmds_post2 %>%
			group_by(type) %>%
			mutate(x_cen = mean(NMDS1, na.rm = TRUE)) %>%
			mutate(y_cen = mean(NMDS2, na.rm = TRUE)) %>%
			ungroup() %>%
			ggplot(., aes(x=NMDS1,y=NMDS2, color=type)) +
				scale_color_manual(name="Type",
                       breaks=c("Control", "Hot_control", "Treatment"),
                       labels=c("Cold-Cold", "Hot-Hot", "Cold-Hot"),
                       values=c("#4477AA","#d57d2c","#76b183")) +
				geom_point(size=2) +
				geom_segment(aes(x=x_cen, y=y_cen, xend=NMDS1, yend=NMDS2), alpha=0.2) +
        labs(y= element_blank (), x="NMDS1 \n Functional beta diversity") +
				theme_classic()+
				theme(legend.position="none")
```


```{.r .script-source}
ggarrange(p0, p1, p2, p3, ncol=2, nrow=2, common.legend = TRUE, legend="right")
```

<img src="_main_files/figure-html/beta_neutral_nmds_plot_post2_final-1.png" width="768" />

<!--chapter:end:05_community_composition_comparisons.Rmd-->

# Functional differences


```{.r .script-source}
load("data/data.Rdata")
```

## Data preparation


```{.r .script-source}
# Aggregate bundle-level GIFTs into the compound level
GIFTs_elements <- to.elements(genome_gifts, GIFT_db)
GIFTs_elements_filtered <- GIFTs_elements[rownames(GIFTs_elements) %in% genome_counts$genome, ]
GIFTs_elements_filtered <- as.data.frame(GIFTs_elements_filtered) %>%
  select_if(~ !is.numeric(.) || sum(.) != 0)

elements <- GIFTs_elements_filtered %>%
  as.data.frame()

# Aggregate element-level GIFTs into the function level
GIFTs_functions <- to.functions(GIFTs_elements_filtered, GIFT_db)
functions <- GIFTs_functions %>%
  as.data.frame()

# Aggregate function-level GIFTs into overall Biosynthesis, Degradation and Structural GIFTs
GIFTs_domains <- to.domains(GIFTs_functions, GIFT_db)
domains <- GIFTs_domains %>%
  as.data.frame()

# Get community-weighed average GIFTs per sample
GIFTs_elements_community <- to.community(GIFTs_elements_filtered, genome_counts_filt %>% column_to_rownames(., "genome") %>% tss(), GIFT_db)
GIFTs_functions_community <- to.community(GIFTs_functions, genome_counts_filt %>% column_to_rownames(., "genome") %>% tss(), GIFT_db)
GIFTs_domains_community <- to.community(GIFTs_domains, genome_counts_filt %>% column_to_rownames(., "genome") %>% tss(), GIFT_db)

uniqueGIFT_db<- unique(GIFT_db[c(2,4,5,6)]) %>% unite("Function",Function:Element, sep= "_", remove=FALSE)
```

## Genomes GIFT profiles


```{.r .script-source}
GIFTs_elements %>%
  as_tibble(., rownames = "MAG") %>%
  reshape2::melt() %>%
  rename(Code_element = variable, GIFT = value) %>%
  inner_join(GIFT_db,by="Code_element") %>%
  ggplot(., aes(x=Code_element, y=MAG, fill=GIFT, group=Code_function))+
    geom_tile()+
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE))+
    scale_fill_gradientn(colours=rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da")))+
    facet_grid(. ~ Code_function, scales = "free", space = "free")+
    theme_grey(base_size=8)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text.x = element_text(angle = 90))
```

<img src="_main_files/figure-html/gift_function_heatmap-1.png" width="960" />

## Function level


```{.r .script-source}
GIFTs_functions_community %>%
    as.data.frame() %>%
    rownames_to_column(var="sample") %>%
    filter(sample!="AD69") %>%
    pivot_longer(!sample,names_to="trait",values_to="gift") %>%
    left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
    ggplot(aes(x=trait,y=time_point,fill=gift)) +
        geom_tile(colour="white", size=0.2)+
        scale_fill_gradientn(colours=rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da")))+
        facet_grid(type ~ ., scales="free",space="free")
```

<img src="_main_files/figure-html/gift_function_heatmap_1-1.png" width="960" />

## Element level


```{.r .script-source}
GIFTs_elements_community_merged<-GIFTs_elements_community %>%
    as.data.frame() %>%
    rownames_to_column(var="sample") %>%
    filter(sample!="AD69") %>%
    pivot_longer(!sample,names_to="trait",values_to="gift") %>%
    left_join(sample_metadata, by = join_by(sample == Tube_code))%>%
    mutate(functionid = substr(trait, 1, 3)) %>%
    mutate(trait = case_when(
      trait %in% GIFT_db$Code_element ~ GIFT_db$Element[match(trait, GIFT_db$Code_element)],
      TRUE ~ trait
    )) %>%
    mutate(functionid = case_when(
      functionid %in% GIFT_db$Code_function ~ GIFT_db$Function[match(functionid, GIFT_db$Code_function)],
      TRUE ~ functionid
    )) %>%
    mutate(trait=factor(trait,levels=unique(GIFT_db$Element))) %>%
    mutate(functionid=factor(functionid,levels=unique(GIFT_db$Function)))

# Create an interaction variable for time_point and sample
GIFTs_elements_community_merged$interaction_var <- interaction(GIFTs_elements_community_merged$sample, GIFTs_elements_community_merged$time_point)
  
ggplot(GIFTs_elements_community_merged,aes(x=interaction_var,y=trait,fill=gift)) +
        geom_tile(colour="white", linewidth=0.2)+
        scale_fill_gradientn(colours=rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da")))+
        facet_grid(functionid ~ type, scales="free",space="free") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),
              strip.text.y = element_text(angle = 0)) + 
        labs(y="Traits",x="Time_point",fill="GIFT")+
  scale_x_discrete(labels = function(x) gsub(".*\\.", "", x))
```

<img src="_main_files/figure-html/gift_element_heatmap-1.png" width="960" />

## Comparison of samples from the 0 Time_point (0_Wild)

### GIFTs Functional community


```{.r .script-source}
GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="0_Wild") %>%
  group_by(species) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 2 × 3
  species             MCI     sd
  <chr>             <dbl>  <dbl>
1 Podarcis_liolepis 0.327 0.0244
2 Podarcis_muralis  0.346 0.0194
```

#### GIFT test visualisation


```{.r .script-source}
GIFTs_functions_community %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="0_Wild") %>%
  select(c(1:21, 27)) %>%
  pivot_longer(-c(sample,type),names_to = "trait", values_to = "value") %>%
  mutate(trait = case_when(
      trait %in% GIFT_db$Code_function ~ GIFT_db$Function[match(trait, GIFT_db$Code_function)],
      TRUE ~ trait
    )) %>%
  mutate(trait=factor(trait,levels=unique(GIFT_db$Function))) %>%
  ggplot(aes(x=value, y=type, group=type, fill=type, color=type)) +
    geom_boxplot() +
    scale_color_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c","#76b183")) +
      scale_fill_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
    facet_grid(trait ~ ., space="free", scales="free") +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.y = element_text(angle = 0)) + 
        labs(y="Traits",x="Metabolic capacity index")
```

<img src="_main_files/figure-html/gift_test_function_plot_wild-1.png" width="960" />

### GIFTs Domain community


```{.r .script-source}
GIFTs_domains_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="0_Wild") %>%
  group_by(species) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 2 × 3
  species             MCI     sd
  <chr>             <dbl>  <dbl>
1 Podarcis_liolepis 0.375 0.0315
2 Podarcis_muralis  0.390 0.0208
```

### GIFTs Elements community


```{.r .script-source}
GIFTs_elements_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="0_Wild") %>%
  group_by(species) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 2 × 3
  species             MCI     sd
  <chr>             <dbl>  <dbl>
1 Podarcis_liolepis 0.313 0.0329
2 Podarcis_muralis  0.345 0.0233
```


```{.r .script-source}
sample_metadata_wild <- sample_metadata%>% 
  filter(time_point == "0_Wild")

element_gift_wild <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(., sample_metadata_wild[c(1,3)], by="Tube_code")
```


```{.r .script-source}
# Find numeric columns
numeric_cols <- sapply(element_gift_wild, is.numeric)

# Calculate column sums for numeric columns only
col_sums_numeric <- colSums(element_gift_wild[, numeric_cols])

# Identify numeric columns with sums not equal to zero
nonzero_numeric_cols <- names(col_sums_numeric)[col_sums_numeric != 0]

# Remove numeric columns with sums not equal to zero
filtered_data <- element_gift_wild[, !numeric_cols | colnames(element_gift_wild) %in% nonzero_numeric_cols]
```


```{.r .script-source}
significant_elements_wild <- filtered_data %>%
  pivot_longer(-c(Tube_code,species), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ species, exact=FALSE)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t <- element_gift_wild  %>% 
  dplyr::select(-c(species))  %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt <- subset(element_gift_t, trait %in% significant_elements_wild$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_wild[c(1,3)], by = join_by(Tube_code == Tube_code))

element_gift_filt %>%
  dplyr::select(-Tube_code)%>%
  group_by(species)  %>%
  summarise(across(everything(), mean))%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))

element_gift_names <- element_gift_filt%>%
  dplyr::select(-species)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  dplyr::select(-Elements)%>%
  dplyr::select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_wild[c(1,4)], by = join_by(Tube_code == Tube_code))
```


```{.r .script-source}
colNames <- names(element_gift_names)[2:30] #always check names(element_gift_names) first to know where your traits finish
for(i in colNames){
  plt <- ggplot(element_gift_names, aes(x=Population, y=.data[[i]], color = Population)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}
```

<img src="_main_files/figure-html/gitfs_plot_1-1.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-2.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-3.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-4.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-5.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-6.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-7.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-8.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-9.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-10.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-11.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-12.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-13.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-14.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-15.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-16.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-17.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-18.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-19.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-20.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-21.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-22.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-23.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-24.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-25.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-26.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-27.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-28.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1-29.png" width="768" />

## Comparison of samples from the 1st Time_point (1_Acclimation)

### GIFTs Functional community


```{.r .script-source}
GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="1_Acclimation") %>%
  group_by(species) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 2 × 3
  species             MCI     sd
  <chr>             <dbl>  <dbl>
1 Podarcis_liolepis 0.348 0.0158
2 Podarcis_muralis  0.331 0.0321
```

#### GIFT test visualisation


```{.r .script-source}
GIFTs_functions_community %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="1_Acclimation") %>%
  select(c(1:21, 27)) %>%
  pivot_longer(-c(sample,type),names_to = "trait", values_to = "value") %>%
  mutate(trait = case_when(
      trait %in% GIFT_db$Code_function ~ GIFT_db$Function[match(trait, GIFT_db$Code_function)],
      TRUE ~ trait
    )) %>%
  mutate(trait=factor(trait,levels=unique(GIFT_db$Function))) %>%
  ggplot(aes(x=value, y=type, group=type, fill=type, color=type)) +
    geom_boxplot() +
    scale_color_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c","#76b183")) +
      scale_fill_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
    facet_grid(trait ~ ., space="free", scales="free") +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.y = element_text(angle = 0)) + 
        labs(y="Traits",x="Metabolic capacity index")
```

<img src="_main_files/figure-html/gift_test_function_plot_accli-1.png" width="960" />

### GIFTs Domain community


```{.r .script-source}
GIFTs_domains_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="1_Acclimation") %>%
  group_by(species) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 2 × 3
  species             MCI     sd
  <chr>             <dbl>  <dbl>
1 Podarcis_liolepis 0.395 0.0211
2 Podarcis_muralis  0.370 0.0307
```

### GIFTs Elements community


```{.r .script-source}
GIFTs_elements_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="1_Acclimation") %>%
  group_by(species) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 2 × 3
  species             MCI     sd
  <chr>             <dbl>  <dbl>
1 Podarcis_liolepis 0.350 0.0225
2 Podarcis_muralis  0.332 0.0316
```


```{.r .script-source}
sample_metadata_accli <- sample_metadata%>% 
  filter(time_point == "1_Acclimation")

element_gift_accli <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(., sample_metadata_accli[c(1,3)], by="Tube_code")
```


```{.r .script-source}
# Find numeric columns
numeric_cols <- sapply(element_gift_accli, is.numeric)

# Calculate column sums for numeric columns only
col_sums_numeric <- colSums(element_gift_accli[, numeric_cols])

# Identify numeric columns with sums not equal to zero
nonzero_numeric_cols <- names(col_sums_numeric)[col_sums_numeric != 0]

# Remove numeric columns with sums not equal to zero
filtered_data <- element_gift_accli[, !numeric_cols | colnames(element_gift_accli) %in% nonzero_numeric_cols]
```


```{.r .script-source}
significant_elements_accli <- filtered_data %>%
  pivot_longer(-c(Tube_code,species), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ species, exact=FALSE)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_adjust < 0.05)  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t <- element_gift_accli  %>% 
  dplyr::select(-c(species))  %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt <- subset(element_gift_t, trait %in% significant_elements_accli$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_accli[c(1,3)], by = join_by(Tube_code == Tube_code))

element_gift_filt %>%
  dplyr::select(-Tube_code)%>%
  group_by(species)  %>%
  summarise(across(everything(), mean))%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))

element_gift_names <- element_gift_filt%>%
  dplyr::select(-species)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  dplyr::select(-Elements)%>%
  dplyr::select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_accli[c(1,4)], by = join_by(Tube_code == Tube_code))
```


```{.r .script-source}
colNames <- names(element_gift_names)[2:14] #always check names(element_gift_names) first to know where your traits finish
for(i in colNames){
  plt <- ggplot(element_gift_names, aes(x=Population, y=.data[[i]], color = Population)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}
```

<img src="_main_files/figure-html/gitfs_plot_1_1-1.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-2.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-3.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-4.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-5.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-6.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-7.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-8.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-9.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-10.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-11.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-12.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1-13.png" width="768" />



## Comparison of samples from the 5th Time_point (5_Post-FMT1)

### GIFTs Functional community


```{.r .script-source}
GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="5_Post-FMT1") %>%
  group_by(type) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 3 × 3
  type          MCI     sd
  <chr>       <dbl>  <dbl>
1 Control     0.373 0.0247
2 Hot_control 0.372 0.0367
3 Treatment   0.353 0.0186
```

#### GIFT test visualisation


```{.r .script-source}
GIFTs_functions_community %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="5_Post-FMT1") %>%
  select(c(1:21, 27)) %>%
  pivot_longer(-c(sample,type),names_to = "trait", values_to = "value") %>%
  mutate(trait = case_when(
      trait %in% GIFT_db$Code_function ~ GIFT_db$Function[match(trait, GIFT_db$Code_function)],
      TRUE ~ trait
    )) %>%
  mutate(trait=factor(trait,levels=unique(GIFT_db$Function))) %>%
  ggplot(aes(x=value, y=type, group=type, fill=type, color=type)) +
    geom_boxplot() +
    scale_color_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c","#76b183")) +
      scale_fill_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
    facet_grid(trait ~ ., space="free", scales="free") +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.y = element_text(angle = 0)) + 
        labs(y="Traits",x="Metabolic capacity index")
```

<img src="_main_files/figure-html/gift_test_function_plot_tm5-1.png" width="960" />

### GIFTs Domain community


```{.r .script-source}
GIFTs_domains_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="5_Post-FMT1") %>%
  group_by(type) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 3 × 3
  type          MCI     sd
  <chr>       <dbl>  <dbl>
1 Control     0.412 0.0213
2 Hot_control 0.415 0.0437
3 Treatment   0.391 0.0263
```

### GIFTs Elements community


```{.r .script-source}
GIFTs_elements_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="5_Post-FMT1") %>%
  group_by(type) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 3 × 3
  type          MCI     sd
  <chr>       <dbl>  <dbl>
1 Control     0.380 0.0280
2 Hot_control 0.379 0.0372
3 Treatment   0.359 0.0214
```


```{.r .script-source}
sample_metadata_tm5 <- sample_metadata%>% 
  filter(time_point == "5_Post-FMT1")%>% 
  filter(type != "Hot_control")

element_gift_tm5 <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(sample_metadata_tm5 %>% select(1, 7), by = "Tube_code")
```


```{.r .script-source}
# Find numeric columns
numeric_cols <- sapply(element_gift_tm5, is.numeric)

# Calculate column sums for numeric columns only
col_sums_numeric <- colSums(element_gift_tm5[, numeric_cols])

# Identify numeric columns with sums not equal to zero
nonzero_numeric_cols <- names(col_sums_numeric)[col_sums_numeric != 0]

# Remove numeric columns with sums not equal to zero
filtered_data <- element_gift_tm5[, !numeric_cols | colnames(element_gift_tm5) %in% nonzero_numeric_cols]
```


```{.r .script-source}
significant_elements_tm5 <- filtered_data %>%
  pivot_longer(-c(Tube_code,type), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ type, exact=FALSE)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_value < 0.05)  %>% #take into account that p_value is used and not p_adjust
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t <- element_gift_tm5  %>% 
  dplyr::select(-c(type))  %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt <- subset(element_gift_t, trait %in% significant_elements_tm5$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_tm5[c(1,7)], by = join_by(Tube_code == Tube_code))

element_gift_filt %>%
  dplyr::select(-Tube_code)%>%
  group_by(type)  %>%
  summarise(across(everything(), mean))%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))

element_gift_names <- element_gift_filt%>%
  dplyr::select(-type)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  dplyr::select(-Elements)%>%
  dplyr::select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_tm5[c(1,7)], by = join_by(Tube_code == Tube_code))
```


```{.r .script-source}
colNames <- names(element_gift_names)[2:14] #always check names(element_gift_names) first to now where your traits finish
for(i in colNames){
  plt <- ggplot(element_gift_names, aes(x=type, y=.data[[i]], color = type)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}
```

<img src="_main_files/figure-html/gitfs_plot-1.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-2.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-3.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-4.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-5.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-6.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-7.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-8.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-9.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-10.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-11.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-12.png" width="768" /><img src="_main_files/figure-html/gitfs_plot-13.png" width="768" />

## Comparison of samples from the 6th Time_point (6_Post-FMT2)

### GIFTs Functional community


```{.r .script-source}
GIFTs_functions_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="6_Post-FMT2") %>%
  group_by(type) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 3 × 3
  type          MCI     sd
  <chr>       <dbl>  <dbl>
1 Control     0.352 0.0223
2 Hot_control 0.350 0.0293
3 Treatment   0.346 0.0255
```

#### GIFT test visualisation


```{.r .script-source}
GIFTs_functions_community %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="6_Post-FMT2") %>%
  select(c(1:21, 27)) %>%
  pivot_longer(-c(sample,type),names_to = "trait", values_to = "value") %>%
  mutate(trait = case_when(
      trait %in% GIFT_db$Code_function ~ GIFT_db$Function[match(trait, GIFT_db$Code_function)],
      TRUE ~ trait
    )) %>%
  mutate(trait=factor(trait,levels=unique(GIFT_db$Function))) %>%
  ggplot(aes(x=value, y=type, group=type, fill=type, color=type)) +
    geom_boxplot() +
    scale_color_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA","#d57d2c","#76b183")) +
      scale_fill_manual(name="type",
          breaks=c("Control","Hot_control", "Treatment"),
          labels=c("Cold-Cold","Hot-Hot", "Cold-Hot"),
          values=c("#4477AA50","#d57d2c50","#76b18350")) +
    facet_grid(trait ~ ., space="free", scales="free") +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              strip.text.y = element_text(angle = 0)) + 
        labs(y="Traits",x="Metabolic capacity index")
```

<img src="_main_files/figure-html/gift_test_function_plot-1.png" width="960" />

### GIFTs Domain community


```{.r .script-source}
GIFTs_domains_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="6_Post-FMT2") %>%
  group_by(type) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 3 × 3
  type          MCI     sd
  <chr>       <dbl>  <dbl>
1 Control     0.399 0.0171
2 Hot_control 0.388 0.0321
3 Treatment   0.392 0.0240
```

### GIFTs Elements community


```{.r .script-source}
GIFTs_elements_community %>%
  rowMeans() %>%
  as_tibble(., rownames = "sample") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>%
  filter(time_point=="6_Post-FMT2") %>%
  group_by(type) %>%
  summarise(MCI = mean(value), sd = sd(value))
```

```{.script-output}
# A tibble: 3 × 3
  type          MCI     sd
  <chr>       <dbl>  <dbl>
1 Control     0.357 0.0215
2 Hot_control 0.347 0.0302
3 Treatment   0.350 0.0293
```


```{.r .script-source}
sample_metadata_TM6 <- sample_metadata%>% 
  filter(time_point == "6_Post-FMT2")%>% 
  filter(type != "Hot_control")

element_gift_TM6 <- GIFTs_elements_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(sample_metadata_TM6 %>% select(1, 7), by = "Tube_code")
```


```{.r .script-source}
# Find numeric columns
numeric_cols <- sapply(element_gift_TM6, is.numeric)

# Calculate column sums for numeric columns only
col_sums_numeric <- colSums(element_gift_TM6[, numeric_cols])

# Identify numeric columns with sums not equal to zero
nonzero_numeric_cols <- names(col_sums_numeric)[col_sums_numeric != 0]

# Remove numeric columns with sums not equal to zero
filtered_data <- element_gift_TM6[, !numeric_cols | colnames(element_gift_TM6) %in% nonzero_numeric_cols]
```


```{.r .script-source}
significant_elements_TM6 <- filtered_data %>%
  pivot_longer(-c(Tube_code,type), names_to = "trait", values_to = "value") %>%
  group_by(trait) %>%
  summarise(p_value = wilcox.test(value ~ type, exact=FALSE)$p.value) %>%
  mutate(p_adjust=p.adjust(p_value, method="BH")) %>%
  filter(p_value < 0.05)  %>% #take into account that p_value is used and not p_adjust
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(trait == Code_element))

element_gift_t <- element_gift_TM6  %>% 
  dplyr::select(-c(type))  %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "trait")

element_gift_filt <- subset(element_gift_t, trait %in% significant_elements_TM6$trait) %>% 
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_TM6[c(1,7)], by = join_by(Tube_code == Tube_code))

element_gift_filt %>%
  dplyr::select(-Tube_code)%>%
  group_by(type)  %>%
  summarise(across(everything(), mean))%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))

element_gift_names <- element_gift_filt%>%
  dplyr::select(-type)%>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Elements")  %>%
  left_join(.,uniqueGIFT_db[c(1,3)],by = join_by(Elements == Code_element))%>%
  dplyr::select(-Elements)%>%
  dplyr::select(Function, everything())%>%
  t()%>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate_if(is.character, as.numeric)  %>%
  rownames_to_column(., "Tube_code")%>% 
  left_join(., sample_metadata_TM6[c(1,7)], by = join_by(Tube_code == Tube_code))
```


```{.r .script-source}
colNames <- names(element_gift_names)[2:20] #always check names(element_gift_names) first to now where your traits finish
for(i in colNames){
  plt <- ggplot(element_gift_names, aes(x=type, y=.data[[i]], color = type)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA, width = 0.3, show.legend = FALSE) +
    geom_jitter(width = 0.1, show.legend = TRUE) +
    theme_minimal() +
    theme(
      axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank())
  print(plt)
}
```

<img src="_main_files/figure-html/gitfs_plot_1_1_1-1.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-2.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-3.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-4.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-5.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-6.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-7.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-8.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-9.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-10.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-11.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-12.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-13.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-14.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-15.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-16.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-17.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-18.png" width="768" /><img src="_main_files/figure-html/gitfs_plot_1_1_1-19.png" width="768" />

## Domain level

### Comparison of samples from the 0 Time_point (0_Wild)


```{.r .script-source}
#Merge the functional domains with the metadata
merge_gift_wild<- GIFTs_domains_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(., sample_metadata_wild, by="Tube_code")
```



```{.r .script-source}
#Biosynthesis
p1 <-merge_gift_wild %>%
  ggplot(aes(x=species,y=Biosynthesis,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Degradation
p2 <-merge_gift_wild %>%
  ggplot(aes(x=species,y=Degradation,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.45, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Structure
p3 <-merge_gift_wild %>%
  ggplot(aes(x=species,y=Structure,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 3, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Overall
p4 <-merge_gift_wild %>%
  ggplot(aes(x=species,y=Overall,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")
```

<img src="_main_files/figure-html/final_plot_1_2-1.png" width="960" />

### Comparison of samples from the 1st Time_point (1_Acclimation)


```{.r .script-source}
#Merge the functional domains with the metadata
merge_gift_accli<- GIFTs_domains_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(., sample_metadata_accli, by="Tube_code")
```



```{.r .script-source}
#Biosynthesis
p1 <-merge_gift_accli %>%
  ggplot(aes(x=species,y=Biosynthesis,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Degradation
p2 <-merge_gift_accli %>%
  ggplot(aes(x=species,y=Degradation,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.45, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Structure
p3 <-merge_gift_accli %>%
  ggplot(aes(x=species,y=Structure,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 3, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Overall
p4 <-merge_gift_accli %>%
  ggplot(aes(x=species,y=Overall,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")
```

<img src="_main_files/figure-html/final_plot_1_1-1.png" width="960" />

### Comparison of samples from the 5th Time_point (5_Post-FMT1)


```{.r .script-source}
#Merge the functional domains with the metadata
merge_gift_tm5<- GIFTs_domains_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(., sample_metadata_tm5, by="Tube_code")
```



```{.r .script-source}
#Biosynthesis
p1 <-merge_gift_tm5 %>%
  ggplot(aes(x=species,y=Biosynthesis,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Degradation
p2 <-merge_gift_tm5 %>%
  ggplot(aes(x=species,y=Degradation,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.45, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Structure
p3 <-merge_gift_tm5 %>%
  ggplot(aes(x=species,y=Structure,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 3, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")

#Overall
p4 <-merge_gift_tm5 %>%
  ggplot(aes(x=species,y=Overall,color=species,fill=species))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Species")
```


```
Warning: Computation failed in `stat_compare_means()`.
Computation failed in `stat_compare_means()`.
Computation failed in `stat_compare_means()`.
Computation failed in `stat_compare_means()`.
Caused by error:
! argument "x" is missing, with no default
```

<img src="_main_files/figure-html/final_plot_1-1.png" width="960" />

### Comparison of samples from the 6th Time_point (6_Post-FMT2)


```{.r .script-source}
#Merge the functional domains with the metadata
merge_gift_TM6 <- GIFTs_domains_community %>% 
  as.data.frame() %>% 
  rownames_to_column(., "Tube_code") %>% 
  inner_join(., sample_metadata_TM6, by="Tube_code")
```


```{.r .script-source}
#Biosynthesis
p1 <-merge_gift_TM6 %>%
  ggplot(aes(x=type,y=Biosynthesis,color=type,fill=type))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Type")

#Degradation
p2 <-merge_gift_TM6 %>%
  ggplot(aes(x=type,y=Degradation,color=type,fill=type))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "type")

#Structure
p3 <-merge_gift_TM6 %>%
  ggplot(aes(x=type,y=Structure,color=type,fill=type))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 3, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "type")

#Overall
p4 <-merge_gift_TM6 %>%
  ggplot(aes(x=type,y=Overall,color=type,fill=type))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "type")
```

<img src="_main_files/figure-html/final_plot-1.png" width="960" />



<!--chapter:end:06_functional_differences.Rmd-->

