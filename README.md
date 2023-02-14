
# Overview of our DArTseq analyses

This git repository explains our analytical workflow of DArTseq data to generate population genetic outputs. Here is a general overview of how to set up files, the type of analyses and outputs we obtain and how we interpret the results.

The analyses are split into individual R files so that more in depth explaining can be done for each analysis. The files have been numbered to make it easier to navigate through. Typically we combine the scripts into a single file to keep track of the analyses that are run on a specific dataset.

## DArTseq data

### How to use?

Our analyses are split into R files which focus on one analysis at a time. This makes particular plots or calculations easy to find in the git repository, however to use them they will need to be combined into a single document. I recommend that you download the first file (`1 data import and qc.Rmd`) and then add the subsequent code that you want to use as chunks. The R files are numbered in the order that I would usually run them, however this can change depending on your data and goals. Below you will find a summary of what these files are.

### Directories

You will have to set up your directories in a specific way for this to work, as per the original ReCER method. For example, if I'm working with the species *Zieria obcordata* I would have the directories set up as:

    --ZierObco
      |-- Analyses.Rmd
      `--ZierObco
        |--dart_raw
        |   `-Report_DZ22-7321_SNP_mapping_2.csv
        |--meta
        |   `-ZierObco_DZ22-7321_meta.csv
        |--popgen
        `--outputs

Where `Report_DZ22-7321_SNP_mapping_2.csv` is my DArT data and `ZierObco_DZ22-7321_meta.xlsx` is my metadata.

### Metadata

The metadata file has four compulsory columns (sample, site, lat, long), plus any additional data you wish to use. I would recommend having a species or "sp" column as well. You do not need to leave any cells -intentionally blank as with the previous ReCER protocol, as any unwanted samples can be removed during the R analysis by whatever condition you're interested in.

## Installation of RRtools

Obtain a copy of the in house made R package from XXX

    library(devtools)

    install("C:/Users/yaps/Desktop/RRtools")

## Main analyses

### Main analyses

#### `1 data import and qc.Rmd`

-   imports DArT data and metadata
-   filters low quality and fixed loci
-   sets working directory and major packages

<img src="https://user-images.githubusercontent.com/67452867/208334441-7383c64e-61ed-43a1-9782-cc9b80ad5eff.png" alt="qc" width="400"/>

#### `2 create sites by distance.R`

-   groups samples my geographic distance into "sites"
-   good for large or complicated data sets
-   not necessary if you already have clear sites

#### `3 find clones.R`

-   finds individuals that are clones (kinship\>0.45) and removes them from the working dataset

#### `4 kinship heatmap.R`

-   calculates kinship between individuals and plots it with ComplexHeatmap
-   both SNPrelate and popkin package methods are covered

<img src="https://user-images.githubusercontent.com/67452867/208334388-9b6707d8-7e51-4b00-9d48-0449307d8165.png" alt="kin_heatmap" width="400"/>

#### `5 PCA and UMAP.R`

-   runs and plots the results of the dimensionality reduction methods PCA and UMAP
-   needed to identify trends in the data and genetic groups
-   accurate identification of genetic groups is important for selecting the grouping variable in kinship, FST, and diversity statistics analyses

#### `6 SplitsTree and SVDquartets.R`

-   makes output file for visualisation of phylogenetic network in SplitsTree desktop
-   makes output file for SVDquartets analysis in PAUP -- bootstrapped tree with SNP data

#### `7 diversity stats.R`

-   calculates heterozygosity (HO, HE), inbreeding (FIS), and allelic richness (AR) of groups of individuals (usually sites)

#### `8 fst and geographic distance.R`

-   calculates FST and geographic distances between groups of individuals (usually sites)
-   tests the significance of the association between FST and distance using mantel test
-   returns heatmaps of FST and distance <img src="https://user-images.githubusercontent.com/67452867/208334306-1589f924-49ab-4eb2-8635-26cfbfe6b918.png" alt="fst_dist_heatmap" width="600"/> <img src="https://user-images.githubusercontent.com/67452867/208334322-1eff5c93-9953-4045-920e-dac3953d590f.png" alt="mantel" width="400"/>

#### `9 admixture analysis.R`

-   predicts the ancestry of individuals
-   returns admixure plot and scatterpie plot <img src="https://user-images.githubusercontent.com/67452867/208334283-70bf4980-0a8e-41ff-b473-d1c27953de73.png" alt="admixture" width="600"/>

#### `10 maps.R`

-   plots a basic grey geographic map with the sites
-   guidance on how to make satellite map with sites

<img src="https://user-images.githubusercontent.com/67452867/208337159-59ee60af-8a68-4a2e-86b9-ab6569408a74.png" alt="sat_map" width="400"/>

### Accessory analyses

#### `impossible progeny loci.R`

-   only usable for data with seedling and mother relationships
-   determines how many alleles are impossible for the seedling to have based on the mother's alleles

#### `ploidy.R`

-   uses readcount data file supplied by dart
-   makes histograms of the readcount data which can sometimes be used to determine if individuals or groups of individuals have poidy variation

