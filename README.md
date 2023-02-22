# Overview of our DArTseq analyses

[This README documentation is adopted from Eilish's README file: <https://github.com/eilishmcmaster/BaseSoSAnalyses>]

This git repository explains our analytical workflow of DArTseq data that was designed to streamline the production of population genetic outputs from a range of common Australian species in New South Wales as part of our flagship Restore and Renew project, and this workflow have since been used by ReCER for a range of large scale multispecies studies.

Here is a general overview of how to set up files, the type of analyses and outputs we obtain and how we interpret the results.

The analyses are split into individual R files so that more in depth explaining can be done for each analysis. The files have been numbered to make it easier to navigate through. Typically we combine the scripts into a single file to keep track of the analyses that are run on a specific dataset.

## DArTseq data

### What do we get?

DArTseq data is returned to us as a zipped folder containing multiple spreadsheets. There are multiple spreadsheets because the data be SNP data formatted in single or double rows or SilicoDArT data. We only use the spreadsheet with SNP data in single row format for our work.

Note that DArT does not label their files in a consistent manner, the file name can be in a .xlsx or .csv format and can end with a "SNP_mapping_2.csv" or "SNP_singlerow.csv".

Always check that you have the correct file:

The sequence data produced by DArT goes through their proprietary filters and the resulting spreadsheets contain loci ID and relevant quality information (e.g. call rate, repeatability, information content). The genotype data is further to the right of the spreadsheet. Refer to the "metada.xlsx" file accompanying the spreadsheets for the type of quality information captured in the spreadsheets and the samples that were successfully sequenced. The SNP data formatted in single rows does not contain duplicated IDs and the data contains 0s (homozygote), 1s (homozygote), 2s (heterozygote).

Note that when a dataset is received, it is up to the person in charge of the study to save the dataset on N drive under: N:\\Evolutionary Ecology\\Restore & Renew\\DArT data

### Analysing the data

An R package called RRtools was developed in house as part of the Restore and Renew analytical workflow. The package was designed to allow the analysis of DArTseq data for a large number of species in a streamline manner (i.e. consistent labelling system and structured working directories).

### Setting up your directories

Prior to analysing, a working directory needs to be set up. The directory is structured in a way that allows RRtools, our in house developed R package to access the DArTseq data and output the data in various formats.

Below is an example working directory for the *Zieria obcordata* dataset.

The name of the working directory is usually a combination of the abbreviated genus and species name because the DArTseq data is mostly used for studying a species. e.g., ZierObco for Zieria obcordata.

In the directory, three folders are to be manually generated (dart_raw, meta, popgen) and two folders (dart_standard, outputs) will appear as RRtools is run:

    |-- ZierObco ........................working directory
        |-- dart_raw ........................where DArTseq data is kept
            |-- Report_DZ22-7321_SNP_mapping_2.csv
        |-- meta ........................where metadata is kept
            |-- ZierObco_DZ22-7321_meta.csv
        |-- popgen ...................where DArTseq data in different formats are kept
        |-- dart_standard ....... where saved objects containing DArTseq data are kept  
        |-- outputs....................  where visual outputs are kept

Where `Report_DZ22-7321_SNP_mapping_2.csv` is the DArTseq data in single row format and `ZierObco_DZ22-7321_meta.xlsx` is themetadata.

### Metadata

The metadata is in a spreadsheet with four compulsory columns (sample, site, lat, long), plus any additional data you wish to use. This can be a species or "sp" column if outgroup species are included in the dataset. A notes column is recommended to keep track of the samples during the analysis.

You do not need to leave any cells -intentionally blank as with the previous ReCER protocol, as any unwanted samples can be removed during the R analysis by whatever condition you're interested in.

## Installing RRtools

The RRtools workflow consist of an initial quality SNP filtering step and subsequent conversion of the genotype data into a range of formats depending on the type of analysis to be conducted.

RRtools can only be installed from a local directory (i.e., not from an online repository). A copy of the R package can be obtained from N:\\Evolutionary Ecology\\R for everyone.

Example of how to install RRtools in R below:

    library(devtools)

    install("C:/Users/yaps/Desktop/RRtools")

## Using RRtools

### Quality filtering

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
