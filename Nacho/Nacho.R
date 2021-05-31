#########
### this script is intented to analyse Nanostring data
#########
##################
### loading required packages
##################

if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("NACHO")) {
  install.packages("NACHO", ask =FALSE)
  library("NACHO")
}
if (!require("GEOquery")) {
  BiocManager::install("GEOquery", dependencies = TRUE)
  library("GEOquery")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library("dplyr")
}
if (!require("tidyr")) {
  BiocManager::install("tidyr", dependencies = TRUE)
  library("tidyr")
}
if (!require("tibble")) {
  BiocManager::install("tibble", dependencies = TRUE)
  library("tibble")
}
if (!require("limma")) {
  BiocManager::install("limma", dependencies = TRUE)
  library("limma")
}




vignette("NACHO-analysis")
Maja_Covid <- load_rcc(data_directory = "/media/rmejia/mountme88/Projects/Maja-covid/Data/Original_RCC_files/", ssheet_csv ="/media/rmejia/mountme88/Projects/Maja-covid/Data/ssheet_csv.csv", id_colname = "uniqueID" )

# visualise( covidnma_01_01)

Maja_Covid_normalized <- normalise(nacho_object=Maja_Covid, housekeeping_norm = TRUE, normalisation_method = "GEO",remove_outliers = TRUE)

Maja_Covid_normalized$pc_sum
Maja_Covid_normalized$nacho
dim(Maja_Covid_normalized$nacho)
#render(Maja_Covid_normalized, output_dir = "")

######



#######

data_directory <- file.path("/media/rmejia/mountme88/Projects/Maja-covid/Data/NACHO_test", "GSE70970", "Data")
dir.create(data_directory,recursive = TRUE)
gse <- getGEO("GSE70970")
targets <- pData(phenoData(gse[[1]]))
getGEOSuppFiles(GEO = "GSE70970", baseDir = "/media/rmejia/mountme88/Projects/Maja-covid/Data/NACHO_test")
untar(
  tarfile = file.path( "/media/rmejia/mountme88/Projects/Maja-covid/Data/NACHO_test", "GSE70970", "GSE70970_RAW.tar"), 
  exdir = data_directory
)
targets$IDFILE <- list.files(data_directory)                
GSE70970 <- load_rcc(data_directory, targets, id_colname = "IDFILE")


### GEt the phenotypes
selected_pheno <- GSE70970[["nacho"]] %>% 
  select(IDFILE, `age:ch1`, `gender:ch1`, `chemo:ch1`, `disease.event:ch1`) %>% 
  distinct() %>% 
  mutate_all(~ na_if(.x, "NA")) %>% 
  drop_na()

## get the normalised counts
expr_counts <- GSE70970[["nacho"]] %>% 
  filter(grepl("Endogenous", CodeClass)) %>% 
  select(IDFILE, Name, Count_Norm) %>% 
  pivot_wider(names_from = "Name", values_from = "Count_Norm") %>% 
  column_to_rownames("IDFILE") %>% 
  t()
## Select the phenotypes and counts
samples_kept <- intersect(selected_pheno[["IDFILE"]], colnames(expr_counts))
expr_counts <- expr_counts[, samples_kept]
selected_pheno <- filter(selected_pheno, IDFILE %in% !!samples_kept)

## Design matrix
design <- model.matrix(~ `disease.event:ch1`, selected_pheno)

# limma
results <- eBayes(lmFit(expr_counts, design))

topTable(results)
