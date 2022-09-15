## Prepare environment
if (!require("shiny")){
  install.packages("shiny")
}
if (!require("devtools")){
  install.packages("devtools")
}
if (!require("ggplot2")){
  install.packages("ggplot2")
}
if (!require("ggrepel")){
  install.packages("ggrepel")
}
if (!require("data.table")){
  install.packages("data.table")
}
if (!require("DT")){
  install.packages("DT")
}
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("GenomeInfoDbData")){
  BiocManager::install("GenomeInfoDbData")
}

# Get CCprofiler package differential branch
if (!require("CCprofiler")){
  devtools::install_github("CCprofiler/CCprofiler", ref = "differential")
}

# load packages
library(shiny)
library(CCprofiler)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)
library(ggrepel)

### Prepare/read data
# generic
pepTraces <- readRDS("data/pepTracesList_filtered.rds")
protTraces <- readRDS("data/protein_traces_list.rds")
designMatrix <- readRDS("data/design_matrix.rds")
calibration <- readRDS("data/calibration.rds")

# protein-centric
proteinFeatures <- readRDS("data/proteinFeaturesFiltered.rda")
diffProteins_stimulated_undifferentiated <- readRDS("data/protein_DiffExprProtein_stimulated_undifferentiated.rda")
diffProteins_stimulated_differentiated <- readRDS("data/protein_DiffExprProtein_stimulated_differentiated.rda")
diffProteins_differentiated_undifferentiated <- readRDS("data/protein_DiffExprProtein_differentiated_undifferentiated.rda")

diffAssemblyState_stimulated_undifferentiated <- fread("data/diffAssemblyState_stimulated_undifferentiated.txt")
diffAssemblyState_stimulated_differentiated <- fread("data/diffAssemblyState_stimulated_differentiated.txt")
diffAssemblyState_differentiated_undifferentiated <- fread("data/diffAssemblyState_differentiated_undifferentiated.txt")

# complex-centric
complexFeatures <- readRDS("data/complexFeaturesCollapsed.rda")
diffComplexes_stimulated_undifferentiated <- readRDS("data/complex_DiffExprComplex_stimulated_undifferentiated.rda")
diffComplexes_stimulated_differentiated <- readRDS("data/complex_DiffExprComplex_stimulated_differentiated.rda")
diffComplexes_differentiated_undifferentiated <- readRDS("data/complex_DiffExprComplex_differentiated_undifferentiated.rda")

# Explore & process into useable overview tables
proteins = unique(rbind(protTraces$undifferentiated_1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$undifferentiated_2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$undifferentiated_3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$differentiated_1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$differentiated_2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$differentiated_3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$stimulated_1$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$stimulated_2$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)],
                           protTraces$stimulated_3$trace_annotation[, .(protein_id, Entry_name, Gene_names, Entry_name, Length, Mass, protein_mw)]))
# remove redundant Entry_name column
proteins[, 2]=NULL

proteins[, in_diffProteins:=protein_id%in%c(diffProteins_differentiated_undifferentiated$feature_id,
                                            diffProteins_stimulated_differentiated$feature_id,
                                            diffProteins_stimulated_undifferentiated$feature_id), .(protein_id)]

proteins[, in_diffAssembly:=protein_id%in%c(diffAssemblyState_differentiated_undifferentiated$protein_id,
                                            diffAssemblyState_stimulated_differentiated$protein_id,
                                            diffAssemblyState_stimulated_undifferentiated$protein_id), .(protein_id)]

# Add proteins info to diffProteins tables
if (!("Gene_names" %in% names(diffProteins_differentiated_undifferentiated))) {
  diffProteins_differentiated_undifferentiated = merge(diffProteins_differentiated_undifferentiated, proteins[, .(protein_id, Entry_name, Gene_names)],
                                                       by.x = "feature_id", by.y = "protein_id", all.x = T, all.y = F)
saveRDS(diffProteins_differentiated_undifferentiated,
        file = "data/protein_DiffExprProtein_differentiated_undifferentiated.rda")
  }

if (!("Gene_names" %in% names(diffProteins_stimulated_differentiated))) {
  diffProteins_stimulated_differentiated = merge(diffProteins_stimulated_differentiated, proteins[, .(protein_id, Entry_name, Gene_names)],
                                                       by.x = "feature_id", by.y = "protein_id", all.x = T, all.y = F)
  saveRDS(diffProteins_stimulated_differentiated,
          file = "data/protein_DiffExprProtein_stimulated_differentiated.rda")
}

if (!("Gene_names" %in% names(diffProteins_stimulated_undifferentiated))) {
  diffProteins_stimulated_undifferentiated = merge(diffProteins_stimulated_undifferentiated, proteins[, .(protein_id, Entry_name, Gene_names)],
                                                 by.x = "feature_id", by.y = "protein_id", all.x = T, all.y = F)
  saveRDS(diffProteins_stimulated_undifferentiated,
          file = "data/protein_DiffExprProtein_stimulated_undifferentiated.rda")
}

# Add proteins info to diffAssembly tables
if (!("Gene_names" %in% names(diffAssemblyState_differentiated_undifferentiated))) {
  diffAssemblyState_differentiated_undifferentiated = merge(diffAssemblyState_differentiated_undifferentiated, proteins[, .(protein_id, Entry_name, Gene_names)],
                                                       by = "protein_id", all.x = T, all.y = F)
  saveRDS(diffAssemblyState_differentiated_undifferentiated,
            file = "data/diffAssemblyState_differentiated_undifferentiated.rda")
}

if (!("Gene_names" %in% names(diffAssemblyState_stimulated_differentiated))) {
  diffAssemblyState_stimulated_differentiated = merge(diffAssemblyState_stimulated_differentiated, proteins[, .(protein_id, Entry_name, Gene_names)],
                                                            by = "protein_id", all.x = T, all.y = F)
  saveRDS(diffAssemblyState_stimulated_differentiated,
          file = "data/diffAssemblyState_stimulated_differentiated.rda")
}

if (!("Gene_names" %in% names(diffAssemblyState_stimulated_undifferentiated))) {
  diffAssemblyState_stimulated_undifferentiated = merge(diffAssemblyState_stimulated_undifferentiated, proteins[, .(protein_id, Entry_name, Gene_names)],
                                                            by = "protein_id", all.x = T, all.y = F)
  saveRDS(diffAssemblyState_stimulated_undifferentiated,
          file = "data/diffAssemblyState_stimulated_undifferentiated.rda")
}


# Annotate parent complexes
proteins[, parent_complex_ids:=paste(unique(complexFeatures[grep(protein_id, complexFeatures$subunits)]$complex_id), collapse = ","), .(protein_id)]
proteins[, parent_complex_names:=paste(unique(complexFeatures[grep(protein_id, complexFeatures$subunits)]$complex_name), collapse = ","), .(protein_id)]
proteins[, parent_complex_feature_identifiers:=paste(unique(complexFeatures[grep(protein_id, complexFeatures$subunits)]$unique_feature_identifier), collapse = ","), .(protein_id)]

saveRDS(proteins, "data/proteins.rda")
saveRDS(names(proteins)[c(1:3,9,10,11)], "data/protidcols.rda")

# Make complex selection table (Differential Tab)
complexes_diff_ids = unique(c(diffComplexes_differentiated_undifferentiated$complex_id,
                              diffComplexes_stimulated_differentiated$complex_id,
                              diffComplexes_stimulated_undifferentiated$complex_id))
saveRDS(complexes_diff_ids, file = "data/complexes_diff_ids.rda")

# Make complex feature selection table
names(complexFeatures)
cfidcols = names(complexFeatures)[c(7,1,2,3)]
saveRDS(cfidcols, file = "data/cfidcols.rda")

# Set bioconductor repository for Shiny deploy to go thru..
library(BiocManager)
options(repos = BiocManager::repositories())
