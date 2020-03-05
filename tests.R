default_proteins <- c("NUP54", "NUP62", "NUP58 KIAA0410 NUPL1")

input = list("fcolumn" = "Gene_names",
             "fvalue" = default_proteins,
             "logscale" = FALSE,
             "collapse_conditions" = FALSE,
             "collapse_replicates" = TRUE,
             "show_monomers" = TRUE,
             "pc_volcano_pvalcutoff" = 0.05,
             "pc_volcano_fccutoff" = 2,
             "cc_volcano_pvalcutoff" = 0.05,
             "cc_volcano_fccutoff" = 2,
             "pc_assemblyScatter_meanDiffcutoff" = 0.3)

target_id_traces = subset(protTraces, trace_subset_ids = input$fvalue,
                          trace_subset_type = input$fcolumn)

plot(target_id_traces,
     #colour_by = input$fcolumn,
     collapse_conditions = input$collapse_conditions,
     aggregateReplicates = input$collapse_replicates,
     name = "",
     monomer_MW = input$show_monomers,
    # log = input$logscale,
     design_matrix = designMatrix,
     plot = FALSE)

# get protein info for selected set independent of filter column type
proteins[get(input$fcolumn) %in% input$fvalue]

## Tabs
# '1) View protein traces'
# '2) Protein differential intensity'
# '3) Protein differential assembly'
# '4) Complex differential intensity'

## Objects/plots
# plot
# anntable
# pc_volcano
# pc_difftable
# pc_assemblyScatter
# pc_assemblyTable
# cc_volcano
# cc_traces

##############
# Volcano tests

diffProteins = diffProteins_differentiated_undifferentiated
 
  selected_protein_ids = proteins[get(input$fcolumn) %in% input$fvalue]$protein_id
  
  dplot <<- plotVolcano(diffProteins,
                        pBHadj_cutoff = input$pc_volcano_pvalcutoff,
                        FC_cutoff = input$pc_volcano_fccutoff,
                        highlight = selected_protein_ids,
                        plot = FALSE)
  ggplotly(dplot)


# Scatter test
diffProteinAssemblyState = copy(diffAssemblyState_differentiated_undifferentiated)

selected_protein_ids = proteins[get(input$fcolumn) %in% input$fvalue]$protein_id

meanDiff_cutoff = input$pc_assemblyScatter_meanDiffcutoff


splot1 <<- ggplot(diffProteinAssemblyState,
                 aes(x=meanAMF1, y=meanAMF2, colour=-log10(betaPval_BHadj), label = paste(protein_id,Entry_name,Gene_names))) +
  geom_abline(intercept = meanDiff_cutoff, slope = 1) +
  geom_abline(intercept = -meanDiff_cutoff, slope = 1) +
  geom_point() +
  theme_bw() 

splot2 = splot1 + geom_point(data = diffProteinAssemblyState[protein_id %in% selected_protein_ids],
             aes(x=meanAMF1, y=meanAMF2, label = paste(protein_id,Entry_name,Gene_names)),
             colour="red", size = 3)
ggplotly(splot2)


## Complex-centric tests
library(UpSetR)
upset(fromList(list("d_ud" = diffComplexes_differentiated_undifferentiated$complex_id,
                    "st_d" = diffComplexes_stimulated_differentiated$complex_id,
                    "st_ud" = diffComplexes_stimulated_undifferentiated$complex_id)))
# All are the same!

# CC volcano
diffComplexes = copy(diffComplexes_differentiated_undifferentiated)
selected_complex_id = "127_corum_corum"

plotVolcano(diffComplexes,
            pBHadj_cutoff = input$cc_volcano_pvalcutoff,
            FC_cutoff = input$cc_volcano_fccutoff,
            highlight = selected_complex_id,
            plot = FALSE)

# debugonce(plotVolcano)
# label choice n/a, workaround: Additional function plotVolcano_c for complex level plot with hover label


## Test cc feature plots
target_complex = selected_complex_id

# Test alternative complex ids that are part of a complexfeatureid group
# target_complex = "2202_corum_corum" -> runs itself dead
# full feature id:
input[["cf_column"]] = "complex_id"
input[["cf_value"]] = complexFeaturesCollapsed$complex_id[11]

names(complexFeaturesCollapsed)

complex_feature_complex_ids = complexFeaturesCollapsed$complex_id[11] # works smoothly

CCprofiler::plotFeatures(feature_table = complexFeaturesCollapsed,
             traces = protTraces,
             feature_id = input$cf_value,
             design_matrix=designMatrix,
             calibration=calibration_functions,
             annotation_label = "Entry_name",
             peak_area = T,
             legend = F,
             onlyBest = F,
             PDF = FALSE,
             monomer_MW=T,
             aggregateReplicates=T)


p = plotFeaturesObject(feature_table = complexFeaturesCollapsed,
             traces = protTraces,
             feature_id = input$cf_value,
             design_matrix=designMatrix,
             calibration=calibration_functions,
             annotation_label = "Entry_name",
             peak_area = T,
             legend = F,
             onlyBest = F,
             PDF = FALSE,
             monomer_MW=T,
             aggregateReplicates=T)
# plot object return won't work; internal methods required:
# Fehler in plotFeaturesObject(feature_table = complexFeaturesCollapsed,  : 
# unbenutzte Argumente (design_matrix = designMatrix, aggregateReplicates = T)





