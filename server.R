####################################################################################
# SECexplorer
####################################################################################
# Server

## prepare environment
if (!require("shiny")){
  install.packages("shiny")
}
if (!require("ggplot2")){
  install.packages("ggplot2")
}
if (!require("plotly")){
  install.packages("plotly")
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
library(ggplot2)
library(plotly)
library(data.table)
library(DT)
library(CCprofiler)

# load modified methods
source("methods.R")

## prepare data
#calibration_functions <- readRDS("www/data/calibration_functions.rda")
calibration_functions <- readRDS("data/calibration.rds")
proteins <- readRDS("data/proteins.rda")
protTraces = readRDS("data/protein_traces_list.rds")
designMatrix = readRDS("data/design_matrix.rds")
# default_proteins <- c("GPS1 COPS1 CSN1", "COPS3 CSN3", "COPS8 CSN8")
default_proteins <- c("NDC80 HEC HEC1 KNTC2",
                      "SPC24 SPBC24",
                      "NUF2 CDCA1 NUF2R",
                      "SPC25 SPBC25 AD024")

default_complexftid <- "127_corum_corum"

# Differentials
# Protein-level
diffProteins_differentiated_undifferentiated <- readRDS("data/protein_DiffExprProtein_differentiated_undifferentiated.rda")
diffProteins_stimulated_undifferentiated <- readRDS("data/protein_DiffExprProtein_stimulated_undifferentiated.rda")
diffProteins_stimulated_differentiated <- readRDS("data/protein_DiffExprProtein_stimulated_differentiated.rda")

diffAssemblyState_stimulated_undifferentiated <- readRDS("data/diffAssemblyState_stimulated_undifferentiated.rda")
diffAssemblyState_stimulated_differentiated <- readRDS("data/diffAssemblyState_stimulated_differentiated.rda")
diffAssemblyState_differentiated_undifferentiated <- readRDS("data/diffAssemblyState_differentiated_undifferentiated.rda")

# Complex-level
complexFeaturesCollapsed = readRDS("data/complexFeaturesCollapsed.rda")
diffComplexes_stimulated_undifferentiated <- readRDS("data/complex_DiffExprComplex_stimulated_undifferentiated.rda")
diffComplexes_stimulated_differentiated <- readRDS("data/complex_DiffExprComplex_stimulated_differentiated.rda")
diffComplexes_differentiated_undifferentiated <- readRDS("data/complex_DiffExprComplex_differentiated_undifferentiated.rda")

## define server roles
#######################

shinyServer(function(input, output, session) {
  
  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  # protein selection
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(proteins[[input$fcolumn]]))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Search and select proteins of interest", values,
                   multiple = TRUE, options = list(maxOptions = 6000),
                   selected = default_proteins)
  })
  
  # complex feature selection
  output$cfcolumnvalues <- renderUI({
    values <- sort(unique(complexFeaturesCollapsed[[input$cfcolumn]]))
    selectizeInput("cfvalue", "Search and select complex features of interest", values,
                   multiple = FALSE, options = list(maxOptions = 6000),
                   selected = default_complexftid)
  })
  
  ############################
  ## Viewer Tab              #
  ############################

  ## generate selected protein SEC traces plot
  # Subset traces
  target_id_traces <- eventReactive(input$fvalue,{
    
    selected_protein_ids = proteins[get(input$fcolumn) %in% input$fvalue]$protein_id
    
    target_id_traces = subset(protTraces, trace_subset_ids = selected_protein_ids,
                             trace_subset_type = "id")
  })

  ## Plot the selected traces
  output$plot <- renderPlotly({
    vplot <<- plot(target_id_traces(),
                   # colour_by = input$fcolumn, ## causes problems in combination with collapsing
                   collapse_conditions = input$collapse_conditions,
                   aggregateReplicates = input$collapse_replicates,
                   name = "",
                   monomer_MW = input$show_monomers,
                   log = input$logscale,
                   design_matrix = designMatrix,
                   plot = FALSE)
    ggplotly(vplot)
  })
  
  ## Download the displayed plot
  output$downloadPlot <- downloadHandler(
    filename = function() { paste("currentPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, width=10, height=6, plot = vplot, device = "pdf")
    }
  )
  
  # Display the annotation table for the selected proteins
  output$anntable <- renderDT({
    proteins[get(input$fcolumn) %in% input$fvalue]
  })
  
  #####################################
  ## Differential protein intensity   #
  #####################################
  
  # Select dataset based on user-defined comparison
  
  # choices= c("Differentiated vs. undifferentiated",
  #            "Stimulated vs. differentiated",
  #            "Stimulated vs. undifferentiated")
  
  diffProteins <- eventReactive(input$comparison_pINT,{
    
    if (input$comparison_pINT == "Differentiated vs. undifferentiated"){
      diffProteins = diffProteins_differentiated_undifferentiated
    } else if (input$comparison_pINT == "Stimulated vs. differentiated"){
      diffProteins = diffProteins_stimulated_differentiated
    } else {
      diffProteins = diffProteins_stimulated_undifferentiated
    }
    
  })
  
  # render pc volcano
  output$pc_volcano <- renderPlotly({
    
    selected_protein_ids = proteins[get(input$fcolumn) %in% input$fvalue]$protein_id
    
    dplot <<- plotVolcano(diffProteins(),
                          pBHadj_cutoff = input$pc_volcano_pvalcutoff,
                          FC_cutoff = input$pc_volcano_fccutoff,
                          highlight = selected_protein_ids,
                          plot = FALSE)
    ggplotly(dplot)
  })
  
  # render pc diff table
  output$pc_difftable <- DT::renderDT({
    diffProteins.s = diffProteins()[, .(feature_id, Entry_name, Gene_names, Npeptides,
                                      apex, pBHadj, medianLog2FC, qVal, global_pBHadj, global_medianLog2FC, global_qVal)]
    if(input$pc_difftable_show_all){
      diffProteins.s
    } else {
      diffProteins.s[pBHadj <= input$pc_volcano_pvalcutoff][abs(medianLog2FC) >= 
                                                          log2(input$pc_volcano_fccutoff)]
    }
  })
  
  #####################################
  ## Differential protein assembly    #
  #####################################
  
  # Select dataset based on comparison
  diffProteinAssemblyState <- eventReactive(input$comparison_pAMF,{
    
    if (input$comparison_pAMF == "Differentiated vs. undifferentiated"){
      diffProteinAssemblyState = diffAssemblyState_differentiated_undifferentiated
    } else if (input$comparison_pAMF == "Stimulated vs. differentiated"){
      diffProteinAssemblyState = diffAssemblyState_stimulated_differentiated
    } else {
      diffProteinAssemblyState = diffAssemblyState_stimulated_undifferentiated
    }
    
  })
  
  # render assembly state scatter plot
  output$pc_assemblyScatter <- renderPlotly({
    
    selected_protein_ids = proteins[get(input$fcolumn) %in% input$fvalue]$protein_id
    
    meanDiff_cutoff = input$pc_assemblyScatter_meanDiffcutoff
    
    splot1 <<- ggplot(diffProteinAssemblyState(),
                      aes(x=meanAMF1, y=meanAMF2, colour=-log10(betaPval_BHadj), label = paste(protein_id,Entry_name,Gene_names))) +
      geom_abline(intercept = meanDiff_cutoff, slope = 1) +
      geom_abline(intercept = -meanDiff_cutoff, slope = 1) +
      geom_point() +
      theme_bw() 
    
    splot2 = splot1 + geom_point(data = diffProteinAssemblyState()[protein_id %in% selected_protein_ids],
                                 aes(x=meanAMF1, y=meanAMF2, label = paste(protein_id,Entry_name,Gene_names)),
                                 colour="red", size = 3)
    ggplotly(splot2)
  })
  
  # render assembly state output table
  output$pc_assemblyTable <- DT::renderDT({
    if(input$pc_diffAssemblytable_show_all){
      dt = diffProteinAssemblyState()
      dt[, more_assembled_in:=NA]
      dt[meanDiff >= input$pc_assemblyScatter_meanDiffcutoff, more_assembled_in:=2]
      dt[meanDiff <= -input$pc_assemblyScatter_meanDiffcutoff, more_assembled_in:=1]
      dt
      
    } else {
      rbind(diffProteinAssemblyState()[meanDiff >= input$pc_assemblyScatter_meanDiffcutoff, more_assembled_in:=2],
            diffProteinAssemblyState()[meanDiff <= -input$pc_assemblyScatter_meanDiffcutoff, more_assembled_in:=1])
    }
  })
  
  #####################################
  ## Complex feature viewer           #
  #####################################
  
  # render complex feature plot
  # Note: This plot is non-interactive as export of plot object from base function
  # plotFeatures (in combination with tracesList input to traces arg) doesn't fly
  output$cf_plot <- renderPlot({
    
    selected_complex_ft_id = complexFeaturesCollapsed[get(input$cfcolumn) %in% input$cfvalue]$complex_id
    
    CCprofiler::plotFeatures(feature_table = complexFeaturesCollapsed,
                 traces = protTraces,
                 feature_id = selected_complex_ft_id,
                 design_matrix=designMatrix,
                 calibration=calibration_functions,
                 annotation_label = "Entry_name",
                 peak_area = T,
                 legend = F,
                 onlyBest = F,
                 PDF = FALSE,
                 monomer_MW=T,
                 aggregateReplicates=T)
    
  })
  
  # re-render plot for download (no plot object available)
  output$downloadPlot_cF <- downloadHandler(
    filename = function() { paste("complexFeatures_", input$cfvalue, '.pdf', sep='') },
    content = function(file) {ggsave(file, width=10, height=6, 
                                     plot = CCprofiler::plotFeatures(feature_table = complexFeaturesCollapsed,
                                                                     traces = protTraces,
                                                                     feature_id = complexFeaturesCollapsed[get(input$cfcolumn) %in% input$cfvalue]$complex_id,
                                                                     design_matrix=designMatrix,
                                                                     calibration=calibration_functions,
                                                                     annotation_label = "Entry_name",
                                                                     peak_area = T,
                                                                     legend = F,
                                                                     onlyBest = F,
                                                                     PDF = FALSE,
                                                                     monomer_MW=T,
                                                                     aggregateReplicates=T), 
                                     device = "pdf")}
    )
  
  # render complex feature table
  output$cf_table <- DT::renderDT({
    if(input$cf_table_show_all){
      complexFeaturesCollapsed
    } else {
      complexFeaturesCollapsed[get(input$cfcolumn) %in% input$cfvalue]
    }
  }
  )
  
  #####################################
  ## Differential complex intensity   #
  #####################################
  
  # Select dataset based on comparison
  diffComplexes <- eventReactive(input$comparison_cINT,{
    
    if (input$comparison_cINT == "Differentiated vs. undifferentiated"){
      diffComplexes = diffComplexes_differentiated_undifferentiated
    } else if (input$comparison_cINT == "Stimulated vs. differentiated"){
      diffComplexes = diffComplexes_stimulated_differentiated
    } else {
      diffComplexes = diffComplexes_stimulated_undifferentiated
    }
    
  })
  
  # Render complex volcano
  output$cc_volcano <- renderPlotly({
    
    selected_complex_id = input$complexid
    
    cvplot <<- plotVolcano_c(diffComplexes(),
                          pBHadj_cutoff = input$cc_volcano_pvalcutoff,
                          FC_cutoff = input$cc_volcano_fccutoff,
                          highlight = selected_complex_id,
                          plot = FALSE)
    ggplotly(cvplot)
  })
  
  # Render complex table
  output$cc_difftable <- DT::renderDT({
    if(input$cc_difftable_show_all){
      diffComplexes()
    } else {
      diffComplexes()[pBHadj <= input$cc_volcano_pvalcutoff][abs(medianLog2FC) >= 
                                                             log2(input$cc_volcano_fccutoff)]
    }
  }
  #, options = list(autoWidth = TRUE,
  #                  columnDefs = list(list(width = '100px', targets = 1))
  #                  )
  # attempt to limit column width failed.. move on
  )
  
})



