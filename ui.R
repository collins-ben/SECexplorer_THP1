####################################################################################
# SECexplorer
####################################################################################
# UI

## prepare environment
if (!require("shiny")){
  install.packages("shiny")
}
if (!require("shinythemes")){
  install.packages("shinythemes")
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
library(shinythemes)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)
library(CCprofiler)

## prepare data
idcols <- readRDS("data/protidcols.rda")
complexes_diff_ids <- readRDS("data/complexes_diff_ids.rda")
cfidcols <- readRDS("data/cfidcols.rda")

## define user interface
########################
shinyUI(fluidPage(
  
  # Theme
  theme = shinytheme("cosmo"),

  # Application title
  headerPanel("SECexplorer-THP1", windowTitle = "SECexplorer-THP1"),

  # Sidebar with a slider input for number of observations
  sidebarLayout(
    sidebarPanel(
      width = 3,
      p(),
      selectInput(inputId = "fcolumn",
                  label = "Select proteins based on:",
                  choices = idcols,
                  selected = "Gene_names"),
      uiOutput("fcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
      p("Delete above entries by backspace and start typing for live search for your target protein(s)"),

      conditionalPanel('input.tabid === "1) View protein traces"',

                       checkboxInput("collapse_conditions", label = "Collapse conditions",
                                     value = FALSE),
                       checkboxInput("collapse_replicates", label = "Collapse replicates",
                                     value = TRUE),
                       checkboxInput("logscale", "LOG-transform Y-Axis", value = FALSE),
                       checkboxInput("show_monomers", "Indicate monomer expected fractions", value = TRUE),
                       # checkboxInput("error_bars", "Plot error area", value = TRUE),
                       uiOutput("errortype")
      ),
      conditionalPanel('input.tabid === "2) Protein differential intensity"',
                       p("Identify outlier proteins from plot hover info. To analyze traces, add the proteins of interest via the fields above."),
                       selectizeInput("comparison_pINT",
                                      label = "Select comparison",
                                      choices= c("Differentiated vs. undifferentiated",
                                                 "Stimulated vs. differentiated",
                                                 "Stimulated vs. undifferentiated"),
                                      selected = "Differentiated vs. undifferentiated",
                                      multiple = FALSE),
                       numericInput("pc_volcano_pvalcutoff", label = "pBHadj cutoff:", value = 0.05, max = 1),
                       numericInput("pc_volcano_fccutoff", label = "pBHadj cutoff:", value = 2),
                       checkboxInput("pc_difftable_show_all", "Show all (i.e. also non-significant)\n entries in table", value = FALSE)
                       
      ),
      conditionalPanel('input.tabid === "3) Protein differential assembly"',
                       p("Identify outlier proteins from plot hover info. To analyze traces, add the proteins of interest via the fields above."),
                       selectizeInput("comparison_pAMF",
                                      label = "Select comparison",
                                      choices= c("Differentiated vs. undifferentiated",
                                                 "Stimulated vs. differentiated",
                                                 "Stimulated vs. undifferentiated"),
                                      selected = "Differentiated vs. undifferentiated",
                                      multiple = FALSE),
                       numericInput("pc_assemblyScatter_meanDiffcutoff", label = "mean Difference cutoff", value = 0.25),
                       checkboxInput("pc_diffAssemblytable_show_all",
                                     "Show all entries in table\n(i.e. not only those passing mean Difference cutoff criterion)", value = FALSE)
                       
      ),
      conditionalPanel('input.tabid === "4) View complex features"',
                      selectInput(inputId = "cfcolumn",
                                  label = "Select complex features based on:",
                                  choices = cfidcols,
                                  selected = "complex_id"),
                      uiOutput("cfcolumnvalues"), #The possible choices of this field are calculated on the server side and passed over by uiOutput
                      p("Select complex feature of interest above"),
                      p(),
                      checkboxInput("cf_table_show_all",
                                    "Show all entries in table\n(i.e. not only the selected feature id)", value = FALSE)
                      
                       
      ),
      conditionalPanel('input.tabid === "5) Complex differential intensity"',
                       selectizeInput("comparison_cINT",
                                      label = "Select comparison",
                                      choices= c("Differentiated vs. undifferentiated",
                                                 "Stimulated vs. differentiated",
                                                 "Stimulated vs. undifferentiated"),
                                      selected = "Stimulated vs. differentiated",
                                      multiple = FALSE),
                       p("Identify regulated complexes from plot hover info.\n
                         To analyze traces, select complex of interest via the field below:"),
                       selectizeInput(inputId = "complexid",
                                      label = "Search and select complex of interest",
                                      choices = complexes_diff_ids,
                                      multiple = FALSE, 
                                      # options = list(maxOptions = 6000),
                                      selected = "1172_string;624_string"),
                       numericInput("cc_volcano_pvalcutoff", label = "maximal pBHadj cutoff:", value = 0.05, max = 1),
                       numericInput("cc_volcano_fccutoff", label = "minimal fold-change (FC) cutoff:", value = 2),
                       checkboxInput("cc_difftable_show_all", "Show all (i.e. also non-significant)\n entries in table", value = FALSE)
                       
      ),
      br(),
      br(),
      br(),
    ),

    # Instruction and analysis panels
    mainPanel(
      tabsetPanel(
        id = 'tabid',
        tabPanel('Welcome & Instructions',
                 h1("Welcome to SECexplorer THP1"),
                 p("SECexplorer offers a dynamic interface for the investigator-driven exploration of mass spectrometric SEC-SWATH-MS datasets describing proteome assembly state dynamics as a function of different perturbed states of the biological system."),
                 p("In this study, we established an integrated experimental and computational workflow for rapid profiling of protein complex re-organization in perturbed systems"),
                 p(),
                 p("Bludau, I., Nicod, C., Martelli, C., Xue, P., Heusel, M., Fossati, A., Uliana, F., Frommelt, F., Aebersold, R. & Collins, B. C. Rapid profiling of protein complex re-organization in perturbed systems. bioRxiv 2021.12.17.473177 (2021) https://doi.org/10.1101/2021.12.17.473177"),
                 p(),
                 p("THP1 cells were investigated in three distinct biological conditions:"),
                 tags$ol(
                   tags$ul("Undifferentiated: THP1 cells were grown in suspension in the monocytic state"),
                   tags$ul("Differentiated: Monocytic cells were differentiated in to a macrophage-like phenotype using phorbol ester treatment"),
                   tags$ul("Stimulated: Macrophage cells were stimulated with lipopolysaccharide")),
                 p(),
                 p("This tool facilitates exploration of this dataset in the following analysis tabs:"),
                 tags$ol(
                   tags$ul("1) View protein traces: Allows interactive viewing of protein SEC fractionation profiles across conditions."),
                   tags$ul("2) Protein differential intensity: View differential scores on protein level across conditions."),
                   tags$ul("3) Protein differential assembly: View proteins differential assembled mass fraction across conditions."),
                   tags$ul("4) View complex features: View detected complex signals across conditions."),
                   tags$ul("5) Complex differential intensity: View differential scores on complex level across conditions.")
                 ),
                 h1("Usage instructions:"),
                 p("To start an analysis, select the respective analysis tab above. Read this first, though or come back here for instructions if you get stuck :-)."),
                 p(),
                 
                
                  h2("1) View protein traces"),
                 img(src='1.png', align = "left", width = "100%"),
                 h3("1. Selection of protein identifier type or attribute on which to select proteins"),
                 tags$ul(
                   tags$li("Entry_name or Gene_Names are most informative and intuitive to search."),
                   tags$li("Alternatively, proteins can be selected based on their participation in a given protein complex feature."),
                   tags$li("Searching for specific gene/protein identifiers, or complex features they are associated with, is then possible in field 2.")
                   ),
                 h3("2. Search and selection of multiple proteins"),
                 tags$ul(
                   tags$li("Searching is achieved by deleting the current entry with backspace and starting to 
                 type. All identifiers of the type selected in (1) will be searched for the string 
                     entered, on-the-fly, with potential results showing up below the field, selectable 
                     by <Enter> or by a left mouse click."),
                   tags$li("NOTE: Individual proteins/selections can be removed again by <Backspace>.")
                 ),
                 h3("3. Protein chromatogram view across conditions"),
                 tags$ul(
                   tags$li("By default, a split graph shows the protein level abundance profiles over the 
                            chromatographic fractions in either condition."),
                   tags$li("The chromatograms are displayed as mean of the three replicates per condition."),
                   tags$li("Options for chromatogram display can be selected in (4) and a more detailed pdf version of the interactive graph
                           can be downloaded via the download button (5).")
                 ),
                 p("Note that the plots are interactive! (Thanks, plotly!)"),
                 h3("4. Options for chromatogram display"),
                 tags$ul(
                   tags$li("Selection whether conditions or replicates should be collapsed in the plot. 
                           \n Collapsing replicates produces mean intensities across replicates.
                           \n Collapsing conditions produces line-type-encoded conditions."),
                   tags$li("Option for LOG10-transformation of the intensity axis to spot low-abundant protein pools"),
                   tags$li("Selection whether monomer expected fraction markers shall be displayed. 
                           \n Only applies to the downloadable .pdf") #,
                   # tags$li("Option whether (and which type of) error areas shall be indicated. Works best on split plot.")
                   ),
                 
                 h3("5. Download current graph as .pdf"),
                 tags$ul(
                   tags$li("Allows download of a more detailed pdf version of the current graph.")
                 ),
                 h3("6. Annotation table."),
                 tags$ul(
                   tags$li("Contains metainformation on the currently selected protein set."),
                   tags$li("Also contains information in which complex id and complex feature id the protein is present.
                           Detected complex features (including co-complex member subunits) can be selected and visualized in tab 4)!")
                 ),
                 p(),
                 
                 
                 
                 h2("2) Protein differential intensity"),
                 img(src='2.png', align = "left", width = "100%"),
                 h3("1. Volcano plot of protein-feature-level differential scores of the comparison selected in (2)"),
                 tags$ul(tags$li("Features of currently selected proteins are highlighted in the plot.
                                 \nFeatures passing the significance cutoffs chosen in (3) are highlighted in red.
                                 \nFeatures with insignificant differences are highlighted in dark grey."),
                         tags$li("Outlier proteins identified from the hover information can be added to the
                         analysis via the regular input fields.")),
                 h3("2. Select pairwise comparison of biological conditions"),
                 tags$ul(tags$li("Choose the comparison of biological conditions for the differential analysis.")),
                 h3("3. Select significance cut-offs"),
                 tags$ul(
                   tags$li("Select max. pBHadj and min. fold-change (FC) cut-offs."),
                   tags$li("Selection will affect marker color in the volcano plot (1) and the contents of the table output in dependence on switch (4).")
                   ),
                 h3("4. Show all entries in protein-level differential table"),
                 tags$ul(tags$li("Select whether or not 'significant' features are reported in the output table (5)")),
                 h3("5. Protein-level differential table"),
                 tags$ul(tags$li("Protein feature level differential scores for the given comparison.
                                 \n Displaying scores of only 'significant' features (default) or all scores (switch (4)).")),
                 p(),
                 
                 
                 
                 h2("3) Protein differential assembly"),
                 img(src='3.png', align = "left", width = "100%"),
                 h3("1. Differential assembled mass fraction graph with scores of selected proteins highlighted."),
                 tags$li("Outlier proteins identified from the hover information can be added to the
                         analysis via the regular input fields."),
                 h3("2. Select pairwise comparison of biological conditions"),
                 h3("3. Select minimal mean Difference cutoff."),
                 h3("4. Option to show all differential scores, also those not passing the mean Difference criterion set in (3)."),
                 h3("5. Output table with differential scores in the given comarison chosen in (2) with criteria set in (3) and (4)."),
                 p(),
                 
                 
                 
                 h2("4) View complex features"),
                 img(src='4.PNG', align = "left", width = "100%"),
                 h3("1. Select ID type based on which to select complex_id or feature_id"),
                 h3("2. Select one complex_id or unique_feature_id to display the chromatogram and feature graph"),
                 h3("3. Complex feature graph displaying for a selected complex the subunits' SEC chromatograms across conditions and detected complex signal ranges"),
                 tags$ul(tags$li("3.1 First feature detected for this complex, details are given in complex feature table below (5)"),
                         tags$li("3.2 Second feature detected for this complex, details are given in complex feature table below (5)")),
                 h3("4. Download current graph as .pdf"),
                 h3("5. Complex feature information table, giving details on the currently selected complex and features"),
                 tags$ul(tags$li("Note that all detected features can be displayed via the checkbox below field (2).")),
                 p(),
                 
                 
                 
                 h2("5) Complex differential intensity"),
                 img(src='5.PNG', align = "left", width = "100%"),
                 h3("1. Select pairwise comparison of biological conditions"),
                 h3("2. Select complex ID to be highlighted in the complex-level volcano plot (4)."),
                 h3("3. Set significance criteria for the complex-level volcano plot (4) and the hits shown in the output table (5)."),
                 h3("4. Complex-level differential volcano plot of the given comparison (set in (3))."),
                 tags$li("Outlier complexes identified from the hover information can be added to the
                         analysis via the respective input fields in Tab 4) and 5)."),
                 h3("5. Output table with complex-level difference scores of the current comparison (3) and significance criteria (4)"),
                 tags$ul(tags$li("Note that all detected features can be displayed via the checkbox below field (3).")),
                 p()
                 
        ),
        tabPanel('1) View protein traces',       
                 plotlyOutput("plot", height = 600),
                 downloadButton("downloadPlot", "Download graph as PDF"),
                 p(),
                 DT::dataTableOutput("anntable")
        ),
        tabPanel('2) Protein differential intensity',
                 plotlyOutput("pc_volcano"),
                 DT::dataTableOutput("pc_difftable")
        ),
        tabPanel('3) Protein differential assembly',
                 plotlyOutput("pc_assemblyScatter", height = 400),
                 DT::dataTableOutput("pc_assemblyTable")
        ),
        tabPanel('4) View complex features',
                 plotOutput("cf_plot", height=600),
                 downloadButton("downloadPlot_cF", "Download graph as PDF"),
                 p(),
                 # Not functional as no plot object is available for output
                 DT::dataTableOutput("cf_table")
        ),
        tabPanel('5) Complex differential intensity',
                 plotlyOutput("cc_volcano", height=400),
                 DT::dataTableOutput("cc_difftable")
        )
      )
    )
  )
))
