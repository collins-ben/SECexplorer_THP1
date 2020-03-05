####################################################################################
# SECexplorer
#####################################################################################
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
  headerPanel("SECexplorer", windowTitle = "SECexplorer"),

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
                       numericInput("pc_assemblyScatter_meanDiffcutoff", label = "mean Difference cutoff", value = 0.3),
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
                                      selected = "Differentiated vs. undifferentiated",
                                      multiple = FALSE),
                       p("Identify regulated complexes from plot hover info.\n
                         To analyze traces, select complex of interest via the field below:"),
                       selectizeInput(inputId = "complexid",
                                      label = "Search and select complex of interest",
                                      choices = complexes_diff_ids,
                                      multiple = FALSE, 
                                      # options = list(maxOptions = 6000),
                                      selected = "127_corum_corum"),
                       numericInput("cc_volcano_pvalcutoff", label = "pBHadj cutoff:", value = 0.05, max = 1),
                       numericInput("cc_volcano_fccutoff", label = "pBHadj cutoff:", value = 2),
                       checkboxInput("cc_difftable_show_all", "Show all (i.e. also non-significant)\n entries in table", value = FALSE)
                       
      ),
      br(),
      br(),
      br(),
      p("Prototype viewer app, MH 2020-02-22")
    ),

    # Instruction and analysis panels
    mainPanel(
      tabsetPanel(
        id = 'tabid',
        tabPanel('Welcome & Instructions',
                 h1("Welcome to SECexplorer"),
                 p("SECexplorer offers a dynamic interface for the investigator-driven exploration of mass spectrometric SEC-SWATH-MS datasets describing proteome assembly state dynamics as a function of different perturbed states of the biological system."),
                 tags$ol(
                   tags$ul("1) View protein traces: Allows interactive viewing of protein SEC fractionation profiles across conditions."),
                   tags$ul("2) Protein differential intensity: View differential scores on protein level across conditions."),
                   tags$ul("3) Protein differential assembly: View proteins differential assembled mass fraction across conditions."),
                   tags$ul("4) View complex features: View detected complex signals across conditions."),
                   tags$ul("5) Complex differential intensity: View differential scores on complex level across conditions.")
                 ),
                 h1("Usage instructions:"),
                 h2("(i) Protein profile viewer"),
                 img(src='1_Viewer.png', align = "left", width = "100%"),
                 h3("1. Selection of ID type"),
                 tags$ul(
                   tags$li("Entry_name or Gene_Names are most informative and intuitive to search."),
                   tags$li("Searching for specific gene/protein identifiers is then possible in field 2.")
                   ),
                 
                 h3("2. Search and selection of multiple proteins"),
                 tags$ul(
                   tags$li("Searching is achieved by deleting the current entry with backspace and starting to 
                 type. All identifiers of the type selected in (1) will be searched for the string 
                     entered, on-the-fly, with potential results showing up below the field, selectable 
                     by <Enter> or by a left mouse click."),
                   tags$li("NOTE: Individual proteins can be removed again by <Backspace>.")
                 ),
                 h3("3. Chromatogram view across conditions "),
                 tags$ul(
                   tags$li("By default, a split graph shows the protein level abundance profiles over the 
                            chromatographic fractions in either condition, interphase and mitosis."),
                   tags$li("The chromatograms are displayed as mean of the three replicates, with shaded areas indicating the standard error of the mean (SEM)"),
                   tags$li("Options for chromatogram display can be selected in (4)")
                 ),
                 h3("4. Options for chromatogram display"),
                 tags$ul(
                   tags$li("Selection whether the plot should be split by condition, or if conditions shall be encoded as line type."),
                   tags$li("Option for LOG10-transformation of the intensity axis to spot low-abundant protein pools"),
                   tags$li("Selection whether monomer expected fraction markers shall be displayed"),
                   tags$li("Option whether (and which type of) error areas shall be indicated. Works best on split plot.")
                   ),
                 h3("Note that the plots are interactive (thanks, plotly!)"),
                 p(),
                 
                 h2("(ii) Search for co-eluting proteins"),
                 img(src='HowtoSearchInput.png', align = "left", width = "100%"),
                 h3("1. Experimental condition"),
                 tags$ul(tags$li("Choose the experimental condition where the search will be performed")),
                 h3("2. Base Protein"),
                 tags$ul(tags$li("Choose the protein Trace that will be searched")),
                 h3("3. Search area"),
                 tags$ul(
                   tags$li("Define in which range co-eluting proteins should be searched for"),
                   tags$li("If no selection is made all fractions will be used")
                   ),
                 h3("4. Perform the search"),
                 p(),
                 h2("(ii) Search for co-eluting proteins: Search Results"),
                 img(src='HowtoSearchResult.png', align = "left", width = "100%"),
                 h3("1. Correlation cutoff"),
                 tags$ul(
                   tags$li("Choose a cutoff for the local and global correlation of proteins to be displayed"),
                   tags$li("Global correllation: The minimum correllation across all selected SEC fractions"),
                   tags$li("Local correllation: The minimum correllation within the defined search area")
                 ),
                 h3("2. Result table"),
                 tags$ul(tags$li("A table of all proteins meeting the search criteria")),
                 h3("3. Reset"),
                 tags$ul(tags$li("Return to the search input window")),
                 p(),
                 h2("(iii) Browse differential association map"),
                 img(src='3_BrowseScoreMap.png', align = "left", width = "100%"),
                 h3("1. Select proteins based on scores"),
                 tags$ul(
                   tags$li("Each point represents one protein feature's differential scores between the conditions. Click the point or draw a box around the scores to select a set of proteins."),
                   tags$li("Proteins are added to those already selected via the viewer. Visit this page to see the differential scores of a selected protein set of interest.")),
                 h3("2. Conditional SEC profiles of the selected proteins. For more visualization options, visit the viewer tab (i)"),
                 p(),
                 h2("(iv) Query StringDB partners"),
                 img(src='4_QueryStringPartners.png', align = "left", width = "100%"),
                 h3("1. Selected proteins"),
                 h3("2. Select one protein as basis for the query to StringDB"),
                 h3("3. Select parameters for the query to StringDB"),
                 h3("4. Add the retrieved partner proteins to the current analysis"),
                 h3("5. StringDB network view of the current query"),
                 h3("6. Preview of the SEC-SWATH-MS chromatograms of the proteins of the current query (MS-detectable subset)")
                 
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
