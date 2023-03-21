library(shiny)
library(plotly)
shinyUI(pageWithSidebar(
  headerPanel("Environmental Enrichment TgDyrk1A analysis"),
  sidebarPanel(
    # checkboxGroupInput("variablegroups", "Choose groups:","",selected="variablegroups"),
    # actionLink("selectall","Select All"),
    width=1
    # #selectInput("variable",label = h5("Variables for PCA"),"",multiple=T)
  ),
  # sidebarPanel(
  #   checkboxGroupInput("variable3", "Choose grouping for colors:",""),position=c('right')
  #   #selectInput("variable",label = h5("Variables for PCA"),"",multiple=T)
  # ),

  mainPanel(
    
    # Output: Tabset w/ plot, summary, and table ----
    tabsetPanel(type = "tabs",
                tabPanel("General",
                         selectInput("variablecomp2", "Choose Comparison","",selected="variablecomp2",width='100%'),
                         selectInput("variablemet", "Choose Metric","",selected="variablemet",width='100%'),
                         plotOutput("metric",height = "600px", width = "100%"),
                         downloadLink("downloadPlot", "Download Plot"),
                         # fluidRow(
                         #   column(width = 12, 
                         #          h3('Statistics'),  
                         #          # tableOutput("table")
                         #   ))
                         ),
                tabPanel("Layer Width",
                         selectInput("variablecomp", "Choose Comparison","",selected="variablecomp",width='100%'),
                         selectInput("variablesubl", "Choose Sublayer","",selected="variablesubl",width='100%'),
                         plotOutput("width",height = "600px", width = "100%"),
                         downloadLink("downloadPlot4", "Download Plot"),
                         # fluidRow(
                         #   column(width = 12, 
                         #          h3('Statistics'),  
                         #          # tableOutput("tablewidth")
                         #   ))
                         ),
                tabPanel("Reconstructions",
                         selectInput("variablecomp3", "Choose Comparison","",selected="variablecomp2",width='100%'),
                         selectInput("variablerecmet", "Choose Metric","",selected="variablemet",width='100%'),
                         plotOutput("recmet",height = "600px", width = "100%"),
                         downloadLink("downloadPlot2", "Download Plot"),
                         # fluidRow(
                         #   column(width = 12, 
                         #          h3('Statistics'),  
                         #          # tableOutput("tablerec")
                           # )
                         # ),
                         plotOutput("projections",height = "600px", width = "100%")
                ),
                tabPanel("Sholl Analysis",  
                         selectInput("variableshollgroup", "Choose comparison","",selected="variableshollgroup",width='100%'),
                         plotOutput("sholl",height = "600px", width = "100%"),
                         # checkboxGroupInput("variablegroups", "Choose groups:","",selected="variablegroups"),
                         # actionLink("selectall","Select All"),
                         downloadLink("downloadPlotsholl", "Download Plot")),
                tabPanel("Input Output frequency", 
                         selectInput("variableinoutgroup", "Choose comparison","",selected="variableinoutgroup",width='100%'),
                         selectInput("variableinoutregion", "Choose input region","",selected="variableinoutregion",width='100%'),
                         plotOutput("freqinout",height = "900px", width = "100%"),
                         downloadLink("downloadPlot5", "Download Plot")),
                tabPanel("LTP electrophysiology",
                         selectInput("variablecomp4", "Choose Comparison","",selected="variablecomp4",width='100%'),
                         plotOutput("slopeltp",height = "500px", width = "100%"),
                         downloadLink("downloadPlot3", "Download Plot"))
    )
  ))
  # mainPanel(
  #   #tableOutput("inputfile"),
  #   plotlyOutput("PCA",height = "800px", width = "100%")
  # )
)