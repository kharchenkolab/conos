ui <- fluidPage(
  fluidRow(
    column(4,
           fluidRow(
             column(12, 
                    sidebarPanel(
                      sliderInput(inputId = "N",
                                  label = "Number of clusters:",
                                  step = 1,
                                  min = 2,
                                  max = attr(dend,"members"),
                                  value = which.max(modularities)+1,
                                  width = "94%"), width=12, 
                      tabPanel("Plot",plotOutput(outputId = "treePlot1", height = "250"),br()),
                      actionButton("do", "Save membership")),
                    div(style = "height:250;"))),
           fluidRow(
             column(12,
                    tabPanel("Plot",plotOutput(outputId = "treePlot6")),
                    div(style = "height:500;")))),
    column(8,                       
           mainPanel(
             tabsetPanel(type = "tabs",
                         tabPanel("Tree", plotOutput(outputId = "treePlot2")),
                         tabPanel("Composition similarity plot 1", d3heatmapOutput(outputId = "treePlot3", width = "1000", height = "800")),
                         tabPanel("Composition similarity plot 2", d3heatmapOutput(outputId = "treePlot4", width = "1000", height = "800")),
                         tabPanel("tSNE plots", plotOutput(outputId = "treePlot5")))),div(style = "height:800;"))))