library(markdown)
library(shiny)
library(gridExtra)
library(dendextend)
library(heatmaply)
library(d3heatmap)
library(igraph)
library(conos)
library(ggplot2)
library(largeVis)
library(entropy)
library(ramify)

source("ShinyApp/helpers.R")

greedy.modularity.cut.result <- conos::greedy.modularity.cut(wt=con$clusters$walktrap$result,N=50,leaf.labels = mapping,minbreadth = 0.1)

leafContent <- greedy.modularity.cut.result$leafContent
hc <- greedy.modularity.cut.result$hc
dend <- as.dendrogram(hc)


#sequential numeration of elements from hclust merge matrix
nodes_labels <- c()
for (i in 1:nrow(hc$merge)){
  pair <- hc$merge[i,]
  nodes_labels <- c(nodes_labels,sequential.numeration(pair))
}
labels_mapping <- 1:length(nodes_labels)
names(labels_mapping) <- rev(nodes_labels)

#correspondence between cells and samples 
mapping <- mapping[rownames(leafContent)] 
sample_names <- unique(mapping)

#correspondence between elements from hclust merge matrix and precalculated breadth values 
breadth_mapping_names <- c(-(1:(nrow(hc$merge)+1)),1:(nrow(hc$merge)-1)) #excluding node corresponding to only one cluster 
breadth_mapping <- greedy.modularity.cut.result$breadth
breadth_mapping <- breadth_mapping[1:(length(breadth_mapping)-1)]
names(breadth_mapping) <- labels_mapping[as.character(breadth_mapping_names)]

dend.colored <- color_nodes(x=dend,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,breadth_mapping=breadth_mapping)

#coordinates for drawing labels 
xy <- get_nodes_xy(dend)
xy_ordered <- cbind(xy,get_nodes_attr(dend.colored$x,"ID"))
xy_ordered <- xy_ordered[order(xy_ordered[,1]),]
order_of_nodes <- xy_ordered[,3]

#data for plots under a slider 
modularities <- cumsum(greedy.modularity.cut.result$deltaM)
breadth <- calculate_breadth(dend=dend.colored$x,hc=hc,breadth_mapping=breadth_mapping,labels_mapping=labels_mapping)
resolution <- calculate_resolution(dend = dend.colored$x,hc=hc,mapping = mapping,sample_names = sample_names,labels_mapping=labels_mapping)

#calculate embedding once 
#con$embedGraph()

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


server <- function(input, output) { #plots under a slider 
  

  dataInput <- reactive({  
    
    current_cut <- labels_mapping[1:((input$N-1)*2)] 
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(input$N-1)]))]
    cut_order_match <- match(current_cut,order_of_nodes)
    cut_order_match <- cut_order_match[order(cut_order_match)]
    current_cut_ordered <- order_of_nodes[cut_order_match] #
    
    colors <- rainbow(input$N)
    names(colors) <-  current_cut_ordered 
    mask <- xy_ordered[,3] %in% current_cut
    dend.colored <- color_nodes(x=dend.colored$x,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,current_cut=current_cut,colors=colors,breadth_mapping=breadth_mapping)
    
    structure_vectors <- get_structure_vectors(dend.colored$x,current_cut_ordered)
    
    heights <- get_nodes_attr(dend.colored$x,"height")
    IDs <- get_nodes_attr(dend.colored$x,"ID")
    plot(dend.colored$x)
    
    height_mask <- (heights > 0) & (IDs %in% current_cut)
    
    if(any(height_mask)){height <- max(heights[height_mask])}
    else{height <- 0}
    
    dend.cut <- cut(dend.colored$x,h = height)
    
    dend.cut <- hang.dendrogram(dend.cut$upper, hang = 1)
    labels(dend.cut) <- current_cut_ordered
    
    memb <- get_memberhip(dend=dend.colored$x,current_cut=current_cut)
    
    return(list(dend.colored=dend.colored, mask=mask,structure_vectors=structure_vectors,dend.cut=dend.cut,memb=memb,colors=colors))
  })
  
  #regulation of number of digits at y ticks 
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  output$treePlot1 <- renderPlot({     #modularity/breadth/resolution pictures under slider 
    mods <- data.frame(cbind(modularities,2:(length(modularities)+1)))
    colnames(mods) <- c("Modularity","Number_of_clusters")
    b <-  data.frame(cbind((breadth-1/12),2:(length(breadth)+1)))
    colnames(b) <- c("Breadth","Number_of_clusters")
    r <- data.frame(cbind(round(resolution,1),2:(length(resolution)+1)))
    colnames(r) <- c("Resolution","Number_of_clusters")
    
    g12 <- subset(mods, Number_of_clusters == input$N)
    g22 <- subset(b, Number_of_clusters == input$N)
    g32 <- subset(r, Number_of_clusters == input$N)
    g1 <- ggplot(dat=mods,aes(x=Number_of_clusters, y=Modularity)) + 
      geom_point() + geom_point(data=g12, colour="red") + 
      scale_y_continuous(position = "right")  + 
      scale_x_continuous(expand = c(0.01, 0.01))  + 
      theme(axis.title.x = element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      geom_vline(xintercept=c(input$N,input$N), linetype="dotted") +
      ylab("Modularity")
    g2 <- ggplot(dat=b,aes(x=Number_of_clusters, y=Breadth)) + 
      geom_point() + geom_point(data=g22, colour="red") + 
      scale_y_continuous(position = "right") + 
      scale_x_continuous(expand = c(0.01, 0.01)) + 
      theme(axis.title.x = element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      geom_vline(xintercept=c(input$N,input$N), linetype="dotted") 
    
    g3 <- ggplot(dat=r,aes(x=Number_of_clusters, y=Resolution)) + 
      geom_point() + geom_point(data=g32, colour="red") + 
      scale_y_continuous(position = "right",labels=scaleFUN) +
      scale_x_continuous(expand = c(0.01, 0.01)) + 
      geom_vline(xintercept=c(input$N,input$N), linetype="dotted")
    glist <- list(g1,g2,g3)
    grid.arrange(grobs=glist,nrow=length(glist))
    
  },height = 250)
  output$treePlot2 <- renderPlot({ #draw colored tree 
    #color the tree and plot it with labels 
    dend.colored.local <- dataInput()$dend.colored
    labels(dend.colored.local$x) <- rep(NA,attr(dend.colored.local$x,"members"))
    mask <- dataInput()$mask
    plot(dend.colored.local$x)
    text(xy_ordered[mask,1]+0.5, xy_ordered[mask,2]+0.5, labels=xy_ordered[mask,3], col="black",cex = 0.8)
    
  },width = 1000, height = 800)
  
  output$treePlot3 <- renderD3heatmap({ #
    struct_similarity <- dataInput()$structure_vectors
    Rowv <- dataInput()$dend.cut
    d3heatmap(t(struct_similarity),colors = Reds(10),Rowv=Rowv, dendrogram = "column" )})#,
  
  output$treePlot4 <- renderD3heatmap({
    struct_similarity <- dataInput()$structure_vectors
    Rowv <- dataInput()$dend.cut
    dists <- dist(t(struct_similarity))
    d3heatmap(dists,colors = rev(Reds((max(dists))%/%0.1)),Rowv=Rowv,Colv=Rowv)}) #,
  
  output$treePlot5 <- renderPlot({
    memb <- dataInput()$memb
    colors <- dataInput()$colors  
    con$plotPanel(groups=memb, adjust.func = function(gg) gg + scale_color_manual(values = setNames(colors, names(colors))))
  },width = 1000, height = 800)
  
  output$treePlot6 <- renderPlot({
    memb <- dataInput()$memb
    colors <- dataInput()$colors 
    con$plotGraph(groups = memb) + scale_color_manual(values = setNames(colors, names(colors)))
  }, height = 320)
  
  observeEvent(input$do, {
    memb <- dataInput()$memb
    gready.cut.groups <- paste("gready.cut.groups.",".clusters",sep=as.character(input$N))
    names <- names(memb)
    con$clusters$walktrap[[gready.cut.groups]] <- as.character(memb)
    names(con$clusters$walktrap[[gready.cut.groups]]) <- names
  })
  
}

shinyApp(ui, server)




