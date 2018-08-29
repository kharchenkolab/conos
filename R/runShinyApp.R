#' @import markdown
#' @import shiny
#' @import gridExtra
#' @import dendextend
#' @import heatmaply
#' @import d3heatmap
#' @import largeVis
#' @import entropy
#' @import ramify

#' Perform ordering of merge matrix values by rows 
sequential.numeration <- function(pair){
  if ((pair[1] < 0) & (pair[2] < 0)){
    return(c(pair[1],pair[2]))
  }
  else if ((pair[1] > 0) & (pair[2] > 0)){
    return(c(min(pair),max(pair)))
  }
  else{
    return(c(max(pair),min(pair)))
  }
}

#' Return row for value in merge matrix
which.row <- function(hc,label){
  if (label %in% hc$merge[,1]){
    return(which(hc$merge[,1]==label))
  }
  else if (label %in% hc$merge[,2]){
    return(which(hc$merge[,2]==label))
  }
}

#' Calculate breadth based on precalculated values 
calculate_breadth <- function(dend=NULL,hc=NULL,breadth_mapping=NULL,labels_mapping=NULL){
  breadth <- c()
  for (i in 2:attr(dend,"members")){
    current_cut <- labels_mapping[1:((i-1)*2)] # 4 = N
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(i-1)]))]
    breadth <- c(breadth,mean(breadth_mapping[as.character(current_cut)]))
  }
  return(breadth)
}

#' Get membership for current cut 
get_memberhip <- function(dend=NULL,current_cut=NULL){
  IDs <- get_nodes_attr(dend,"ID")
  memberships <- get_nodes_attr(dend,"membership")
  membership <- memberships[,IDs %in% current_cut]
  if (ncol(membership)>1){
    membership <- rowSums(membership)
  }
  return(membership)
}

#' Get composition vectors for each cluster - proportion of cells from each sample
get_structure_vectors <- function(dend=NULL,current_cut=NULL){
  IDs <- get_nodes_attr(dend,"ID")
  structure_vectors <- get_nodes_attr(dend,"structure_vector")
  
  cut_IDs_match <- match(current_cut,IDs)
  cut_IDs_match <- cut_IDs_match[order(cut_IDs_match)]
  
  structure_vectors <- structure_vectors[,cut_IDs_match]
  colnames(structure_vectors) <- current_cut
  return(structure_vectors)
}

#' Calculate resolution for each cut
calculate_resolution <- function(dend=NULL,hc=NULL,mapping=NULL,sample_names=NULL,labels_mapping=NULL){
  resolution <- c()
  for (i in 2:attr(dend,"members")){
    current_cut <- labels_mapping[1:((i-1)*2)] 
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(i-1)]))]
    membership <- get_memberhip(dend,current_cut)
    local_resolution <- sapply(sample_names,function(name) length(unique(membership[mapping==name])))
    resolution <- c(resolution,mean(local_resolution))
  }
  return(resolution )
}

#' Color node and adjacent branches based on label of cluster, number of cells and breadth
color_node <- function(x=NULL,labels_mapping=NULL,hc=NULL,leafContent=NULL,mapping=NULL,sample_names=NULL,current_cut=NULL,colors=NULL,parent_color=NULL,breadth_mapping=NULL){
  if (is.null(attr(x, "merge_ID"))){
    if (attr(x[[1]],"height") > attr(x[[2]],"height")) {
      attr(x[[1]],"merge_ID") = max(hc$merge[nrow(hc$merge),])
      attr(x[[2]],"merge_ID") = min(hc$merge[nrow(hc$merge),])
      attr(x[[1]],"ID") = as.numeric(labels_mapping[as.character(max(hc$merge[nrow(hc$merge),]))])
      attr(x[[2]],"ID") = as.numeric(labels_mapping[as.character(min(hc$merge[nrow(hc$merge),]))])
    }
    else{
      attr(x[[1]],"merge_ID") = min(hc$merge[nrow(hc$merge),])
      attr(x[[2]],"merge_ID") = max(hc$merge[nrow(hc$merge),])
      attr(x[[1]],"ID") = as.numeric(labels_mapping[as.character(min(hc$merge[nrow(hc$merge),]))])
      attr(x[[2]],"ID") = as.numeric(labels_mapping[as.character(max(hc$merge[nrow(hc$merge),]))])
    }
    attr(x,"ID") <- 0
  }
  else{
    if (!is.leaf(x)){
      if (attr(x[[1]],"height") > attr(x[[2]],"height")) {
        attr(x[[1]],"merge_ID") = max(hc$merge[attr(x,"merge_ID"),])
        attr(x[[2]],"merge_ID") = min(hc$merge[attr(x,"merge_ID"),])
        attr(x[[1]],"ID") = as.numeric(labels_mapping[as.character(max(hc$merge[attr(x,"merge_ID"),]))])
        attr(x[[2]],"ID") = as.numeric(labels_mapping[as.character(min(hc$merge[attr(x,"merge_ID"),]))])
      }
      else if (attr(x[[1]],"height") < attr(x[[2]],"height")) {
        attr(x[[1]],"merge_ID") = min(hc$merge[attr(x,"merge_ID"),])
        attr(x[[2]],"merge_ID") = max(hc$merge[attr(x,"merge_ID"),])
        attr(x[[1]],"ID") = as.numeric(labels_mapping[as.character(min(hc$merge[attr(x,"merge_ID"),]))])
        attr(x[[2]],"ID") = as.numeric(labels_mapping[as.character(max(hc$merge[attr(x,"merge_ID"),]))])
      }
      else{
        attr(x[[1]],"ID") = as.numeric(labels_mapping[as.character(-labels(x[[1]]))])
        attr(x[[2]],"ID") = as.numeric(labels_mapping[as.character(-labels(x[[2]]))])
        attr(x[[1]],"merge_ID") = 0
        attr(x[[2]],"merge_ID") = 0
      }
    }
    
  }
  if (is.null(mapping)){ #just index and get memberships
    labs <- labels(x)
    if (length(labs)>1){
      memb <- rowSums(leafContent[,labs])
    }
    else{
      memb <- leafContent[,labs]
    }
    attr(x, "membership") <- memb*attr(x, "ID")
  }
  else{
    if (is.null(attr(x, "entropy"))){
      labs <- labels(x)
      if (length(labs)>1){
        memb <- rowSums(leafContent[,labs])
      }
      else{
        memb <- leafContent[,labs]
      }
      attr(x, "membership") <- memb*attr(x, "ID")
      structure_vector <- rep(0,length(sample_names))
      names(structure_vector) <- sample_names
      number_of_cells <- sum(memb)
      attr(x, "number_of_cells") <- number_of_cells
      
      structure_tab <- table(mapping[as.logical(memb)])
      for (name in names(structure_tab)){
        structure_vector[name] <- structure_tab[name]
      }
      structure_vector <- structure_vector/sum(structure_vector)
      attr(x, "structure_vector") <- structure_vector
      #entropy_measure <- entropy(structure_vector,unit = "log2")/log(length(sample_names),base = 2) + 1/length(sample_names)
      entropy_measure <- as.numeric(breadth_mapping[as.character(attr(x,"ID"))]) + 1/length(sample_names)
      attr(x, "entropy") <- entropy_measure  #/length(sample_names)
    }
    else{
      number_of_cells <- attr(x, "number_of_cells")
      entropy_measure <- attr(x, "entropy")
    }
    if (attr(x,"ID") %in% current_cut){
      x <- x %>% 
        assign_values_to_nodes_nodePar(value = 19, nodePar = "pch") %>% 
        assign_values_to_nodes_nodePar(value = number_of_cells/1000, nodePar = "cex")  %>% 
        set("branches_lty", 1) %>%
        set("branches_lwd", entropy_measure*5)# %>%
      parent_color <- as.character(colors[as.character(attr(x,"ID"))])
      x <- assign_values_to_nodes_nodePar(dend = x, value = rep(parent_color,attr(x,"members")), nodePar = "col")
      x <- assign_values_to_branches_edgePar(dend = x, value = rep(parent_color,attr(x,"members")), edgePar = "col")
    }
    else{
      x <- x %>% 
        assign_values_to_nodes_nodePar(value = NA, nodePar = "pch") %>% 
        assign_values_to_nodes_nodePar(value = NA, nodePar = "cex")  %>% 
        assign_values_to_nodes_nodePar(value = NA, nodePar = "col")  %>%
        set("branches_lty", 1) %>%
        set("branches_lwd", entropy_measure*5) #%>%
      
      if (!is.leaf(x) & !(is.null(parent_color))){
        x <- assign_values_to_branches_edgePar(dend = x, value = rep(parent_color,attr(x,"members")), edgePar = "col")
      }
      if (is.null(parent_color)){
        x <- x %>% set("branches_col", "grey80")
      }
    }
  }
  return(list(x=x,parent_color=parent_color))
}

#' Apply color_node function for the entire tree
color_nodes <- function (x=NULL,labels_mapping=NULL,hc=NULL,leafContent=NULL,mapping=NULL,sample_names=NULL,current_cut=NULL,colors=NULL,parent_color=NULL,breadth_mapping=NULL){
  if (!is.leaf(x)){
    res <- color_node(x=x,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,current_cut=current_cut,colors=colors,parent_color=parent_color,breadth_mapping=breadth_mapping)
    x <- res$x
    parent_color <- res$parent_color
    for (j in 1:2){
      res <- color_nodes(x=x[[j]],labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,current_cut=current_cut,colors=colors,parent_color=parent_color,breadth_mapping=breadth_mapping)
      x[[j]] <- res$x
    }
  }
  else{
    res <- color_node(x=x,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,current_cut=current_cut,colors=colors,parent_color=parent_color,breadth_mapping=breadth_mapping)
    x <- res$x
    parent_color <- res$parent_color
  }
  return(list(x=x,parent_color=parent_color))
}

#get list of children nodes of a node
get.list.of.descendants <- function(hc=NULL,parent.ID=NULL,list.of.descendants=list()){
  children.ID <- hc$merge[parent.ID,]
  if (length(list.of.descendants)>0){
    list.of.descendants[length(list.of.descendants)+1] <- children.ID[1]
    list.of.descendants[length(list.of.descendants)+1] <- children.ID[2]
  }
  else{
    list.of.descendants[[1]] <- children.ID[1]
    list.of.descendants[[2]] <- children.ID[2]
  }
  
  if (children.ID[1] > 0){
    list.of.descendants <- get.list.of.descendants(hc, children.ID[1], list.of.descendants)
  }
  
  if (children.ID[2] > 0){
    list.of.descendants <- get.list.of.descendants(hc, children.ID[2], list.of.descendants)
  }
  return(list.of.descendants)
}

##' get memerships vector from greedy.modularity.cut function results without running shiny app 
##'
##' @param no_clusters number of clusters less or equal to N + 1 number of splits in greedy.modularity.cut results
##' @param greedy.modularity.cut.result greedy.modularity.cut function results 
##' @export
get.greedy.cut.groups <- function(no_clusters=NULL,greedy.modularity.cut.result=NULL){
  hc <- greedy.modularity.cut.result$hc
  dend <- as.dendrogram(hc)
  leafContent <- greedy.modularity.cut.result$leafContent
  max_no_clusters <- attr(dend,"members")
  if (no_clusters <= max_no_clusters){
    nodes_labels <- c()
    
    for (i in 1:nrow(hc$merge)){
      pair <- hc$merge[i,]
      nodes_labels <- c(nodes_labels,sequential.numeration(pair))
    }
    labels_mapping <- 1:length(nodes_labels)
    names(labels_mapping) <- rev(nodes_labels)
    
    dend.indexed <- color_nodes(x=dend,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent)
    
    current_cut <- labels_mapping[1:((no_clusters-1)*2)] 
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(no_clusters-1)]))]
    
    memb <- get_memberhip(dend=dend.indexed$x,current_cut=current_cut)
    return(memb)
  }
  else{
    stop("no_clusters exceeds a maximum number of clusters in greedy.modularity.cut results (N + 1)")
  }
}

##' deploys Shiny application to visualize waltktrap tree and select cut level 
##'
##' @param con conos object with walktrap results
##' @param N number of top greedy splits to take
##' @param leaf.labels leaf sample label factor, for breadth calculations - must be a named factor containing all wt$names, or if wt$names is null, a factor listing cells in the same order as wt leafs
##' @param minsize minimum size of the branch (in number of leafs) 
##' @param minbreadth minimum allowed breadth of a branch (measured as normalized entropy)
##' @export
runShinyApp <- function(con=NULL, N=NULL, leaf.labels=NULL, minsize=0, minbreadth=0) {
  #.GlobalEnv$con <- con
  greedy.modularity.cut.result <- conos::greedy.modularity.cut(wt=con$clusters$walktrap$result,N=N,leaf.labels=leaf.labels,minsize=minsize,minbreadth=minbreadth)
  
  leafContent <- greedy.modularity.cut.result$leafContent
  hc <- greedy.modularity.cut.result$hc
  dend <- as.dendrogram(hc)
  
  #sequential numeration of elements from hclust merge matrix
  nodes_labels <- c()
  for (i in 1:nrow(hc$merge)){
    pair <- hc$merge[i,]
    nodes_labels <- c(nodes_labels,sequential.numeration(pair))
  }
  
  #labels_mapping <- 0:length(nodes_labels)
  #names(labels_mapping) <- c(nrow(hc$merge),rev(nodes_labels))
  labels_mapping <- 1:length(nodes_labels)
  names(labels_mapping) <- rev(nodes_labels)
  
  #correspondence between cells and samples 
  mapping <- leaf.labels[rownames(leafContent)] 
  sample_names <- unique(mapping)
  
  #correspondence between elements from hclust merge matrix and precalculated breadth values 
  breadth_mapping_names <- c(-(1:(nrow(hc$merge)+1)),1:(nrow(hc$merge)-1)) #excluding node corresponding to only one cluster 
  breadth_mapping <- greedy.modularity.cut.result$breadth
  breadth_mapping <- breadth_mapping[1:(length(breadth_mapping)-1)]
  names(breadth_mapping) <- labels_mapping[as.character(breadth_mapping_names)]
  
  dend.colored.source <- color_nodes(x=dend,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,breadth_mapping=breadth_mapping)
  
  #coordinates for drawing labels 
  xy <- get_nodes_xy(dend)
  xy_ordered <- cbind(xy,get_nodes_attr(dend.colored.source$x,"ID"))
  xy_ordered <- xy_ordered[order(xy_ordered[,1]),]
  order_of_nodes <- xy_ordered[,3]
  
  #data for plots under a slider 
  modularities <- cumsum(greedy.modularity.cut.result$deltaM)
  breadth <- calculate_breadth(dend=dend.colored.source$x,hc=hc,breadth_mapping=breadth_mapping,labels_mapping=labels_mapping)
  resolution <- calculate_resolution(dend = dend.colored.source$x,hc=hc,mapping = mapping,sample_names = sample_names,labels_mapping=labels_mapping)
  
  #calculate embedding once 
  #con$embedGraph()
  
  colors_source <- sample(rainbow(length(labels_mapping)))
  names(colors_source) <- sample(labels_mapping)
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
                        actionButton("do", "Save membership"),
                        #actionButton("plot_graph", "Plot the graph")),
                        checkboxInput("plot_graph", "Plot graph", value = TRUE)),
                      div(style = "height:250;"))),
             fluidRow(
               column(12,
                      tabPanel("Plot",plotOutput(outputId = "treePlot6")),
                      div(style = "height:500;")))),
      column(8,                       
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Tree", plotOutput(outputId = "treePlot2",click = "plot_click",dblclick="plot_db_click")),
                           tabPanel("Composition similarity plot 1", d3heatmapOutput(outputId = "treePlot3", width = "1000", height = "800")),
                           tabPanel("Composition similarity plot 2", d3heatmapOutput(outputId = "treePlot4", width = "1000", height = "800")),
                           tabPanel("tSNE plots", plotOutput(outputId = "treePlot5")))),div(style = "height:800;"))))
  
  
  server <- function(input, output) { #plots under a slider 
    click_value <- reactiveVal(NULL)
    observeEvent(input$plot_click,{click_value(input$plot_click)})
    observeEvent(input$N,{click_value(NULL)})
    
    db_click_value <- reactiveVal(NULL)
    observeEvent(input$plot_db_click,{db_click_value(input$plot_db_click)})
    observeEvent(input$N,{db_click_value(NULL)})
    
    observeEvent(input$plot_click,{db_click_value(NULL)})
    observeEvent(input$plot_db_click,{click_value(NULL)})
    
    current_cut_saved <- reactiveVal(NULL)
    observeEvent(input$N,{current_cut_saved(NULL)})
    changed <- reactiveVal(FALSE)
    observeEvent(input$N,{changed(TRUE)})
    
    dataInput <- reactive({  
      if (changed()){
        current_cut <- NULL
        current_cut <- labels_mapping[1:((input$N-1)*2)] 
        current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(input$N-1)]))]
        isolate(current_cut_saved(current_cut))
        changed(FALSE)
      }
      else{
        current_cut <- isolate(current_cut_saved())
      }
      
      if (!is.null(click_value())){
        euc_dist <- (xy_ordered[,1]-click_value()$x)^2 + (xy_ordered[,2]-click_value()$y)^2 
        ID <- xy_ordered[which(euc_dist==min(euc_dist)),3]
        if (ID > 0){
          merge_mtrx_ID <- names(labels_mapping[labels_mapping==ID])
          num_ID <- as.numeric(merge_mtrx_ID)
          if (num_ID > 0){ 
            if (!(any(labels_mapping[as.character(unlist(get.list.of.descendants(hc,num_ID)))] %in% current_cut)) & (ID %in% current_cut)){
              current_cut <- current_cut[current_cut!=ID]
              child1 <- hc$merge[num_ID,1]
              child2 <- hc$merge[num_ID,2]
              current_cut_names <- c(names(current_cut),c(child1,child2))
              current_cut <- c(current_cut,c(labels_mapping[as.character(child1)],labels_mapping[as.character(child2)]))
              names(current_cut) <- current_cut_names
              isolate(current_cut_saved(current_cut))
            }
          }
        }    
      }
      
      if (!is.null(db_click_value())){
        euc_dist <- (xy_ordered[,1]-db_click_value()$x)^2 + (xy_ordered[,2]-db_click_value()$y)^2 
        ID <- xy_ordered[which(euc_dist==min(euc_dist)),3]
        if (ID > 0){
          merge_mtrx_ID <- names(labels_mapping[labels_mapping==ID])
          num_ID <- as.numeric(merge_mtrx_ID)
          if (num_ID > 0){ 
            array_of_parents <- as.numeric(names(labels_mapping[labels_mapping %in% current_cut]))
            array_of_parents <- array_of_parents[array_of_parents > 0]
            array_of_children <- unlist(sapply(array_of_parents,function(num_ID){labels_mapping[as.character(unlist(get.list.of.descendants(hc,num_ID)))]}))
            array_of_children <- array_of_children[as.numeric(names(array_of_children)) > 0]
            if (!(ID %in% array_of_children)){
              current_cut <- current_cut[!(current_cut %in% labels_mapping[as.character(unlist(get.list.of.descendants(hc,num_ID)))])]
              current_cut_names <- c(names(current_cut),as.character(num_ID))
              current_cut <- c(current_cut,ID)
              names(current_cut) <- current_cut_names
              isolate(current_cut_saved(current_cut))
            }
          }
        }
      }
      
      cut_order_match <- match(current_cut,order_of_nodes)
      cut_order_match <- cut_order_match[order(cut_order_match)]
      current_cut_ordered <- order_of_nodes[cut_order_match]
      
      #colors <- rainbow(length(current_cut))
      #names(colors) <-  sample(current_cut)
      colors <- colors_source[names(colors_source) %in% current_cut]
      mask <- xy_ordered[,3] %in% current_cut
      dend.colored <- color_nodes(x=dend.colored.source$x,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,current_cut=current_cut,colors=colors,breadth_mapping=breadth_mapping)
      
      structure_vectors <- get_structure_vectors(dend.colored$x,current_cut_ordered)
      
      heights <- get_nodes_attr(dend.colored$x,"height")
      IDs <- get_nodes_attr(dend.colored$x,"ID")
      #plot(dend.colored$x)
      
      height_mask <- (heights > 0) & (IDs %in% current_cut)
      
      if(any(height_mask)){height <- max(heights[height_mask])}
      else{height <- 0}
      
      dend.cut <- cut(dend.colored$x,h = height)
      
      dend.cut <- hang.dendrogram(dend.cut$upper, hang = 1)
      labels(dend.cut) <- current_cut_ordered
      
      memb <- get_memberhip(dend=dend.colored$x,current_cut=current_cut)
      return(list(dend.colored=dend.colored, mask=mask,structure_vectors=structure_vectors,dend.cut=dend.cut,memb=memb,colors=colors,current_cut=current_cut))
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
      #plot(value())
      dend.colored.local <- dataInput()$dend.colored
      labels(dend.colored.local$x) <- rep(NA,attr(dend.colored.local$x,"members"))
      mask <- dataInput()$mask
      
      
      plot(dend.colored.local$x)
      text(xy_ordered[mask,1]+0.5, xy_ordered[mask,2]+0.5, labels=xy_ordered[mask,3], col="black",cex = 0.8)
      
    },width = 1000, height = 800)
    
    output$treePlot3 <- renderD3heatmap({ #
      struct_similarity <- dataInput()$structure_vectors
      if (is.null(click_value()) & is.null(db_click_value())){
        Rowv <- dataInput()$dend.cut
        d3heatmap(t(struct_similarity),colors = Reds(10),Rowv=Rowv, dendrogram = "column" )}
      else{
        d3heatmap(t(struct_similarity),colors = Reds(10), dendrogram = "column" )}
    })
    
    output$treePlot4 <- renderD3heatmap({
      struct_similarity <- dataInput()$structure_vectors
      dists <- dist(t(struct_similarity))
      if (is.null(click_value()) & is.null(db_click_value())){
        Rowv <- dataInput()$dend.cut
        d3heatmap(dists,colors = rev(Reds((max(dists))%/%0.1)),Rowv=Rowv,Colv=Rowv)}
      else{
        d3heatmap(dists,colors = rev(Reds((max(dists))%/%0.1)),dendrogram = "none")
      }
    })
    
    output$treePlot5 <- renderPlot({
      #plot(input$plot_graph)
      #if (plot_graph()){
      memb <- dataInput()$memb
      colors <- dataInput()$colors  
      con$plotPanel(groups=memb, adjust.func = function(gg) gg + scale_color_manual(values = colors))
      
    },width = 1000, height = 800)
    
    output$treePlot6 <- renderPlot({
      if (input$plot_graph){
        #if (plot_graph()){
        memb <- dataInput()$memb
        colors <- dataInput()$colors
        con$plotGraph(groups = memb) + scale_color_manual(values = colors)
        #}
      }
    }, height = 320)
    
    observeEvent(input$do, {
      memb <- dataInput()$memb
      
      if (is.null(click_value()) & is.null(db_click_value())){
        gready.cut.groups <- paste("gready.cut.groups.",".clusters",sep=as.character(input$N))
        names <- names(memb)
      }
      else{
        current_cut <- dataInput()$current_cut
        gready.cut.groups <- paste(paste("gready.cut.groups.",".clusters",sep=as.character(length(current_cut))),paste(current_cut,collapse="_"),sep="_")
        names <- names(memb)
      }
      con$clusters$walktrap[[gready.cut.groups]] <- as.character(memb)
      names(con$clusters$walktrap[[gready.cut.groups]]) <- names
    })
    
  }
  
  shinyApp(ui, server)
}



