#' @import shiny
#' @import ggplot2
#' @import dendextend
NULL

#' Perform ordering of merge matrix values by rows
sequentialNumeration <- function(pair){
  if ((pair[1] < 0) & (pair[2] < 0)){
    return(c(pair[1], pair[2]))
  }

  if ((pair[1] > 0) & (pair[2] > 0)){
    return(c(min(pair), max(pair)))
  }

  return(c(max(pair), min(pair)))
}

#' Calculate breadth based on precalculated values
calculateBreadth <- function(dend=NULL,hc=NULL,breadth_mapping=NULL,labels_mapping=NULL){
  breadth <- c()
  for (i in 2:attr(dend,"members")){
    current_cut <- labels_mapping[1:((i-1)*2)] # 4 = N
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(i-1)]))]
    breadth <- c(breadth, mean(breadth_mapping[as.character(current_cut)]))
  }

  return(breadth)
}

#' Get membership for current cut
getMemberhip <- function(dend=NULL,current_cut=NULL){
  ids <- get_nodes_attr(dend,"ID")
  memberships <- get_nodes_attr(dend,"membership")
  membership <- memberships[,ids %in% current_cut]
  if (ncol(membership)>1){
    membership <- rowSums(membership)
  }
  return(membership)
}

#' Get composition vectors for each cluster - proportion of cells from each sample
getStructureVectors <- function(dend=NULL,current_cut=NULL){
  IDs <- get_nodes_attr(dend,"ID")
  structure_vectors <- get_nodes_attr(dend,"structure_vector")

  cut_ids_match <- match(current_cut,IDs)
  cut_ids_match <- cut_ids_match[order(cut_ids_match)]

  structure_vectors <- structure_vectors[,cut_ids_match]
  colnames(structure_vectors) <- current_cut
  return(structure_vectors)
}

#' Calculate resolution for each cut
calculateResolution <- function(dend=NULL,hc=NULL,mapping=NULL,sample_names=NULL,labels_mapping=NULL){
  resolution <- c()
  for (i in 2:attr(dend,"members")){
    current_cut <- labels_mapping[1:((i-1)*2)]
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(i-1)]))]
    membership <- getMemberhip(dend,current_cut)
    local_resolution <- sapply(sample_names,function(name) length(unique(membership[mapping==name])))
    resolution <- c(resolution, mean(local_resolution))
  }

  return(resolution)
}

#' Color node and adjacent branches based on label of cluster, number of cells and breadth
colorNode <- function(x=NULL,labels_mapping=NULL,hc=NULL,leafContent=NULL,mapping=NULL,sample_names=NULL,current_cut=NULL,
                       colors=NULL,parent_color=NULL,breadth_mapping=NULL,binar_mapping=NULL,binar_mapping2=NULL,binar_mapping3=NULL){
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

      local_binar_mapping <- binar_mapping[as.logical(memb)]
      attr(x, "cell_type_ratio1") <- sum(as.double(local_binar_mapping))/length(local_binar_mapping)

      local_binar_mapping2 <- binar_mapping2[as.logical(memb)]
      attr(x, "cell_type_ratio2") <- sum(as.double(local_binar_mapping2))/length(local_binar_mapping2)

      if (!(is.null(binar_mapping3))){
        local_binar_mapping3 <- binar_mapping3[as.logical(memb)]
        attr(x, "cell_type_ratio3") <- sum(as.double(local_binar_mapping3))/length(local_binar_mapping3)
      }
      structure_tab <- table(mapping[as.logical(memb)])

      for (name in names(structure_tab)){
        structure_vector[name] <- structure_tab[name]
      }
      structure_vector <- structure_vector/sum(structure_vector)
      attr(x, "structure_vector") <- structure_vector
      entropy_measure <- as.numeric(breadth_mapping[as.character(attr(x,"ID"))]) + 1/length(sample_names)
      attr(x, "entropy") <- entropy_measure
    }
    else{
      number_of_cells <- attr(x, "number_of_cells")
      entropy_measure <- attr(x, "entropy")
    }
    if (attr(x,"ID") %in% current_cut){
      x <- x %>%
        assign_values_to_nodes_nodePar(value = 19, nodePar = "pch") %>%
        assign_values_to_nodes_nodePar(value = 2, nodePar = "cex")  %>%
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
        x <- assign_values_to_branches_edgePar(dend = x, value = rep(parent_color, attr(x,"members")), edgePar = "col")
      }
      if (is.null(parent_color)){
        x <- set(x, "branches_col", "grey80")
        parent_color <- "grey80"
      }
    }
  }
  return(list(x=x,parent_color=parent_color))
}

#' Apply colorNode function for the entire tree
colorNodes <- function (x=NULL, ...){
  if (is.leaf(x)) {
    return(colorNode(x=x, ...))
  }

  res <- colorNode(x=x, ...)
  x <- res$x
  parent_color <- res$parent_color
  for (j in 1:2) {
    res <- colorNodes(x=x[[j]], ...)
    x[[j]] <- res$x
  }

  return(list(x=x,parent_color=parent_color))
}

colorNodeByCt <- function(x=NULL,palette=NULL,current_cut=NULL){
  entropy_measure <- attr(x, "entropy")
  cell_type_ratio1 <- attr(x, "cell_type_ratio1")
  cell_type_ratio2 <- attr(x, "cell_type_ratio2")
  if (!(is.null(attr(x, "cell_type_ratio3")))){
    cell_type_ratio3 <- attr(x, "cell_type_ratio3")
  }
  number_of_cells <- attr(x, "number_of_cells")
  if (attr(x,"ID") %in% current_cut){
    x <- x %>%
      assign_values_to_nodes_nodePar(value = 19, nodePar = "pch") %>%
      assign_values_to_nodes_nodePar(value = number_of_cells/1000, nodePar = "cex")
  }
  else{
    x <- x %>%
      assign_values_to_nodes_nodePar(value = NA, nodePar = "pch") %>%
      assign_values_to_nodes_nodePar(value = NA, nodePar = "cex")
  }
  x <- x %>%
    set("branches_lty", 1) %>%
    set("branches_lwd", entropy_measure*5)# %>%
  if (!(is.null(attr(x, "cell_type_ratio3")))){
    col <- rgb(cell_type_ratio1, cell_type_ratio2, cell_type_ratio3)
  } else {
    col <- rgb(cell_type_ratio1,0,cell_type_ratio2)
  }
  attr(x, "ct_color") <- col
  x <- assign_values_to_branches_edgePar(dend = x, value = rep(col,attr(x,"members")), edgePar = "col")

  return(list(x=x))
}

colorNodesByCt <- function (x=NULL,parent_color=NULL,palette=NULL,current_cut=NULL){
  if (is.leaf(x)){
    return(list(x=colorNodeByCt(x=x,palette=palette,current_cut=current_cut)))
  }

  x <- colorNodeByCt(x=x,palette=palette,current_cut=current_cut)$x
  for (j in 1:2){
    x[[j]] <- colorNodesByCt(x=x[[j]],palette=palette,current_cut=current_cut)$x
  }

  return(list(x=x))
}

#get list of children nodes of a node
getListOfDescendants <- function(hc=NULL,parent.ID=NULL,list.of.descendants=list()){
  children.id <- hc$merge[parent.ID,]
  if (length(list.of.descendants)>0){
    list.of.descendants[length(list.of.descendants)+1] <- children.id[1]
    list.of.descendants[length(list.of.descendants)+1] <- children.id[2]
  }
  else{
    list.of.descendants[[1]] <- children.id[1]
    list.of.descendants[[2]] <- children.id[2]
  }

  if (children.id[1] > 0){
    list.of.descendants <- getListOfDescendants(hc, children.id[1], list.of.descendants)
  }

  if (children.id[2] > 0){
    list.of.descendants <- getListOfDescendants(hc, children.id[2], list.of.descendants)
  }
  return(list.of.descendants)
}

##' get memerships vector from greedyModularityCut function results without running shiny app
##'
##' @param n.clusters number of clusters less or equal to N + 1 number of splits in greedyModularityCut results
##' @param greedy.modularity.cut.result greedyModularityCut function results
##' @export
getGreedyCutGroups <- function(n.clusters=NULL,greedy.modularity.cut.result=NULL){
  hc <- greedy.modularity.cut.result$hc
  dend <- as.dendrogram(hc)
  leafContent <- greedy.modularity.cut.result$leafContentArray
  max_no_clusters <- attr(dend,"members")
  if (n.clusters > max_no_clusters) {
    stop("n.clusters exceeds a maximum number of clusters in greedy.modularity.cut results (N + 1)")
  }

  nodes_labels <- c()

  for (i in 1:nrow(hc$merge)){
    pair <- hc$merge[i,]
    nodes_labels <- c(nodes_labels,sequentialNumeration(pair))
  }
  labels_mapping <- 1:length(nodes_labels)
  names(labels_mapping) <- rev(nodes_labels)

  dend.indexed <- colorNodes(x=dend,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent)

  current_cut <- labels_mapping[1:((n.clusters-1)*2)]
  current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(n.clusters-1)]))]

  memb <- getMemberhip(dend=dend.indexed$x,current_cut=current_cut)

  return(memb)
}

##' runs Shiny application to visualize waltktrap tree and select cut level
##'
##' @param con conos object with walktrap results
##' @param N number of top greedy splits to take
##' @param leaf.labels leaf sample label factor, for breadth calculations - must be a named factor containing all wt$names, or if wt$names is null, a factor listing cells in the same order as wt leafs
##' @param minsize minimum size of the branch (in number of leafs)
##' @param minbreadth minimum allowed breadth of a branch (measured as normalized entropy)
##' @param flat.cut whether to use a flat cut instead of a dynamic one
##' @export
conosShinyApp <- function(con, N=30, leaf.labels=NULL, tissue_mapping=NULL, tissue_factors=NULL, minsize=0, minbreadth=0, flat.cut=TRUE) {
  if (!requireNamespace("shinycssloaders", quietly = TRUE)) {
    stop("Package 'shinycssloaders' is required to run shiny app")
  }

  if (!requireNamespace("d3heatmap", quietly = TRUE)) {
    stop("Package 'd3heatmap' is required to run shiny app")
  }

  if(is.null(con$clusters$walktrap)) stop("please run findCommunities(method=walktrap.communities) to calculate walktrap clustering first")
  if(is.null(leaf.labels)) {
    # get sample labels for cells
    cl <- lapply(con$samples,getCellNames);
    leaf.labels <- as.factor(setNames(rep(1:length(cl),unlist(lapply(cl,length))),unlist(cl)))
  }

  #.GlobalEnv$con <- con
  greedy.modularity.cut.result <- conos::greedyModularityCut(wt=con$clusters$walktrap$result,N=N,leaf.labels=leaf.labels,minsize=minsize,minbreadth=minbreadth,flat.cut=flat.cut)

  leafContent <- greedy.modularity.cut.result$leafContentArray
  hc <- greedy.modularity.cut.result$hc
  dend <- as.dendrogram(hc)

  #sequential numeration of elements from hclust merge matrix
  nodes_labels <- c()
  for (i in 1:nrow(hc$merge)){
    pair <- hc$merge[i,]
    nodes_labels <- c(nodes_labels,sequentialNumeration(pair))
  }

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

  if (!is.null(tissue_factors)){
    if (!is.null(tissue_mapping)){
      if ((length(tissue_factors) == 2) & all(tissue_factors %in% unique(tissue_mapping))){
        binar_mapping1 <- tissue_mapping
        binar_mapping1[binar_mapping1==tissue_factors[1]] = 1
        binar_mapping1[binar_mapping1!="1"] = 0
        binar_mapping1 <- as.numeric(binar_mapping1)
        names(binar_mapping1) <- names(tissue_mapping)

        binar_mapping2 <- tissue_mapping
        binar_mapping2[binar_mapping2==tissue_factors[2]] = 1
        binar_mapping2[binar_mapping2!="1"] = 0
        binar_mapping2 <- as.numeric(binar_mapping2)
        names(binar_mapping2) <- names(tissue_mapping)

        binar_mapping1 <- binar_mapping1[rownames(leafContent)]
        binar_mapping2 <- binar_mapping2[rownames(leafContent)]
        binar_mapping3 <- NULL
      }
      else if ((length(tissue_factors) == 3) & all(tissue_factors %in% unique(tissue_mapping))){
        binar_mapping1 <- tissue_mapping
        binar_mapping1[binar_mapping1==tissue_factors[1]] = 1
        binar_mapping1[binar_mapping1!="1"] = 0
        binar_mapping1 <- as.numeric(binar_mapping1)
        names(binar_mapping1) <- names(tissue_mapping)

        binar_mapping2 <- tissue_mapping
        binar_mapping2[binar_mapping2==tissue_factors[2]] = 1
        binar_mapping2[binar_mapping2!="1"] = 0
        binar_mapping2 <- as.numeric(binar_mapping2)
        names(binar_mapping2) <- names(tissue_mapping)

        binar_mapping3 <- tissue_mapping
        binar_mapping3[binar_mapping3==tissue_factors[3]] = 1
        binar_mapping3[binar_mapping3!="1"] = 0
        binar_mapping3 <- as.numeric(binar_mapping3)
        names(binar_mapping3) <- names(tissue_mapping)

        binar_mapping1 <- binar_mapping1[rownames(leafContent)]
        binar_mapping2 <- binar_mapping2[rownames(leafContent)]
        binar_mapping3 <- binar_mapping3[rownames(leafContent)]
      }
      else{stop("Length of tissue_factors is not 2 or 3 or factors are not in tissue_mapping")}}
    else{stop("Lengh of tissue_mapping is zero")}
  }
  else{
    binar_mapping1 <- NULL
    binar_mapping2 <- NULL
    binar_mapping3 <- NULL
  }
  dend.colored.source <- colorNodes(x=dend,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,breadth_mapping=breadth_mapping,binar_mapping=binar_mapping1,binar_mapping2=binar_mapping2,binar_mapping3=binar_mapping3)

  #coordinates for drawing labels
  xy <- get_nodes_xy(dend)
  xy_ordered <- cbind(xy,get_nodes_attr(dend.colored.source$x,"ID"))
  xy_ordered <- xy_ordered[order(xy_ordered[,1]),]
  order_of_nodes <- xy_ordered[,3]

  #data for plots under a slider
  modularities <- cumsum(greedy.modularity.cut.result$deltaM)
  breadth <- calculateBreadth(dend=dend.colored.source$x,hc=hc,breadth_mapping=breadth_mapping,labels_mapping=labels_mapping)
  resolution <- calculateResolution(dend = dend.colored.source$x,hc=hc,mapping = mapping,sample_names = sample_names,labels_mapping=labels_mapping)

  #calculate embedding once
  #con$embedGraph()

  colors_source <- sample(rainbow(length(labels_mapping)))
  names(colors_source) <- sample(labels_mapping)

  number_of_cells_mapping <- get_nodes_attr(dend.colored.source$x,"number_of_cells")
  names(number_of_cells_mapping) <- as.character(get_nodes_attr(dend.colored.source$x,"ID"))

  reds_palette <- grDevices::colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272",
                                                "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))

  ui <- fillPage(
    tags$style(type="text/css", "body { overflow-y: scroll; }"),
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
                        tabPanel("Plot",plotOutput(outputId = "treePlot1",height = 'auto'),br()),
                        fluidRow(column(6,actionButton("do", HTML("Save<br/>membership"))),column(6,checkboxInput("color_by_ct", HTML("Gradient<br/>tissue coloring"), value = FALSE)))),
                      div(style = "height: 100%;"))),
             fluidRow(
               column(12,
                      tabPanel("Plot",align="center",plotOutput(outputId = "treePlot6")),
                      div(style = "height: 100%;")))),
      column(8,
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Tree", br(),verbatimTextOutput("txt"),shinycssloaders::withSpinner(plotOutput(outputId = "treePlot2",click = "plot_click",dblclick="plot_db_click",hover = "plot_hover"))),
                           tabPanel("Composition of clusters", shinycssloaders::withSpinner(uiOutput(outputId = "treePlot3"))),
                           tabPanel("Composition similarity", shinycssloaders::withSpinner(uiOutput(outputId = "treePlot4"))),
                           tabPanel("tSNE plots", shinycssloaders::withSpinner(plotOutput(outputId = "treePlot5")))),
               tags$head(tags$script('
                                     var dimension = [0, 0];
                                     $(document).on("shiny:connected", function(e) {
                                     dimension[0] = window.innerWidth;
                                     dimension[1] = window.innerHeight;
                                     Shiny.onInputChange("dimension", dimension);
                                     });
                                     $(window).resize(function(e) {
                                     dimension[0] = window.innerWidth;
                                     dimension[1] = window.innerHeight;
                                     Shiny.onInputChange("dimension", dimension);
                                     });'))),div(style = "width: 70%; height: 100%;"))))


  server <- function(input, output, session) { #plots under a slider
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
            if (!(any(labels_mapping[as.character(unlist(getListOfDescendants(hc,num_ID)))] %in% current_cut)) & (ID %in% current_cut)){
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
            array_of_children <- unlist(sapply(array_of_parents,function(num_ID){labels_mapping[as.character(unlist(getListOfDescendants(hc,num_ID)))]}))
            array_of_children <- array_of_children[as.numeric(names(array_of_children)) > 0]
            if (!(ID %in% array_of_children)){
              current_cut <- current_cut[!(current_cut %in% labels_mapping[as.character(unlist(getListOfDescendants(hc,num_ID)))])]
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

      colors <- colors_source[names(colors_source) %in% current_cut]
      mask <- xy_ordered[,3] %in% current_cut
      ct_mapping <- NULL
      if (!input$color_by_ct){
        dend.colored <- colorNodes(x=dend.colored.source$x,labels_mapping=labels_mapping,hc=hc,leafContent=leafContent,mapping=mapping,sample_names=sample_names,current_cut=current_cut,colors=colors,breadth_mapping=breadth_mapping,binar_mapping=binar_mapping1,binar_mapping2=binar_mapping2,binar_mapping3=binar_mapping3)
      }
      else{
        pal <- rev(viridis::viridis(11))
        dend.colored <- colorNodesByCt(x=dend.colored.source$x,palette=pal,current_cut=current_cut)
        ct_mapping <- get_nodes_attr(dend.colored$x,"ct_color")
        names(ct_mapping) <- as.character(get_nodes_attr(dend.colored$x,"ID"))
      }

      structure_vectors <- getStructureVectors(dend.colored$x,current_cut_ordered)

      heights <- get_nodes_attr(dend.colored$x,"height")
      IDs <- get_nodes_attr(dend.colored$x,"ID")

      height_mask <- (heights > 0) & (IDs %in% current_cut)

      if(any(height_mask)){height <- max(heights[height_mask])}
      else{height <- 0}

      dend.cut <- cut(dend.colored$x,h = height)

      dend.cut <- hang.dendrogram(dend.cut$upper, hang = 1)
      labels(dend.cut) <- current_cut_ordered

      memb <- getMemberhip(dend=dend.colored$x,current_cut=current_cut)
      return(list(dend.colored=dend.colored, mask=mask,structure_vectors=structure_vectors,dend.cut=dend.cut,memb=memb,colors=colors,current_cut=current_cut,ct_mapping=ct_mapping,height=height))
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
      gridExtra::grid.arrange(grobs=glist,nrow=length(glist))

    },height = 280)

    output$treePlot2 <- renderPlot({ #draw colored tree
      #color the tree and plot it with labels
      height <- dataInput()$height
      dend.colored.local <- dataInput()$dend.colored
      labels(dend.colored.local$x) <- rep(NA,attr(dend.colored.local$x,"members"))
      mask <- dataInput()$mask
      par(mar = c(0,0,0,1))
      plot(dend.colored.local$x,axes=FALSE)
      abline(h = height+0.25, col = 2, lty = 2)
      text(xy_ordered[mask,1]+0.5, xy_ordered[mask,2]+0.5, labels=xy_ordered[mask,3], col="black",cex = 0.8)
    },width = reactive(input$dimension[1]*2/3), height = reactive(input$dimension[2]*0.80)) #0.85

    output$heatmap1 <- d3heatmap::renderD3heatmap({
      struct_similarity <- dataInput()$structure_vectors
      if (is.null(click_value()) & is.null(db_click_value())){
        Rowv <- dataInput()$dend.cut
        d3heatmap::d3heatmap(t(struct_similarity),colors = reds_palette(10),Rowv=Rowv, dendrogram = "column" )}
      else{
        d3heatmap::d3heatmap(t(struct_similarity),colors = reds_palette(10), dendrogram = "column" )}
    })

    output$treePlot3 <- renderUI({
      d3heatmap::d3heatmapOutput("heatmap1", height = paste0(as.character(input$dimension[2]*0.9), "px"), width = paste0(as.character(input$dimension[1]*2/3), "px"))
    })

    output$heatmap2 <- d3heatmap::renderD3heatmap({
      struct_similarity <- dataInput()$structure_vectors
      dists <- dist(t(struct_similarity))
      if (is.null(click_value()) & is.null(db_click_value())){
        Rowv <- dataInput()$dend.cut
        d3heatmap::d3heatmap(dists,colors = rev(reds_palette((max(dists))%/%0.1)),Rowv=Rowv,Colv=Rowv)}
      else{
        d3heatmap::d3heatmap(dists,colors = rev(reds_palette((max(dists))%/%0.1)),dendrogram = "none")
      }
    })

    output$treePlot4 <- renderUI({
      d3heatmap::d3heatmapOutput("heatmap2", height = paste0(input$dimension[2]*0.9, "px"), width = paste0(input$dimension[1]*2/3, "px"))
    })

    output$treePlot5 <- renderPlot({
      memb <- dataInput()$memb
      colors <- dataInput()$colors
      suppressMessages(con$plotPanel(groups=memb, adjust.func = function(gg) gg + scale_color_manual(values = colors)))

    }, height = function() {session$clientData$output_treePlot5_width*1.2})

    output$treePlot6 <- renderPlot({
      memb <- dataInput()$memb
      colors <- dataInput()$colors
      if (input$color_by_ct){
        ct_mapping <- dataInput()$ct_mapping
        colors <- ct_mapping[names(colors)]
      }
      suppressMessages((con$plotGraph(groups = memb) + scale_color_manual(values = colors)))
    }, height = function() {session$clientData$output_treePlot6_width*0.5}, width = function() {session$clientData$output_treePlot6_width*0.91} )

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

    output$txt <- renderText({
      if(is.null(input$plot_hover)){ return("Hover a cluster\n")}
      else{
        euc_dist <- (xy_ordered[,1]-input$plot_hover$x)^2 + (xy_ordered[,2]-input$plot_hover$y)^2
        ID <- xy_ordered[which(euc_dist==min(euc_dist)),3]
        n_of_cells <- number_of_cells_mapping[as.character(ID)]
        paste(ID," cluster: ",n_of_cells," cells")
      }
    })
  }
  shinyApp(ui, server)
}


