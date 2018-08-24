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

which.row <- function(hc,label){
  if (label %in% hc$merge[,1]){
    return(which(hc$merge[,1]==label))
  }
  else if (label %in% hc$merge[,2]){
    return(which(hc$merge[,2]==label))
  }
}



calculate_breadth <- function(dend=NULL,hc=NULL,breadth_mapping=NULL,labels_mapping=NULL){
  breadth <- c()
  for (i in 2:attr(dend,"members")){
    current_cut <- labels_mapping[1:((i-1)*2)] # 4 = N
    current_cut <- current_cut[!(names(current_cut) %in% as.character((nrow(hc$merge):1)[1:(i-1)]))]
    breadth <- c(breadth,mean(breadth_mapping[as.character(current_cut)]))
  }
  return(breadth)
}

get_memberhip <- function(dend=NULL,current_cut=NULL){
  IDs <- get_nodes_attr(dend,"ID")
  memberships <- get_nodes_attr(dend,"membership")
  membership <- memberships[,IDs %in% current_cut]
  if (ncol(membership)>1){
    membership <- rowSums(membership)
  }
  return(membership)
}

get_structure_vectors <- function(dend=NULL,current_cut=NULL){
  IDs <- get_nodes_attr(dend,"ID")
  structure_vectors <- get_nodes_attr(dend,"structure_vector")
  
  cut_IDs_match <- match(current_cut,IDs)
  cut_IDs_match <- cut_IDs_match[order(cut_IDs_match)]
  
  structure_vectors <- structure_vectors[,cut_IDs_match]
  colnames(structure_vectors) <- current_cut
  return(structure_vectors)
}

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
  return(list(x=x,parent_color=parent_color))
}



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

