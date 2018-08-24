#' @export
runShinyApp <- function(wt=NULL, N=NULL, leaf.labels=NULL, minsize=0, minbreadth=0) {
  appDir <- system.file("app", package = "conos")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  .GlobalEnv$.wt <- wt
  .GlobalEnv$.N <- N
  .GlobalEnv$.leaf.labels <- leaf.labels
  .GlobalEnv$.minsize <- minsize
  .GlobalEnv$.minbreadth <- minbreadth
  .GlobalEnv$.appDir <- appDir
  on.exit(rm(list=c(".wt", ".N", ".leaf.labels", ".minsize", ".minbreadth", ".appDir"),envir=.GlobalEnv))
  shiny::runApp(appDir, display.mode = "normal")
}