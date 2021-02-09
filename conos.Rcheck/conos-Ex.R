pkgname <- "conos"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "conos-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('conos')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("projectKNNs")
### * projectKNNs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: projectKNNs
### Title: Project a distance matrix into a lower-dimensional space.
### Aliases: projectKNNs

### ** Examples

## Not run: 
##D data(CO2)
##D CO2$Plant <- as.integer(CO2$Plant)
##D CO2$Type <- as.integer(CO2$Type)
##D CO2$Treatment <- as.integer(CO2$Treatment)
##D co <- scale(as.matrix(CO2))
##D # Very small datasets often produce a warning regarding the alias table.  This is safely ignored.
##D suppressWarnings(vis <- largeVis(t(co), K = 20, sgd_batches = 1, threads = 2))
##D suppressWarnings(coords <- projectKNNs(vis$wij, threads = 2))
##D plot(t(coords))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("projectKNNs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sgdBatches")
### * sgdBatches

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sgdBatches
### Title: Calculate the default number of batches for a given number of
###   vertices and edges.
### Aliases: sgdBatches

### ** Examples

# Observe that increasing K has no effect on processing time
N <- 70000 # MNIST
K <- 10:250
plot(K, conos:::sgdBatches(rep(N, length(K)), N * K / 2))

# Observe that processing time scales linarly with N
N <- c(seq(from = 1, to = 10000, by = 100), seq(from = 10000, to = 10000000, by = 1000))
plot(N, conos:::sgdBatches(N))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sgdBatches", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
