#!/usr/bin/env Rscript

library("optparse")
library("biomformat")

errQuit <- function(mesg, status = 1) {
  message("Error: ", mesg)
  q(status = status)
}

option_list <- list(
  make_option(c("-i", "--input-file"),
    action = "store", type = "character",
    help = "File path to biom file to be processed"
  ),
  make_option(c("-o", "--output-file"),
    action = "store", type = "character",
    help = "File path to output file. If already exists, will be overwritten"
  ),
  make_option(c("-m", "--method"),
    action = "store", default = NULL, type = "character",
    help = "method that will pass into spieceasi, available options mb, glasso"
  ),
  make_option(c("-l", "--lambda-min-ratio"),
    action = "store", default = NULL, type = "numeric",
    help = "the scaling factor that determines the minimum sparsity/lambda parameter"
  ),
  make_option(c("-d", "--nlambda"),
    action = "store", default = NULL, type = "numeric",
    help = "nlambda"
  ),
  make_option(c("-r", "--rep-num"),
    action = "store", default = NULL, type = "integer",
    help = "pulsar-params"
  ),
  make_option(c("-c", "--ncores"),
    action = "store", default = 1, type = "integer",
    help = "pulsar-params"
  ),
  make_option(c("-t", "--thresh"),
    action = "store", default = 0.05, type = "numeric",
    help = "Threshold for StARS criterion"
  ),
  make_option(c("-b", "--subsample-ratio"),
    action = "store", default = 0.8, type = "numeric",
    help = "Subsample size for StARS"
  ),
  make_option(c("-s", "--seed"),
    action = "store", default = NULL, type = "numeric",
    help = "Set the random seed for subsample set"
  ),
  make_option("--sel-criterion",
    action = "store", default = "stars", type = "character",
    help = "Specifying criterion/method for model selection"
  ),
  make_option(c("-v", "--verbose"),
    action = "store_true", default = FALSE,
    help = "Print extra output"
  ),
  make_option("--not-lambda-log",
    action = "store_true", default = FALSE,
    help = "lambda-log"
  ),
  make_option("--not-pulsar-select",
    action = "store_true", default = FALSE,
    help = "pulsar select flag"
  ),
  make_option("--lambda-min",
    action = "store", default = NULL, type = "numeric",
    help = "lower lambda limit"
  ),
  make_option("--lambda-max",
    action = "store", default = NULL, type = "numeric",
    help = "upper lambda limit"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

inp.file <- opt$`input-file`
out.file <- opt$`output-file`
method <- opt$method
lambda.min.ratio <- opt$`lambda-min-ratio`
nlambda <- opt$nlambda
rep.num <- opt$`rep-num`
ncores <- opt$ncores
thresh <- opt$thresh
subsample.ratio <- opt$`subsample-ratio`
seed <- opt$seed
sel.criterion <- opt$`sel-criterion`
verbose <- opt$verbose
lambda.log <- !(opt$`not-lambda-log`)
lambda.min <- opt$lambda_min
lambda.max <- opt$lambda_max
pulsar.select <- !(opt$`not-pulsar-select`)

### LOAD DATA ###

raw_bioms <- lapply(strsplit(inp.file, ", ")[[1]], read_biom)
data <- lapply(raw_bioms, function(x) {
  biom_matrix <- as.matrix(biom_data(x))
  feature_ids <- rownames(biom_matrix)

  # Transposition changes features from rows to columns. Reassigning the
  # column names explicitly prevents R/SpiecEasi from replacing the original
  # BIOM observation IDs with positional names.
  result <- t(biom_matrix)
  colnames(result) <- feature_ids
  result
})

if (dir.exists(out.file)) {
  errQuit("Output filename ", out.file, " is a directory.")
} else if (file.exists(out.file)) {
  invisible(file.remove(out.file))
}

### LOAD LIBRARIES ###

suppressWarnings(library(SpiecEasi))
suppressWarnings(library(igraph))
suppressWarnings(library(Matrix))

se.out <- spiec.easi(
  data,
  method = method,
  lambda.min.ratio = lambda.min.ratio,
  nlambda = nlambda,
  pulsar.select = pulsar.select,
  pulsar.params = list(
    thresh = thresh,
    subsample.ratio = subsample.ratio,
    rep.num = rep.num,
    seed = seed,
    ncores = ncores
  )
)

# These are the original BIOM observation IDs, in the same order as the
# columns supplied to SpiecEasi. They are the source of truth for graph nodes.
feature_names <- unlist(lapply(data, colnames), use.names = FALSE)
if (length(feature_names) == 0 || any(is.na(feature_names))) {
  stop("Could not recover feature IDs from the input BIOM tables.")
}
if (anyDuplicated(feature_names)) {
  stop("Feature IDs must be unique across all input BIOM tables.")
}

features <- list(
  name = as.character(feature_names),
  Feature = as.character(feature_names)
)

if (method == "mb") {
  bm <- symBeta(getOptBeta(se.out), mode = "maxabs")
  diag(bm) <- 0
  weights <- Matrix::summary(t(bm))[, 3]
  network <- adj2igraph(Matrix::drop0(getRefit(se.out)),
    edge.attr = list(weight = weights),
    vertex.attr = features
  )
} else if (method == "glasso") {
  secor <- cov2cor(getOptCov(se.out))
  bm <- summary(triu(secor * getRefit(se.out), k = 1))
  network <- adj2igraph(getRefit(se.out),
    edge.attr = list(weight = bm[, 3]),
    vertex.attr = features
  )
} else {
  network <- adj2igraph(getRefit(se.out), vertex.attr = features)
}

# igraph may use internal node IDs when writing GraphML. Store the original
# IDs as explicit attributes; Python uses Feature rather than those internals.
V(network)$name <- as.character(feature_names)
V(network)$Feature <- as.character(feature_names)

write_graph(network, out.file, "graphml")