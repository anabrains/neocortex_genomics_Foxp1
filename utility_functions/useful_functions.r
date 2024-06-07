#### functions

# Allows us to save objects with a new name defined by a string
saveit <- function(..., string, file) {
  x <- list(...)
  names(x) <- string
  save(list=names(x), file=file, envir=list2env(x))
}


# Basic function to convert mouse to human gene names
convertMouseToHuman <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://feb2021.archive.ensembl.org", path = "/biomart/martservice")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://feb2021.archive.ensembl.org", path = "/biomart/martservice")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx),2)
  return(humanx)
}

# Basic function to convert human to mouse gene names
convertHumanToMouse <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  mousex <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(mousex),2)
  return(mousex)
}
## Ex usage
#genes <- convertMouseGeneList(humGenes)
#genes <- convertHumanToMouse(humGenes)

########## modified from rio to fast write 10x counts from snATAC signac object ##   h5 version not tested #takes in custom name
write10xCountsAtac <- function(path, x, barcodes=colnames(x), gene.id=rownames(x), gene.symbol=gene.id,
                               overwrite=FALSE, type=c("auto", "sparse", "HDF5"), genome="unknown", version=c("2", "3"),
                               chemistry="snATAC", original.gem.groups=1L, library.ids="custom",
                               custom_name = NULL) {
  .type_chooser <- function(path, type) {
    if (type=="auto") {
      type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
    }
    type
  }
  
  # Doing all the work on a temporary location next to 'path', as we have permissions there.
  # This avoids problems with 'path' already existing.
  temp.path <- tempfile(tmpdir=dirname(path)) 
  on.exit({ 
    if (file.exists(temp.path)) { unlink(temp.path, recursive=TRUE) } 
  })
  
  # Checking the values.
  if (length(gene.id)!=length(gene.symbol) || length(gene.id)!=nrow(x)) {
    stop("lengths of 'gene.id' and 'gene.symbol' must be equal to 'nrow(x)'")
  }
  if (ncol(x)!=length(barcodes)) { 
    stop("'barcodes' must of of the same length as 'ncol(x)'")
  }
  
  # Determining what format to save in.
  version <- match.arg(version)
  type <- .type_chooser(path, match.arg(type))
  if (type=="sparse") {
    .write_sparse(temp.path, x, barcodes, gene.id, gene.symbol, version=version, custom_name = custom_name)
  } else {
    .write_hdf5(temp.path, genome, x, barcodes, gene.id, gene.symbol, version=version, custom_name = custom_name)
  }
  
  # We don't put this at the top as the write functions might fail; 
  # in which case, we would have deleted the existing 'path' for nothing.
  if (overwrite) {
    unlink(path, recursive=TRUE)
  } else if (file.exists(path)) { 
    stop("specified 'path' already exists")
  }
  file.rename(temp.path, path)
  return(invisible(TRUE))
}

.write_sparse <- function(path, x, barcodes, gene.id, gene.symbol, version="2", custom_name = NULL) {
  dir.create(path, showWarnings=FALSE)
  gene.info <- data.frame(gene.id, stringsAsFactors=FALSE)
  
  if (version=="3") {
    mhandle <- if (!is.null(custom_name)) file.path(path, paste0(custom_name, "_matrix.mtx")) else file.path(path, "matrix.mtx")
    bhandle <- gzfile(if (!is.null(custom_name)) file.path(path, paste0(custom_name, "_barcodes.tsv.gz")) else file.path(path, "barcodes.tsv.gz"), open="wb")
    fhandle <- gzfile(if (!is.null(custom_name)) file.path(path, paste0(custom_name, "_features.tsv.gz")) else file.path(path, "features.tsv.gz"), open="wb")
    on.exit({
      close(bhandle)
      close(fhandle)
    })
  } else {
    mhandle <- if (!is.null(custom_name)) file.path(path, paste0(custom_name, "_matrix.mtx")) else file.path(path, "matrix.mtx")
    bhandle <- if (!is.null(custom_name)) file.path(path, paste0(custom_name, "_barcodes.tsv")) else file.path(path, "barcodes.tsv")
    fhandle <- if (!is.null(custom_name)) file.path(path, paste0(custom_name, "_features.tsv")) else file.path(path, "features.tsv")
  }
  
  writeMM(x, file=mhandle)
  write(barcodes, file=bhandle)
  write.table(gene.info, file=fhandle, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  if (version=="3") {
    # Annoyingly, writeMM doesn't take connection objects.
    R.utils::gzip(mhandle)
  }
  
  return(NULL)
}

.write_hdf5 <- function(path, genome, x, barcodes, gene.id, gene.symbol, version="3",
                        chemistry="Single Cell 3' v3", original.gem.groups=1L, library.ids="custom",
                        custom_name = NULL)
{
  path <- path.expand(path) # protect against tilde's.
  h5createFile(path)
  
  if (version=="3") {
    group <- "matrix"
  } else {
    group <- genome
  }
  h5createGroup(path, group)
  
  h5write(barcodes, file=path, name=paste0(group, "/", if (!is.null(custom_name)) paste0(custom_name, "_barcodes") else "barcodes"))
  
  if (version=="3") {
    h5createGroup(path, file.path(group, "features"))
    
    h5write(gene.id, file=path, name=paste0(group, "/features/", if (!is.null(custom_name)) paste0(custom_name, "_id") else "id"))
    h5write(gene.symbol, file=path, name=paste0(group, "/features/", if (!is.null(custom_name)) paste0(custom_name, "_name") else "name"))
    h5writeAttribute(chemistry, h5obj=h5g, name="chemistry_description")
    h5writeAttribute("matrix", h5obj=h5g, name="filetype")
    h5writeAttribute(library.ids, h5obj=h5g, name="library_ids")
    h5writeAttribute(original.gem.groups, h5obj=h5g, name="original_gem_groups")
    h5writeAttribute(as.integer(version) - 1L, h5obj=h5g, name="version") # this is probably correct.
    
  } else {
    h5write(gene.id, file=path, name=paste0(group, "/genes"))
    h5write(gene.symbol, file=path, name=paste0(group, "/gene_names"))
  }
  
  # Saving matrix information.
  x <- as(x, "dgCMatrix")
  h5write(x@x, file=path, name=paste0(group, "/data"))
  h5write(dim(x), file=path, name=paste0(group, "/shape"))
  h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
  h5write(x@p, file=path, name=paste0(group, "/indptr"))
  
  return(NULL)
}
