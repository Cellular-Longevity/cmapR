########################################
### mGCT class and method definitions ##
########################################

#' An S4 class to represent a mGCT object
#' 
#' @slot meth_mat a numeric matrix of methylation percentage
#' @slot cov_mat a numeric matrix of methylation coverage
#' @slot rid a character vector of row ids
#' @slot cid a character vector of column ids
#' @slot rdesc a \code{data.frame} of row descriptors
#' @slot rdesc a \code{data.frame} of column descriptors
#' @slot src a character indicating the source (usually file path) of the data
#' 
#' @description The mGCT class is an extension of the normal GCT class, which
#'   contains slots for two data matrices, and additional data matrix parsing
#'   utilities. The \code{meth_mat} slot contains methylation percentages, 
#'   the \code{cov_mat} slot contains methylation coverage,
#'   \code{rdesc} and \code{cdesc} slots contain data frames with
#'   annotations about the rows and columns, respectively
#'   
#' @seealso \code{\link{parse_gctx}},
#' \code{\link{write_gctx}}, \code{\link{read_gctx_meta}},
#' \code{\link{read_gctx_ids}}
#' @seealso visit \url{http://clue.io/help} for more information on the
#'   mGCT format
methods::setClass("mGCT",
                  methods::representation(
                    meth_mat = "matrix",
                    cov_mat = "matrix",
                    rid = "character",
                    cid = "character",
                    rdesc = "data.frame",
                    cdesc = "data.frame",
                    version = "character",
                    src = "character"
                  )
)


## ----set up methods for checking mGCT validity----
methods::setValidity("mGCT",
                     function(object) {
                       # check whether dimensions of various
                       # slots are in sync
                       m_m <- meth_mat(object)
                       m_c <- cov_mat(object)
                       rid <- ids(object)
                       cid <- ids(object, dim="column")
                       rdesc <- meta(object)
                       cdesc <- meta(object, dim="column")
                       nrows_m <- nrow(m_m)
                       nrows_c <- nrow(m_c)
                       ncols_m <- ncol(m_m)
                       ncols_c <- ncol(m_c)
                       if (nrows_m != nrows_c){
                        return(
                          "methylation and coverage matrix must have same number of rows")
                       }
                       if (ncols_m != ncols_c){
                        return(
                          "methylation and coverage matrix must have same number of cols")
                       }
                       if (nrows_m != length(rid)) {
                         return(
                      "rid must be the same length as number of matrix rows")
                       }
                       if (ncols_m != length(cid)) {
                         return(
                    "cid must be the same length as number of matrix columns")
                       }
                       if (any(duplicated(cid))) {
                         return("cid must be unique")
                       }
                       if (any(duplicated(rid))) {
                         return("rid must be unique")
                       }
                       if (nrow(cdesc) != ncols_m & nrow(cdesc) != 0) {
                         return(paste(
                           "cdesc must either have 0 rows or the",
                           "same number of rows as matrix has columns"))
                       }
                       if (nrow(rdesc) != nrows_m & nrow(rdesc) != 0) {
                         return(paste(
                           "rdesc must either have 0 rows or the same number",
                           "of rows as matrix has rows"))
                       }
                       else {
                         return(T)
                       }
                     }
)

## ----define the initialization method for the mGCT class----
methods::setMethod("initialize",
                   signature = "mGCT",
                   definition = function(.Object, meth_mat=NULL,
                                         cov_mat=NULL, rdesc=NULL,
                                         cdesc=NULL,
                                         src=NULL, rid=NULL, cid=NULL,
                                         matrix_only=FALSE) {
                     # if we were supplied a matrix and annotations, use them
                     if (!is.null(meth_mat)) {
                       .Object@meth_mat <- meth_mat
                       .Object@cov_mat <- cov_mat
                       # if given rid and cid, use those as well
                       if (!is.null(rid)) {
                         .Object@rid <- rid
                       } else {
                         .Object@rid <- rownames(meth_mat)
                       }
                       if (!is.null(cid)) {
                         .Object@cid <- cid
                       } else {
                         .Object@cid <- colnames(meth_mat)
                       }
                       if (!is.null(rdesc)) {
                         .Object@rdesc <- rdesc
                       }
                       if (!is.null(cdesc)) {
                         .Object@cdesc <- cdesc
                       }
                       # make sure rid, cid and dimnames of mat in sync
                       dimnames(.Object@meth_mat) <- list(.Object@rid, .Object@cid)
                       dimnames(.Object@cov_mat) <- list(.Object@rid, .Object@cid)
                     } else if (!is.null(src)) {
                       # we were not given a matrix, were we given a src file?
                       # check to make sure it's a .gctx
                       if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
                         stop("Either a .gct or .gctx file must be given")
                       if (grepl(".gct$", src)) {
                         stop('GCTX file must be used with this implementation')
                       } else { 
                         # parse the .gctx
                         message("reading ", src)
                         .Object@src <- src
                         # if the rid's or column id's are .grp files,
                         # read them in
                         if ( length(rid) == 1 && grepl(".grp$", rid) )
                           rid <- parse_grp(rid)
                         if ( length(cid) == 1 && grepl(".grp$", cid) )
                           cid <- parse_grp(cid)
                         # get all the row and column ids
                         all_rid <- read_gctx_ids(src, dim="row")
                         all_cid <- read_gctx_ids(src, dim="col")
                         # if rid or cid specified, read only those rows/columns
                         # if already numeric, use as is
                         # else convert to numeric indices
                         processed_rids <- process_ids(rid, all_rid, type="rid")
                         processed_cids <- process_ids(cid, all_cid, type="cid")
                         
                         # get number of dimensions in the hdf5 matrix
                         lsdf <- rhdf5::h5ls(src)
                         ndim <- length(strsplit(lsdf[lsdf$group == '/0/DATA/0' & lsdf$name=='matrix','dim'], split=' x ')[[1]])
                         
                         # catch case of 2D matrix at this location in the hdf5 file
                         # in which we place the data in the meth_mat slot
                         # and fill the cov_mat slot with NA
                         if (ndim == 2){
                             warning('Reading old formatted gctx, matrix information will be in object@meth_mat')
                             # read the 2D array first, then assign to the respective slots
                             meth_array <- rhdf5::h5read(
                                 src, name="0/DATA/0/matrix",
                                 index=list(processed_rids$idx, processed_cids$idx))
                             print(dim(meth_array))
                            .Object@meth_mat <- meth_array
                            .Object@cov_mat <- matrix(NA, nrow = nrow(.Object@meth_mat), ncol=ncol(.Object@meth_mat))
                         } else if (ndim == 3){
                             # read the 3D array first, then assign to the respective slots
                             meth_array <- rhdf5::h5read(
                                 src, name="0/DATA/0/matrix",
                                 index=list(processed_rids$idx, processed_cids$idx, NULL))
                             .Object@meth_mat <- matrix(meth_array[,,1,drop=T], nrow=dim(meth_array)[1])
                             .Object@cov_mat <- matrix(meth_array[,,2,drop=T], nrow=dim(meth_array)[1])

                         } else { 
                            stop('Dimension mismatch at 0/DATA/0/matrix in hdf5 file')
                         }
                         # set the row and column ids, casting as characters
                         .Object@rid <- processed_rids$ids
                         .Object@cid <- processed_cids$ids
                         rownames(.Object@meth_mat) <- processed_rids$ids
                         colnames(.Object@meth_mat) <- processed_cids$ids
                         rownames(.Object@cov_mat) <- processed_rids$ids
                         colnames(.Object@cov_mat) <- processed_cids$ids
                         # get the meta data
                         if (!matrix_only) {
                           .Object@rdesc <- read_gctx_meta(
                             src, dim="row",
                             ids=processed_rids$ids)
                           .Object@cdesc <- read_gctx_meta(
                             src, dim="col",
                             ids=processed_cids$ids)
                         }
                         else {
                           .Object@rdesc <- data.frame(id=.Object@rid,
                                                       stringsAsFactors=FALSE)
                           .Object@cdesc <- data.frame(id=.Object@cid,
                                                       stringsAsFactors=FALSE)
                         }
                         # close any open handles and return the object
                         if(utils::packageVersion('rhdf5') < "2.23.0") {
                           rhdf5::H5close()
                         } else {
                           rhdf5::h5closeAll()
                         }
                         message("done")
                       }
                     }
                     # finally, make sure object is valid before returning
                     ok <- methods::validObject(.Object)
                     return(.Object)
                   }
)

#' Initialize an object of class \code{mGCT}
#' @param meth_mat a matrix of methylation percentages
#' @param cov_mat a matrix of sequencing coverages
#' @param rdesc a \code{data.frame} of row metadata
#' @param cdesc a \code{data.frame} of column metadata
#' @param src path to a mGCT file to read
#' @param rid vector of character identifiers for rows
#' @param cid vector of character identifiers for columns
#' @param matrix_only logical indicating whether to read just the matrix
#'   data from \code{src}
#'   
#' @details 
#'   If \code{meth_mat} is provided, \code{rid} and \code{cid} are treated as
#'   the row and column identifiers for the matrix and are assigned to the
#'   \code{rid} and \code{cid} slots of the \code{mGCT} object.
#'   
#'   If \code{mat} is not provided but \code{src} is provided,
#'   \code{rid} and \code{cid} are treated as filters. Data will be read from
#'   the file path provided to \code{src} and will then be restricted to the
#'   character ids or integer indices provided to \code{rid} and \code{cid}.
#'   In a similar manner, \code{matrix_only} controls whether the
#'   row and column metadata are also read from the \code{src} file path.
#'   
#' @returns a \code{mGCT} object
#' @examples 
#' # an empty object
#' (g <- mGCT())
#' # with a matrix
#' # note we must specify row and column ids
#' (g <- mGCT(mat=matrix(rnorm(100), nrow=10),
#'           rid=letters[1:10], cid=letters[1:10]))
#' # from file
#' gct_file <- system.file("extdata", "modzs_n25x50.gctx", package="cmapR")
#' (g <- mGCT(src=gct_file))
#' @family mGCTX parsing functions
#' @importFrom methods new
#' @export
mGCT <- function(meth_mat=NULL, cov_mat=NULL, rdesc=NULL, cdesc=NULL,
                src=NULL, rid=NULL, cid=NULL,
                matrix_only=FALSE) {
  methods::new("mGCT", meth_mat=meth_mat, cov_mat=cov_mat, rdesc=rdesc, cdesc=cdesc,
      src=src, rid=rid, cid=cid, matrix_only=matrix_only)
}

###########################################
### accessor functions for mGCT objects  ###
###########################################

# set method for displaying a mGCT object
# just use the 'str' function to show its structure
setMethod("show", methods::signature("mGCT"), function(object) {
  utils::str(object)
})

#' Extract or set the matrix of mGCT object
#' @param g the mGCT object
#' @param value a numeric matrix
#' @return a matrix
#' @examples 
#' # get the matrix
#' m <- mat(ds)
#' # set the matrix
#' mat(ds) <- matrix(0, nrow=nrow(m), ncol=ncol(m))
#' @family mGCT accessor methods
#' @export
methods::setGeneric("meth_mat", function(g) {
  standardGeneric("meth_mat")
})
methods::setGeneric("cov_mat", function(g) {
  standardGeneric("cov_mat")
})
#' @rdname meth_mat
methods::setMethod("meth_mat", "mGCT", function(g) g@meth_mat)
#' @export
#' @rdname cov_mat
methods::setMethod("cov_mat", "mGCT", function(g) g@cov_mat)
#' @export
#' @rdname meth_mat
methods::setGeneric("meth_mat<-", function(g, value) {
  standardGeneric("meth_mat<-")
})
#' @rdname cov_mat
methods::setGeneric("cov_mat<-", function(g, value) {
  standardGeneric("cov_mat<-")
})
#' @rdname meth_mat
methods::setMethod("meth_mat<-", "mGCT", function(g, value) {
  g@meth_mat <- value
  methods::validObject(g)
  return(g)
})
#' @rdname cov_mat
methods::setMethod("cov_mat<-", "mGCT", function(g, value) {
  g@cov_mat <- value
  methods::validObject(g)
  return(g)
})

#' Extract the or set row or column ids of a mGCT object
#' @param g the mGCT object
#' @param dimension the dimension to extract/update ['row' or 'column']
#' @param value a character vector
#' @return a vector of row ids
#' @examples 
#' # extract rids
#' rids <- ids(ds)
#' # extract column ids
#' cids <- ids(ds, "column")
#' # set rids
#' ids(ds) <- as.character(1:length(rids))
#' # set cids
#' ids(ds, "column") <- as.character(1:length(cids))
#' @family mGCT accessor methods
#' @export
methods::setGeneric("ids", function(g, dimension="row")  {
  standardGeneric("ids")
})
#' @rdname ids
methods::setMethod("ids", "mGCT", function(g, dimension="row") {
  dimension <- tolower(dimension)
  if (dimension == "col") dimension <- "column"
  stopifnot(dimension %in% c("row", "column"))
  switch(dimension, row=g@rid, column=g@cid)
})
#' @export
#' @rdname ids
methods::setGeneric("ids<-", function(g, dimension="row", value)  {
  standardGeneric("ids<-")
})
#' @rdname ids
methods::setMethod("ids<-", "mGCT", function(g, dimension="row", value) {
  dimension <- tolower(dimension)
  if (dimension == "col") dimension <- "column"
  stopifnot(dimension %in% c("row", "column"))
  if (dimension == "row") {
    g@rid <- value
  } else {
    g@cid <- value
  }
  methods::validObject(g)
  return(g)
})

#' Extract the or set metadata of a mGCT object
#' @param g the mGCT object
#' @param dimension the dimension to extract/update ['row' or 'column']
#' @param value a data.frame
#' @return a data.frame
#' @examples 
#' # extract rdesc
#' rdesc <- meta(ds)
#' # extract cdesc
#' cdesc <- meta(ds, dim="column")
#' # set rdesc
#' meta(ds) <- data.frame(x=sample(letters, nrow(rdesc), replace=TRUE))
#' # set cdesc
#' meta(ds, dim="column") <- data.frame(x=sample(letters, nrow(cdesc),
#'   replace=TRUE))
#' @family mGCT accessor methods
#' @export
methods::setGeneric("meta", function(g, dimension="row")  {
  standardGeneric("meta")
})
#' @rdname meta
methods::setMethod("meta", "mGCT", function(g, dimension="row") {
  dimension <- tolower(dimension)
  if (dimension == "col") dimension <- "column"
  stopifnot(dimension %in% c("row", "column"))
  switch(dimension, row=g@rdesc, column=g@cdesc)
})
#' @export
#' @rdname meta
methods::setGeneric("meta<-", function(g, dimension="row", value)  {
  standardGeneric("meta<-")
})
#' @rdname meta
methods::setMethod("meta<-", "mGCT", function(g, dimension="row", value) {
  dimension <- tolower(dimension)
  if (dimension == "col") dimension <- "column"
  stopifnot(dimension %in% c("row", "column"))
  if (dimension == "row") {
    g@rdesc <- value
  } else {
    g@cdesc <- value
  }
  methods::validObject(g)
  return(g)
})


###########################################
###  cast GCT as SummarizedExperiment   ###
###########################################

#' as("GCT", "SummarizedExperiment")
#' 
#' Create SummarizedExperiment object from GCT object.
#' 
#' @examples
#'
#' se <- as(ds, "SummarizedExperiment")
#' 
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods validObject
setAs("mGCT", "SummarizedExperiment", function(from) {
  stop('TODO')
  stopifnot(methods::validObject(from))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat(from)), 
    colData = meta(from, dimension="column"),
    rowData = meta(from))
})
