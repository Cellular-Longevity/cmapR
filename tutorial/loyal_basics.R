# basic tutorial for using the CMapR package internally at Loyal

# to access files from S3 you need the aws.s3 package
if(!require(aws.s3)) install.packages('aws.s3')
# load the cmapR package
library(cmapR)

# the basic unit of this package is the 'mGCT' object, which is basically
# composed of a few different parts
# meth_mat: matrix of methylation values
# cov_mat: matrix of coverage values
# rid: row identifiers (CpG sites or genomic locations)
# cid: row identifiers (CpG sites or genomic locations)

# you may need to configure credentials if they're not available, 
# see: https://www.gormanalysis.com/blog/connecting-to-aws-s3-with-r/

# load the methylation data from the healthspan/novogene_pilot dataset
# this is a small dataset and easy to manipulate
# download is ~30Mb
object.loc <- "s3://bioinformatics-loyal/processed_methylation_data/HEALTHSPAN/novogene_pilot_RRBS/matrices_processed/methylation_filtered.gctx"
mgct <- s3read_using(FUN = parse_gctx, object = object.loc)

# to see the available data in the object, use str(mgct) or slotNames(mgct)
slotNames(mgct)
# [1] "meth_mat" "cov_mat"  "rid"      "cid"      "rdesc"    "cdesc"   
# [7] "version"  "src"  

# methylation matrix is available in mgct@meth_mat
dim(mgct@meth_mat)
# [1] 2928752       5
mgct@meth_mat[1:5,1:5]
hist(mgct@meth_mat[1:1000,], breaks=25)

# coverage matrix is available in mgct@cov_mat
dim(mgct@meth_mat)
# [1] 2928752       5
mgct@cov_mat[1:5,1:5]
hist(mgct@cov_mat[1:1000,], breaks=25)

# row (CpG site) metadata available in mgct@rdesc
# same number of rows as the matrices
dim(mgct@rdesc)
# contains CpG site information
head(mgct@rdesc)

# column (sample) metadata available in mgct@cdesc
# same number of rows as columns of the matrices
dim(mgct@cdesc)
# contains sample information 
head(mgct@rdesc)

# rid: row names. These are the same as rownames(mgct@meth_mat) and 
# as the id field in rdesc
print(mgct@rid[1:5])
identical(mgct@rid, rownames(mgct@meth_mat))
identical(mgct@rid, mgct@rdesc$id)

# cid: column names. These are the same as colnames(mgct@meth_mat) and 
# as the id field in cdesc
print(mgct@cid[1:5])
identical(mgct@cid, colnames(mgct@meth_mat))
identical(mgct@cid, mgct@cdesc$id)

# some functions for working with mGCT objects

# subsetting operations on rows or columns 
# can operate on integer indices
mgct.sub.row <- subset_gct(mgct, rid=1:500)
dim(mgct.sub.row@meth_mat)
# or on character vectors 
mgct.sub.row2 <- subset_gct(mgct, rid=mgct@rid[1:500])
dim(mgct.sub.row2@meth_mat)
identical(mgct.sub.row, mgct.sub.row2)

# columns works the same way
mgct.sub.col <- subset_gct(mgct, cid=1:3)
dim(mgct.sub.col@meth_mat)
mgct.sub.col2 <- subset_gct(mgct, cid=mgct@cid[1:3])
dim(mgct.sub.col2@meth_mat)
identical(mgct.sub.col, mgct.sub.col2)

# subset to a genomic coordinate set
#   TODO: These functions need to be incorporated into the package 
#         install, right now they're in a custom script
source('~/projects/cmapR/R/mGCT_functions.R')
mgct.chr1 <- subset_mgct_genomic(mgct, chr='chr1')
dim(mgct.chr1@meth_mat)
unique(mgct.chr1@rdesc$chrom)

# create a mGCT object from a new dataset
meth.mat <- matrix(rnorm(50), nrow=5, dimnames = list(as.character(1:5),as.character(1:10)))
cov.mat <- matrix(rnorm(50), nrow=5, dimnames = list(as.character(1:5),as.character(1:10)))
sample.metadata <-  data.frame(id=as.character(1:10),row.names = as.character(1:10))
chr.df <- data.frame(id=as.character(1:5), row.names = as.character(1:5))

mgct.new <- make_methylation_gct(meth.mat=meth.mat,
                                 cov.mat=cov.mat,
                                 sample.metadata=sample.metadata,
                                 chr.df=chr.df)
dim(mgct.new@meth_mat)
