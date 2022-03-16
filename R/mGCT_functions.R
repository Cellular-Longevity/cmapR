# helper functions for working with mGCT obejects

# subset to a particular genomic interval
# intervals will be inclusive of start and end
# currently only a single chromosome supported
subset_mgct_genomic <- function(mgct, chr, start=0, end=1e99){
    rids <- mgct@rid[mgct@rdesc$chr == chr & mgct@rdesc$start >= start &
                         mgct@rdesc$end <= end]
    return(subset_gct(mgct, rid = rids))
}


# construct new mGCT from four pieces of data
# matrix of methylation percentage, matrix of methylation coverage, sample metadata and position metadata
make_methylation_gct <- function(meth.mat, cov.mat, sample.metadata, chr.df){
    # chr.df must have index column
    if (! (('index' %in% colnames(chr.df))| ('id' %in% colnames(chr.df)))){
        chr.df$id <- paste(chr.df$chrom, chr.df$start, chr.df$end, sep='_')
    }
    
    # ensure dimensions match up 
    if (nrow(chr.df) != nrow(meth.mat) | nrow(chr.df) != nrow(cov.mat) | 
        ncol(meth.mat) != ncol(cov.mat)){
        stop('Dimension mismatch in creating mGCT')
    }
    
    if((!identical(colnames(meth.mat), colnames(cov.mat))) | 
       (!identical(colnames(meth.mat), rownames(sample.metadata)))){
        stop('colnames of meth.mat, cov.mat and rownames of sample.metadata must be identical')
    }
    
    mgct <- new('mGCT')
    mgct@rdesc <- chr.df
    rownames(mgct@rdesc) <- chr.df$id
    mgct@rid <- chr.df$id
    mgct@meth_mat <- meth.mat
    rownames(mgct@meth_mat) <- chr.df$id
    mgct@cov_mat <- cov.mat
    rownames(mgct@cov_mat) <- chr.df$id
    
    mgct@cid <- colnames(meth.mat)
    mgct@cdesc <- sample.metadata
    return(mgct)
}
