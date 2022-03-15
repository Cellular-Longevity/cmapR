for (file in list.files('~/projects/cmapR/R/', full.names = T)) source(file)

# create test data for package tests
mgct <- make_methylation_gct(matrix(rnorm(50), nrow=5, dimnames = list(as.character(1:5),as.character(1:10))),
                             matrix(rnorm(50), nrow=5, dimnames = list(as.character(1:5),as.character(1:10))),
                             sample.metadata = data.frame(id=as.character(1:10),row.names = as.character(1:10)),
                             chr.df = data.frame(id=as.character(1:5), row.names = as.character(1:5)))
write_gctx(mgct, '~/projects/cmapR/tests/testthat/test_n5x10.gctx')
# test read
mgct2 <- parse_gctx('~/projects/cmapR/tests/testthat/test_n5x10.gctx')

saveRDS(mgct, '~/projects/cmapR/data/ds.RData')


identical(mgct@meth_mat, mgct2@meth_mat)
