context("Testing mGCT class and accessor methods")

test_that("mGCT constructor works properly", {
  # initialize empty object
  g <- mGCT()
  expect_true(is(g, "mGCT"))
  expect_equal(dim(g@meth_mat), c(0, 0))
  expect_equal(dim(g@cov_mat), c(0, 0))
  expect_equal(nrow(g@rdesc), 0)
  expect_equal(nrow(g@cdesc), 0)
  # try with matrix
  g <- mGCT(meth_mat=matrix(rnorm(100), nrow=10),
           cov_mat=matrix(rnorm(100), nrow=10),
           rid=letters[1:10], cid=LETTERS[1:10])
  expect_true(is(g, "mGCT"))
  expect_equal(dim(g@meth_mat), c(10, 10))
  expect_equal(dim(g@cov_mat), c(10, 10))
  expect_equal(nrow(g@rdesc), 0)
  expect_equal(nrow(g@cdesc), 0)
})

test_that("GCT accessor methods work properly", {
    
  ds <- readRDS('~/projects/cmapR/data/ds.RData')
  # get the matrix
  
  expect_equal(ds@meth_mat, meth_mat(ds))
  expect_equal(ds@cov_mat, cov_mat(ds))
  m <- meth_mat(ds)
  # and reassign it
  tmp <- matrix(0, nrow=nrow(m), ncol=ncol(m))
  meth_mat(ds) <- tmp
  expect_equal(tmp, ds@meth_mat)
  
  # extract row ids
  expect_equal(ds@rid, ids(ds))
  x <- ids(ds)
  # and reassign them
  tmp <- as.character(seq_along(x))
  ids(ds) <- tmp
  expect_equal(tmp, ds@rid)
  
  # extract column ids
  expect_equal(ds@cid, ids(ds, dim="col"))
  x <- ids(ds, dim="col")
  # and reassign them
  tmp <- as.character(seq_along(x))
  ids(ds, dim="col") <- tmp
  expect_equal(tmp, ds@cid)
  
  # extract row meta
  expect_equal(ds@rdesc, meta(ds))
  x <- meta(ds)
  # and reassign it
  tmp <- data.frame(x=sample(letters, nrow(x), replace=TRUE))
  meta(ds) <- tmp
  expect_equal(tmp, ds@rdesc)
  
  # extract column meta
  expect_equal(ds@cdesc, meta(ds, dim="col"))
  x <- meta(ds, dim="col")
  # and reassign it
  tmp <- data.frame(x=sample(letters, nrow(x), replace=TRUE))
  meta(ds, dim="col") <- tmp
  expect_equal(tmp, ds@cdesc)
})

