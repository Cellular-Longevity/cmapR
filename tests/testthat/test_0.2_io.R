context("Testing io methods")
# setwd('~/projects/cmapR/tests/testthat/')

test_that("Single column or row parsing works", {
    ds <- parse_gctx("test_n5x10.gctx", rid="1")
    expect_equal(nrow(ds@meth_mat), 1)
    expect_equal(ncol(ds@meth_mat), 10)
    expect_equal(nrow(ds@cov_mat), 1)
    expect_equal(ncol(ds@cov_mat), 10)
    expect_equal(nrow(ds@rdesc), 1)
    expect_equal(nrow(ds@cdesc), 10)
    
    ds <- parse_gctx("test_n5x10.gctx", cid="1")
    expect_equal(nrow(ds@meth_mat), 5)
    expect_equal(ncol(ds@meth_mat), 1)
    expect_equal(nrow(ds@cov_mat), 5)
    expect_equal(ncol(ds@cov_mat), 1)
    expect_equal(nrow(ds@rdesc), 5)
    expect_equal(nrow(ds@cdesc), 1)
})

test_that("GCTX parsing works", {
	ds <- parse_gctx("test_n5x10.gctx")
	expect_equal(nrow(ds@meth_mat), 5)
	expect_equal(ncol(ds@meth_mat), 10)
	expect_equal(nrow(ds@cov_mat), 5)
	expect_equal(ncol(ds@cov_mat), 10)
	expect_true(is.data.frame(ds@cdesc))
	expect_true(is.data.frame(ds@rdesc))
	expect_equal(length(ds@rid), 5)
	expect_equal(length(ds@cid), 10)
	})

test_that(
  "GCTX parsing correctly handles rid or cid that do not exist in dataset", {
  # handle the case when a subset of requested rid / cid are bogus
  expect_warning(
    ds <- parse_gctx("test_n5x10.gctx",
                      rid = c("foo", "1","2"),
                      cid = c("foo",
                              "2")))
  expect_equal(nrow(ds@meth_mat), 2)
  expect_equal(ncol(ds@meth_mat), 1)
  expect_true(is.data.frame(ds@cdesc))
  expect_true(is.data.frame(ds@rdesc))
  expect_equal(length(ds@rid), 2)
  expect_equal(length(ds@cid), 1)
  
  # fail when they're all bogus
  expect_error(ds <- parse_gctx("test_n5x10.gctx",
                                  rid = c("foo", "bar"),
                                  cid = c("foo", "bar")))
  
  # same as above, using numeric indices
  expect_warning(ds <- parse_gctx("test_n5x10.gctx",
                                  rid = c(0, 5),
                                  cid = c(0, 5)))
  expect_equal(nrow(ds@meth_mat), 1)
  expect_equal(ncol(ds@meth_mat), 1)
  expect_true(is.data.frame(ds@cdesc))
  expect_true(is.data.frame(ds@rdesc))
  expect_equal(length(ds@rid), 1)
  expect_equal(length(ds@cid), 1)
})

test_that(
  "Writing GCTX works when row or column descriptors have just one column", {
  ds <- parse_gctx("test_n5x10.gctx")
  # set rdesc and cdesc to single-column data.frames
  ds@rdesc <- data.frame("id"=ds@rdesc[, 1])
  ds@cdesc <- data.frame("id"=ds@cdesc[, 1])
  write_gctx(ds, "foo.gctx", appenddim = FALSE)
  # remove the file
  file.remove("foo.gctx")
})

test_that("fix_datatypes correctly handles variety of data types", {
  # read a table of annotations and force all classes to be character initially
  cdesc <- read.delim("test_cdesc.txt", colClasses = "character")
  # run the fixing
  fixed <- fix_datatypes(cdesc)
  # make sure certain columns are of certain types
  # these fields should be characters
  expect_true(is.character(fixed$pert_id))
  expect_true(is.character(fixed$pert_iname))
  expect_true(is.character(fixed$pert_type))
  # these should all be ints
  expect_true(is.integer(fixed$pert_time))
  expect_true(is.integer(fixed$qc_slope))
  # these should be numeric
  expect_true(is.numeric(fixed$qc_f_logp))
  expect_true(is.numeric(fixed$qc_iqr))
  # sci_note is stored on disk in exponential format, which
  # should be converted to numeric. 
  expect_true(is.numeric(fixed$sci_note))
})

test_that("update_gctx works correctly", {
  stop('THIS FEATURE TODO')
  # make a copy of the example dataset
  fpath <- "test_copy_n5x10.gctx"
  if (file.exists(fpath)) file.remove(fpath)
  file.copy("test_n5x10.gctx", fpath)
  # modify rows 3-7, columns 2-4 to contain all zeros
  m <- matrix(0, nrow=5, ncol=3)
  # update using integer indices
  update_gctx(m, ofile=fpath, rid=3:7, cid=2:4)
  tmp <- parse_gctx(fpath)
  tmp_m <- tmp@meth_mat[3:7, 2:4]
  dimnames(tmp_m) <- NULL
  expect_identical(m, tmp_m)
  # update using character ids
  m2 <- matrix(1, nrow=5, ncol=3)
  rid <- read_gctx_ids("test_n5x10.gctx", dim="row")
  cid <- read_gctx_ids("test_n5x10.gctx", dim="col")
  update_gctx(m2, ofile=fpath, rid=rid[3:7], cid=cid[2:4])
  tmp2 <- parse_gctx(fpath)
  tmp_m2 <- tmp2@meth_mat[3:7, 2:4]
  dimnames(tmp_m2) <- NULL
  expect_identical(m2, tmp_m2)
  # try updating indices that don't exist in the dataset
  # should produce an error
  expect_error(update_gctx(m2, ofile=fpath, rid=3:7, cid=20:30))
  # try updating indices that don't correspond to dims of array
  # should produce an error
  expect_error(update_gctx(m2, ofile=fpath, rid=3:7, cid=1:2))
  expect_error(update_gctx(rep(0, 10), ofile=fpath, rid=3:7, cid=1:2))
  if (file.exists(fpath)) file.remove(fpath)
})
