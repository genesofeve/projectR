context("projectR")

test_that("data is proper",{

	expect_that(p.ESepiGen4c1l$mRNA.Seq,is_a('matrix'))
	expect_that(p.ESepiGen4c1l, is_a('list'))
	expect_that(AP.RNAseq6l3c3t$Amean, is_a('matrix'))
	expect_that(dim(AP.RNAseq6l3c3t),equals(c(108,5)))
	})