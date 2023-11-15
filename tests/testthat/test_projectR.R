context("projectR")

test_that("data is proper",{

	expect_that(p.ESepiGen4c1l$mRNA.Seq,is_a('matrix'))
	expect_that(p.ESepiGen4c1l, is_a('list'))
	expect_that(length(p.ESepiGen4c1l),equals(6))
	expect_that(map.ESepiGen4c1l, is_a('data.frame'))
	expect_true(all(dim(map.ESepiGen4c1l) == c(93,9)))
	expect_that(map.ESepiGen4c1l$GeneSymbols, is_a('character'))
	expect_that(AP.RNAseq6l3c3t, is_a(c('list','CoGAPS')))
	expect_that(length(AP.RNAseq6l3c3t),equals(12))
	expect_that(AP.RNAseq6l3c3t$Amean, is_a('matrix'))
	expect_true(all(dim(AP.RNAseq6l3c3t$Amean) == c(108,5)))
	expect_that(pd.ESepiGen4c1l,is_a('data.frame'))
	expect_true(all(dim(pd.ESepiGen4c1l) == c(9,2)))
	expect_that(pd.RNAseq6l3c3t,is_a('data.frame'))
	expect_true(all(dim(pd.RNAseq6l3c3t) == c(54,38)))
	expect_true(all(dim(CR.RNAseq6l3c3t) == c(54,5)))
	expect_that(multivariateAnalysisR_seurat_test, is_a('Seurat'))
	})

test_that("results are as expected",{
	#CoGAPS check
	library("CoGAPS")
	# CR.RNAseq6l3c3t <- CoGAPS(p.RNAseq6l3c3t, params = new("CogapsParams", nPatterns=5))
	pr_cgps <- projectR(data=p.ESepiGen4c1l$mRNA.Seq,loadings=CR.RNAseq6l3c3t,
		dataNames=map.ESepiGen4c1l[["GeneSymbols"]])
	expect_that(pr_cgps, is_a('matrix'))
	expect_true(all(dim(pr_cgps) == c(5,9)))
	expect_true(all(pr_cgps != 0))
	expect_true(all(!is.na(pr_cgps)))


	#cluster2patter check
	k.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,
		loadings=kmeans(p.RNAseq6l3c3t, 4),
		dataNames=map.ESepiGen4c1l$GeneSymbols,
		loadingsNames=rownames(p.RNAseq6l3c3t),
		full=FALSE, sourceData=p.RNAseq6l3c3t)
	expect_true(all(dim(k.ESepiGen4c1l) == c(4,9)))
	expect_true(all(k.ESepiGen4c1l != 0))
	expect_true(all(!is.na(k.ESepiGen4c1l)))

	#pclust check
	pca.RNAseq6l3c3t<-prcomp(t(p.RNAseq6l3c3t))
	pca.ESepiGen4c1l<-projectR(data=p.ESepiGen4c1l$mRNA.Seq,
		loadings=pca.RNAseq6l3c3t,dataNames=map.ESepiGen4c1l[["GeneSymbols"]])
	expect_true(all(dim(pca.ESepiGen4c1l) == c(54,9)))
	expect_true(all(pca.ESepiGen4c1l != 0))
	expect_true(all(!is.na(pca.ESepiGen4c1l)))
	
	#multivariateAnalysisR check
	output <- multivariateAnalysisR(seuratobj = multivariateAnalysisR_seurat_test, 
	                                patternKeys = list("Pattern_1", "Pattern_2"), 
	                                dictionaries = list(
	                                  list("stage" = "E18"), 
	                                  list("stage" = "Adult")
	                                  )
	                                )
	expect_is(output, "list")
	expect_length(output, 2)
	expect_true("patternKey" %in% names(output[[1]]))
	expect_true("ANOVA" %in% names(output[[1]]))
	expect_true("CI" %in% names(output[[1]]))
	
	})
