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

#projectionDriveR check
#test that expected output is present and in correct format

test_that("results are correctly formatted for confidence interval mode",{
  
  pattern_to_weight <- "Pattern.24"
  drivers <- projectionDriveR(microglial_counts, #expression matrix
                              glial_counts, #expression matrix
                              loadings = retinal_patterns, #feature x pattern dataframe
                              loadingsNames = NULL,
                              pattern_name = pattern_to_weight, #column name
                              pvalue = 1e-5, #pvalue before bonferroni correction
                              display = T,
                              normalize_pattern = T,  #normalize feature weights
                              mode = "CI") #set to confidence interval mode
#check output is in list format
expect_is(drivers, "list")

#check length of dfs
expect_length(drivers, 6)
expect_length(drivers$mean_ci, 3)
expect_length(drivers$weighted_mean_ci, 3)

#check that genes used for calculations overlap both datasets and loadings
expect_true(unique(drivers$mean_ci$gene %in% rownames(microglial_counts)))
expect_true(unique(drivers$mean_ci$gene %in% rownames(glial_counts)))
expect_true(unique(drivers$mean_ci$gene %in% rownames(retinal_patterns)))
expect_true(unique(drivers$weighted_mean_ci$gene %in% rownames(microglial_counts)))
expect_true(unique(drivers$weighted_mean_ci$gene %in% rownames(glial_counts)))
expect_true(unique(drivers$weighted_mean_ci$gene %in% rownames(retinal_patterns)))

#name and class checks
expect_true("mean_ci" %in% names(drivers))
expect_is(drivers$mean_ci, "data.frame")

expect_true("weighted_mean_ci" %in% names(drivers))
expect_is(drivers$mean_ci, "data.frame")

expect_true("normalized_weights" %in% names(drivers))
expect_is(drivers$normalized_weights, "numeric")

expect_true("sig_genes" %in% names(drivers))
expect_is(drivers$sig_genes, "list")
expect_length(drivers$sig_genes, 3)

expect_true(unique(c("unweighted_sig_genes", "weighted_sig_genes", "significant_shared_genes") %in% names(drivers$sig_genes)))

expect_is(drivers$sig_genes$unweighted_sig_genes, "character")

expect_is(drivers$sig_genes$weighted_sig_genes, "character")

expect_is(drivers$sig_genes$significant_shared_genes, "character")

expect_true("meta_data" %in% names(drivers))
expect_is(drivers$meta_data, "list")
expect_length(drivers$meta_data, 2)

#check that matrix names are proper and match source names
expect_true(deparse(substitute(microglial_counts)) == drivers$meta_data$test_matrix)
expect_type(drivers$meta_data$test_matrix, "character")

expect_true(deparse(substitute(glial_counts)) == drivers$meta_data$reference_matrix)
expect_type(drivers$meta_data$reference_matrix, "character")

#check that plot length is correct
expect_true("plotted_ci" %in% names(drivers))
expect_length(drivers$plotted_ci, 2)
})

test_that("results are correctly formatted for P value mode",{
  
  pattern_to_weight <- "Pattern.24"
  drivers <- projectionDriveR(microglial_counts, #expression matrix
                              glial_counts, #expression matrix
                              loadings = retinal_patterns, #feature x pattern dataframe
                              loadingsNames = NULL,
                              pattern_name = pattern_to_weight, #column name
                              pvalue = 1e-5, #pvalue before bonferroni correction
                              display = T,
                              normalize_pattern = T,  #normalize feature weights
                              mode = "PV") #set to p value mode
  #check output is in list format
  expect_is(drivers, "list")
  
  #check length of dfs
  expect_length(drivers, 9)
  expect_length(drivers$mean_stats, 10)
  expect_length(drivers$weighted_mean_stats, 10)
  
  #check that genes used for calculations overlap both datasets and loadings
  expect_true(unique(drivers$mean_stats$gene %in% rownames(microglial_counts)))
  expect_true(unique(drivers$mean_stats$gene %in% rownames(glial_counts)))
  expect_true(unique(drivers$mean_stats$gene %in% rownames(retinal_patterns)))
  expect_true(unique(drivers$weighted_mean_stats$gene %in% rownames(microglial_counts)))
  expect_true(unique(drivers$weighted_mean_stats$gene %in% rownames(glial_counts)))
  expect_true(unique(drivers$weighted_mean_stats$gene %in% rownames(retinal_patterns)))
  
  #name and class checks
  expect_true("mean_stats" %in% names(drivers))
  expect_is(drivers$mean_stats, "data.frame")
  
  expect_true("weighted_mean_stats" %in% names(drivers))
  expect_is(drivers$mean_stats, "data.frame")
  
  expect_true("normalized_weights" %in% names(drivers))
  expect_is(drivers$normalized_weights, "numeric")
  
  expect_true("sig_genes" %in% names(drivers))
  expect_is(drivers$sig_genes, "list")
  expect_length(drivers$sig_genes, 3)
  
  expect_true(unique(c("PV_sig", "weighted_PV_sig", "PV_significant_shared_genes") %in% names(drivers$sig_genes)))
  
  expect_is(drivers$sig_genes$PV_sig, "character")
  
  expect_is(drivers$sig_genes$weighted_PV_sig, "character")
  
  expect_is(drivers$sig_genes$PV_significant_shared_genes, "character")
  
  expect_true("meta_data" %in% names(drivers))
  expect_is(drivers$meta_data, "list")
  expect_length(drivers$meta_data, 3)
  expect_true("pvalue" %in% names(drivers$meta_data))
  expect_is(drivers$meta_data$pvalue, "numeric")
  
  #check that matrix names are proper and match source names
  expect_true(deparse(substitute(microglial_counts)) == drivers$meta_data$test_matrix)
  expect_type(drivers$meta_data$test_matrix, "character")
  
  expect_true(deparse(substitute(glial_counts)) == drivers$meta_data$reference_matrix)
  expect_type(drivers$meta_data$reference_matrix, "character")
  
})


