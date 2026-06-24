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
	withr::local_dir(withr::local_tempdir())
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
	expect_true(file.exists("multivariateAnalysisR_ANOVA.png"))
	expect_true(file.exists("multivariateAnalysisR_ANOVA.csv"))
	expect_true(file.exists("multivariateAnalysisR_CI.png"))
	expect_true(file.exists("multivariateAnalysisR_CI.csv"))
	
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

test_that("projection works on sparse data matrix", {
  dense <- as.matrix(p.RNAseq6l3c3t)
  sparse <- as(p.RNAseq6l3c3t, "sparseMatrix")
  loadings <- CR.RNAseq6l3c3t@featureLoadings

  expect_true("dgCMatrix" %in% class(sparse))
  expect_no_error(projectR(sparse, loadings))

  pdense <- projectR(dense, loadings)
  psparse <- projectR(sparse, loadings)
  expect_identical(pdense, psparse)
})

test_that("projection works on sparse data matrix with full=TRUE", {
  dense <- as.matrix(p.RNAseq6l3c3t)
  sparse <- as(p.RNAseq6l3c3t, "sparseMatrix")
  loadings <- CR.RNAseq6l3c3t@featureLoadings

  expect_true("dgCMatrix" %in% class(sparse))
  expect_no_error(projectR(sparse, loadings))

  pdense <- projectR(dense, loadings, full=TRUE)
  #case with default number of chopBy, when it's > ncol
  psparse <- projectR(sparse, loadings, full=TRUE)
  expect_identical(pdense, psparse)

  #case with chopBy < ncol
  psparse_chunked <- projectR(sparse, loadings, full=TRUE, chopBy=10)
  expect_identical(pdense, psparse_chunked)
})

test_that("prcomp dispatch full=TRUE returns named list with correct slots", {
  pca <- prcomp(t(p.RNAseq6l3c3t), center = TRUE)
  result_full <- projectR(
    data = p.ESepiGen4c1l$mRNA.Seq, loadings = pca,
    dataNames = map.ESepiGen4c1l[["GeneSymbols"]], full = TRUE
  )
  expect_type(result_full, "list")
  expect_named(result_full, c("projection", "pvar", "r_squared"))
  expect_equal(
    dim(result_full$projection),
    dim(projectR(p.ESepiGen4c1l$mRNA.Seq, pca,
      dataNames = map.ESepiGen4c1l[["GeneSymbols"]]
    ))
  )
  expect_true(all(result_full$r_squared <= 1))
  expect_named(result_full$r_squared)
})

test_that("matrix dispatch full=TRUE returns named list with correct slots", {
  result_full <- projectR(
    data = p.ESepiGen4c1l$mRNA.Seq, loadings = AP.RNAseq6l3c3t$Amean,
    dataNames = map.ESepiGen4c1l[["GeneSymbols"]], full = TRUE
  )
  expect_type(result_full, "list")
  expect_named(result_full, c("projection", "pval", "r_squared"))
  expect_equal(
    dim(result_full$projection),
    dim(projectR(p.ESepiGen4c1l$mRNA.Seq, AP.RNAseq6l3c3t$Amean,
      dataNames = map.ESepiGen4c1l[["GeneSymbols"]]
    ))
  )
  expect_true(all(result_full$pval >= 0 & result_full$pval <= 1))
  expect_true(all(result_full$r_squared <= 1))
  expect_named(result_full$r_squared)
})

test_that("center_by_loadings produces different projection than per-sample centering", {
  pca <- prcomp(t(p.RNAseq6l3c3t), center = TRUE)
  proj_default <- projectR(
    data = p.ESepiGen4c1l$mRNA.Seq, loadings = pca,
    dataNames = map.ESepiGen4c1l[["GeneSymbols"]], center_by_loadings = FALSE
  )
  proj_centered <- projectR(
    data = p.ESepiGen4c1l$mRNA.Seq, loadings = pca,
    dataNames = map.ESepiGen4c1l[["GeneSymbols"]], center_by_loadings = TRUE
  )
  expect_false(identical(proj_default, proj_centered))
})

test_that("center_by_loadings=TRUE errors when loadings$center is absent", {
  pca_no_center <- prcomp(t(p.RNAseq6l3c3t), center = FALSE)
  expect_error(
    projectR(
      data = p.ESepiGen4c1l$mRNA.Seq, loadings = pca_no_center,
      dataNames = map.ESepiGen4c1l[["GeneSymbols"]], center_by_loadings = TRUE
    ),
    regexp = "loadings\\$center is missing or FALSE"
  )
})

test_that("matrix dispatch include_intercept=FALSE returns no intercept slot", {
  result <- projectR(
    data      = p.ESepiGen4c1l$mRNA.Seq,
    loadings  = AP.RNAseq6l3c3t$Amean,
    dataNames = map.ESepiGen4c1l[["GeneSymbols"]],
    full      = TRUE
  )
  expect_null(result$intercept)
})

test_that("matrix dispatch include_intercept=TRUE adds intercept slot and raises R2", {
  r0 <- projectR(
    data              = p.ESepiGen4c1l$mRNA.Seq,
    loadings          = AP.RNAseq6l3c3t$Amean,
    dataNames         = map.ESepiGen4c1l[["GeneSymbols"]],
    full              = TRUE
  )
  r1 <- projectR(
    data              = p.ESepiGen4c1l$mRNA.Seq,
    loadings          = AP.RNAseq6l3c3t$Amean,
    dataNames         = map.ESepiGen4c1l[["GeneSymbols"]],
    full              = TRUE,
    include_intercept = TRUE
  )
  expect_length(r1$intercept, ncol(p.ESepiGen4c1l$mRNA.Seq))
  expect_named(r1$intercept)
  expect_equal(dim(r0$projection), dim(r1$projection))
  expect_equal(dim(r0$pval), dim(r1$pval))
  expect_true(all(r1$r_squared >= r0$r_squared))
})

test_that("prcomp dispatch with scalar NP does not drop matrix dimensions", {
  pca <- prcomp(t(p.RNAseq6l3c3t), center = TRUE)
  result <- projectR(
    data               = p.ESepiGen4c1l$mRNA.Seq,
    loadings           = pca,
    dataNames          = map.ESepiGen4c1l[["GeneSymbols"]],
    NP                 = 3L,
    center_by_loadings = TRUE,
    full               = TRUE
  )
  expect_type(result, "list")
  expect_equal(nrow(result$projection), 1L)
  expect_equal(ncol(result$projection), ncol(p.ESepiGen4c1l$mRNA.Seq))
})