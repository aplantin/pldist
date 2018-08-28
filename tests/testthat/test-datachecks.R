context("data checking")

test_that("function stops when it should", {
  data("paired.meta")
  data("paired.otus") 
  
  # Unknown method
  expect_error( pldist(otus = paired.otus, metadata = paired.meta, paired = FALSE, binary = FALSE, method = "clark", tree = NULL, gam = c(0, 0.5, 1) ), "Method does not match" )
  
  # Different # observations in OTU table and metadata 
  expect_error( pldist(otus = rbind(paired.otus, c(1:10)), metadata = paired.meta, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "Number of rows" )
  
  # Sample with all zero values  
  ext.otus <- rbind(paired.otus, rep(0, 10), rep(1:10))
  rownames(ext.otus)[11:12] <- c("subj6a", "subj6b")
  ext.meta <- rbind(paired.meta, c("subj6", "subj6a", 1), c("subj6", "subj6b", 2))
  expect_error( pldist(otus = ext.otus, metadata = ext.meta, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "uniformly zero OTU counts" )
  
  # OTU with all zero values 
  expect_warning( pldist(otus = cbind(paired.otus, rep(0, 10)), metadata = paired.meta, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "OTUs have count zero" )
  
  # Metadata is out of order 
  expect_error(pldist(otus = paired.otus, metadata = paired.meta[,c(2,1,3)], paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "format metadata" )
  
  # OTU matrix rownames don't match sample IDs 
  ext.otus <- paired.otus 
  rownames(ext.otus) <- paste("samp", 1:10, sep = "")
  expect_error(pldist(otus = ext.otus, metadata = paired.meta, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "rownames of OTU matrix")
  
  # Asked for paired dissimilarity on unpaired data 
  data("bal.long.meta")
  data("bal.long.otus")
  expect_error(pldist(otus = bal.long.otus, metadata = bal.long.meta, paired = TRUE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), ">2 unique time points/groups were provided")
  
  # Asked for paired dissimilarity, but missing subject 
  expect_error(pldist(otus = paired.otus[-10,], metadata = paired.meta[-10,], paired = TRUE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "some groups/subjects do not have 2 observations")
  
  # Asked for longitudinal dissimilarity on unbalanced data 
  data("unbal.long.meta")
  data("unbal.long.otus")
  expect_warning(pldist(otus = unbal.long.otus, metadata = unbal.long.meta, paired = FALSE, binary = FALSE, method = "bray", tree = NULL, gam = c(0, 0.5, 1) ), "unbalanced designs")
})

