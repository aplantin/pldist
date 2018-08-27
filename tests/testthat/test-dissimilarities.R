context("all dissimilarities")

test_that("dissimilarities give expected results", {
  otus <- matrix(nrow = 4, ncol = 3) 
  rownames(otus) <- paste(rep(paste("subj", 1:2, sep = ""), each = 2), 
                          rep(c("a","b"), 2), sep = "")
  metadata <- data.frame(subjID = rep(paste("subj", 1:2, sep = ""), each = 2), 
                         sampID = paste(rep(paste("subj", 1:2, sep = ""), each = 2), 
                                        rep(c("a","b"), 2), sep = ""), 
                         time = rep(1:2, 2)) 
  otus[1, ] <- c(0, 0.2, 0.8)
  otus[2, ] <- c(0.4, 0, 0.6)
  otus[3, ] <- c(0, 0.2, 0.8)
  otus[4, ] <- c(0.4, 0, 0.6)
  colnames(otus) <- paste("otu", 1:3, sep = "")
  sim.tree <- rtree(3, tip.label = paste("otu", 1:3, sep = ""))
  
  ## Paired 
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "bray")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "bray")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "jac")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "jac")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "kul")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "kul")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "gow")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "gow")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "unifrac", tree = sim.tree)$D[1,2,"d_0"], 0)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "unifrac", tree = sim.tree)$D[1,2,"d_1"], 0)
  
  ## Not paired 
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "bray")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "bray")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "jac")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "jac")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "kul")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "kul")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "gow")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "gow")$D[1,2], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "unifrac", tree = sim.tree)$D[1,2,"d_0"], 0)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "unifrac", tree = sim.tree)$D[1,2,"d_1"], 0)
})