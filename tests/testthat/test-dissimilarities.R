context("all dissimilarities")

test_that("dissimilarities give expected results", {
  otus <- matrix(nrow = 6, ncol = 3) 
  rownames(otus) <- paste(rep(paste("subj", 1:3, sep = ""), each = 2), 
                          rep(c("a","b"), 3), sep = "")
  metadata <- data.frame(subjID = rep(paste("subj", 1:3, sep = ""), each = 2), 
                         sampID = paste(rep(paste("subj", 1:3, sep = ""), each = 2), 
                                        rep(c("a","b"), 3), sep = ""), 
                         time = rep(1:2, 3)) 
  otus[1, ] <- c(0, 0.2, 0.8)
  otus[2, ] <- c(0.1, 0.3, 0.6)
  otus[3, ] <- c(0.4, 0.4, 0.2) 
  otus[4, ] <- c(0.2, 0.8, 0) 
  otus[5, ] <- c(0.2, 0, 0.8) 
  otus[6, ] <- c(0, 0.4, 0.6) 
  colnames(otus) <- paste("otu", 1:3, sep = "")
  set.seed(5); sim.tree <- rtree(3, tip.label = paste("otu", 1:3, sep = ""))
  
  paired.dat <- pltransform(otus, metadata, paired = TRUE)
  longit.dat <- pltransform(otus, metadata, paired = FALSE) 
  
  ## Bray-Curtis 
  bray.pb <- 1/3
  bray.pq <- sum(abs(c(0.5+1/6, 0.1-1/6, -1/14+1/2)))/3 
  bray.lb <- 2/3
  bray.lq <- sum(abs(c(2/3, 1/3-0.2, 6/7)))/3

  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "bray")$D[1,2], bray.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "Bray")$D[1,2], bray.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "br")$D[1,2], bray.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "bray")$D[1,2], bray.lq)
  
  ## Jaccard 
  jac.pb <- 1
  jac.pq <- 1 - sum(0, 0.1, 1/14)/sum(0.5, 1/6, 0.5)
  jac.lb <- 1 
  jac.lq <- 1 - (1/3 + 0.2 + 1/7)/(1 + 1/3 + 1)
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "jaccard")$D[1,2], jac.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "jac")$D[1,2], jac.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "j")$D[1,2], jac.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "Jacc")$D[1,2], jac.lq)
  
  ## Kulczynski 
  kul.pb <- 1 - 1/3
  kul.pq <- 1 - 0.5 * (1/(0.5+0.1+1/14) + 1/(1/6+1/6+0.5)) * (0.1 + 1/14)
  kul.lb <- 1 
  kul.lq <- 1 - 0.5 * (1/(1 + 0.2 + 1/7) + 1/(1 + 2/3)) * sum(c(1/3, 0.2, 1/7))
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "kul")$D[1,2], kul.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "kul")$D[1,2], kul.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "kul")$D[1,2], kul.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "kul")$D[1,2], kul.lq)
  
  
  ## Gower 
  gow.pb <- 1/3 * (0.5/1 + 0 + 0.5/0.5) 
  gow.pq <- 1/3 * ((2/3)/1 + (1/6-0.1)/0.4 + (0.5+1/14)/(0.5+1/14)) 
  gow.lb <- 1/3 * (1 + 0 + 1)
  gow.lq <- 1/3 * (1 + (1/3-0.2) / 0.8 + 1)
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "gower")$D[1,2], gow.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "Gower")$D[1,2], gow.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "gow")$D[1,2], gow.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "Gow")$D[1,2], gow.lq)
  
  
  ## UniFrac
  otu1.12 <- sim.tree$edge.length[2] 
  otu2.12 <- sim.tree$edge.length[3] 
  otu12.123 <- sim.tree$edge.length[1]
  otu3.123 <- sim.tree$edge.length[4] 
  uf.pb <- (otu1.12*0.5 + otu3.123*0.5) / sum(sim.tree$edge.length)
  uf.pq <- (otu1.12*0.35*2/3 + otu2.12*0.85*(1/6-0.1) + otu3.123*0.8*3/7 + otu12.123*1.2*1/9) / 
    sum(sim.tree$edge.length * c(1.2, 0.35, 0.85, 0.8))
  uf.lb <- (otu1.12 + otu3.123)/sum(sim.tree$edge.length)
  uf.lq <- (otu1.12*0.35*2/3 + otu2.12*0.85*(1/3-0.2) + otu3.123*0.8*6/7 + otu12.123*1.2*2/9) / 
    sum(sim.tree$edge.length * c(1.2, 0.35, 0.85, 0.8))
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "unifrac", tree = sim.tree)$D[1,2,"d_UW"], uf.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "unifrac", tree = sim.tree)$D[1,2,"d_1"], uf.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "unifrac", tree = sim.tree)$D[1,2,"d_UW"], uf.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "unifrac", tree = sim.tree)$D[1,2,"d_1"], uf.lq)
  
})
