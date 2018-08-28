context("transformations")

test_that("transformations give expected result", {
  otus <- matrix(nrow = 4, ncol = 3) 
  rownames(otus) <- paste(rep(paste("subj", 1:2, sep = ""), each = 2), 
                          rep(c("a","b"), 2), sep = "")
  metadata <- data.frame(subjID = rep(paste("subj", 1:2, sep = ""), each = 2), 
                         sampID = paste(rep(paste("subj", 1:2, sep = ""), each = 2), 
                                        rep(c("a","b"), 2), sep = ""), 
                         time = rep(1:2, 2)) 
  otus[1, ] <- c(0, 0.2, 0.8)
  otus[2, ] <- c(0.4, 0, 0.6)
  otus[3, ] <- c(0.4, 0.4, 0.2) 
  otus[4, ] <- c(0.2, 0.8, 0) 
  
  ## Paired 
  paired.tsf <- pltransform(otus, metadata, paired = TRUE)
  
  exp.binary.subj1 <- c(0.5, -0.5, 0) 
  exp.quant.subj1 <- c(0.5, -0.5, -1/14)
  exp.avg.subj1 <- (otus[1,] + otus[2,])/2
  
  expect_equal(pltransform(otus, metadata, paired = TRUE)$tsf.data$dat.binary[1,], exp.binary.subj1)
  expect_equal(pltransform(otus, metadata, paired = TRUE)$tsf.data$dat.quant[1,], exp.quant.subj1)
  expect_equal(pltransform(otus, metadata, paired = TRUE)$tsf.data$avg.prop[1,], exp.avg.subj1) 
  
  ## Longitudinal (balanced, time points 1 and 2)
  long.tsf <- pltransform(otus, metadata, paired = FALSE) 
  
  exp.binary.subj1 <- c(1, 1, 0) 
  exp.quant.subj1 <- c(1, 1, 1/7)
  exp.avg.subj1 <- (otus[1,] + otus[2,])/2
  
  expect_equal(pltransform(otus, metadata, paired = FALSE)$tsf.data$dat.binary[1,], exp.binary.subj1)
  expect_equal(pltransform(otus, metadata, paired = FALSE)$tsf.data$dat.quant[1,], exp.quant.subj1)
  expect_equal(pltransform(otus, metadata, paired = FALSE)$tsf.data$avg.prop[1,], exp.avg.subj1) 
  
  ## Longitudinal (unbalanced -- time points 1, 5, and 7) 
  otus <- matrix(nrow = 6, ncol = 3) 
  rownames(otus) <- paste(rep(paste("subj", 1:2, sep = ""), each = 3), 
                          rep(c("a","b","c"), 2), sep = "")
  metadata <- data.frame(subjID = rep(paste("subj", 1:2, sep = ""), each = 3), 
                         sampID = paste(rep(paste("subj", 1:2, sep = ""), each = 3), 
                                        rep(c("a","b","c"), 2), sep = ""), 
                         time = rep(c(1,5,7), 2)) 
  otus[1, ] <- c(0, 0.2, 0.8)
  otus[2, ] <- c(0.4, 0, 0.6)
  otus[3, ] <- c(0.4, 0.4, 0.2) 
  otus[4, ] <- c(0.2, 0.8, 0) 
  otus[5, ] <- c(0.4, 0.6, 0) 
  otus[6, ] <- c(0.6, 0.2, 0.2) 
  
  long.tsf.2 <- pltransform(otus, metadata, paired = FALSE) 
  
  exp.binary.subj1 <- (c(1, 1, 0)/4 + c(0, 1, 0)/2)/2 
  exp.quant.subj1 <- (c(1, 1, 1/7)/4 + c(0, 1, 1/2)/2)/2
  exp.avg.subj1 <- (otus[1,] + otus[2,] + otus[3,])/3 
  
  expect_equal(pltransform(otus, metadata, paired = FALSE)$tsf.data$dat.binary[1,], exp.binary.subj1)
  expect_equal(pltransform(otus, metadata, paired = FALSE)$tsf.data$dat.quant[1,], exp.quant.subj1)
  expect_equal(pltransform(otus, metadata, paired = FALSE)$tsf.data$avg.prop[1,], exp.avg.subj1) 
})

