suppressWarnings({

	library(dplyr)
	library(readxl)
	library(tidyverse)
	library(writexl)
	library(glmnet)

	path = "C:/Users/Josep/Documents/Work/Research/HP Study - Journey's End/Data/"

	fit_df <- read_csv(paste(path, "Processed/fit_am.csv", sep=""), show_col_types = FALSE)
	test_df <- read_csv(paste(path, "Processed/test_am.csv", sep=""), show_col_types = FALSE)

	basis = unique(fit_df$basis)

	error_rate = c()
	
	formulas <- vector(mode="list", length=length(basis))

	for (b in basis){
	set.seed(5)
	basis_b_df = subset(fit_df, basis == b)
	basis_b_df = select_if(basis_b_df, is.numeric)
	basis_b_df %>% select_if(colSums(.) != 0)

	basis_test_df = subset(test_df, basis == b)
	basis_test_df = select_if(basis_test_df, is.numeric)
	basis_test_df %>% select_if(colSums(.) != 0)
	
	for (i in colnames(basis_b_df)) {
		if (!(i %in% colnames(basis_test_df))) {
		newcol = c(rep(0, nrow(basis_test_df)))
		basis_test_df <- cbind(newcol, basis_test_df)
		colnames(basis_test_df)[1] <- i
			}}
  
	basis_test_df<-basis_test_df[names(basis_b_df)]
	
	x <- as.matrix(basis_b_df[,1:ncol(basis_b_df)-1])
	y <- as.matrix(basis_b_df[,ncol(basis_b_df)])
	
	test_x <- as.matrix(basis_test_df[,1:ncol(basis_test_df)-1])
	test_y <- as.matrix(basis_test_df[,ncol(basis_test_df)])

	cvfit <- cv.glmnet(x, y, family = 'binomial', type.measure = "class")
	
	formulas[[match(c(b),basis)]] <- cvfit
	
	predictions <- predict(cvfit, newx=test_x, s="lambda.min", type="class")
	
	z = round(sum(abs(test_y - as.numeric(unlist(predictions))))/length(test_y),4)
	
	error_rate = c(error_rate, z)
	
  print(b)
  print(z*100)
    
	}

	min_error <- min(error_rate)
	best_basis <- basis[which.min(error_rate)]
	best_basis_mod <- formulas[which.min(error_rate)]
	
	error_df = data.frame(cbind(best_basis, min_error))
	colnames(error_df) <- c('basis', 'error')
	
  write_csv(error_df, paste(path, "Processed/am_basis.csv", sep=""))
  saveRDS(best_basis_mod, file = paste(path, "Processed/am_mod.rda", sep=""))
  
})
