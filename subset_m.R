library(dplyr)
library(readxl)
library(tidyverse)
library(writexl)
library(glmnet)

set.seed(5)

path = "C:/Users/Josep/Documents/Work/Research/HP Study - Journey's End/Data/"

m_mod <- readRDS(paste(path, "Processed/m_mod.rda", sep=""))

fit_df <- read_csv(paste(path, "Processed/fit_m.csv", sep=""), show_col_types = FALSE)
test_df <- read_csv(paste(path, "Processed/predict_m.csv", sep=""), show_col_types = FALSE)

fit_df = select_if(fit_df, is.numeric)
fit_df %>% select_if(colSums(.) != 0)

test_df = select_if(test_df, is.numeric)
test_df %>% select_if(colSums(.) != 0)

for (i in colnames(fit_df)) {
	if (!(i %in% colnames(test_df))) {
	newcol = c(rep(0, nrow(test_df)))
	test_df <- cbind(newcol, test_df)
	colnames(test_df)[1] <- i
		}}

test_df<-test_df[names(fit_df)]

test_x <- as.matrix(test_df[,1:ncol(test_df)-1])
test_y <- as.matrix(test_df[,ncol(test_df)])

predictions <- predict(m_mod[[1]], newx=test_x, s="lambda.min", type="class")

prediction_df <- data.frame(cbind(predictions))

colnames(prediction_df) <- c('prediction')

write_csv(prediction_df, paste(path, "Processed/predicted_m.csv", sep=""))
