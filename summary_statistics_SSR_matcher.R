library(dplyr)
library(caret)

#Create a list
confusionMatrices <- list()
genotypePredictionResults <- list()

# Read in the confusion matrices into a list of dataframes
for(i in 1:5) {
  confusionMatrices[[i]] <- read.table(paste("/Users/Allison/git/SSRmatch_analysis/K1_sample", i, "_matrix.txt", sep=""), header=TRUE, encoding="UTF-8", sep = ",")
}

#
bind_rows(confusionMatrices, .id = "column_label") %>%
  group_by(X) %>%
  summarise(across(c(BR_1.:US_17.), sum)) %>%
  write.csv(paste("/Users/Allison/git/SSRmatch_analysis/K1_summary_confusion_matrix.csv", sep=","), fileEncoding="UTF-8")

# Read in the prediction results into a list of dataframes
for(i in 1:5) {
  genotypePredictionResults[[i]] <- read.table(paste("/Users/Allison/git/SSRmatch_analysis/K1_sample", i, "_results.txt", sep=""), header=TRUE, encoding="UTF-8", sep = "\t")
}

# Bind all of the results into one big table
genotypePredictionTable <- bind_rows(genotypePredictionResults)

# Coerce the columns to factors and give them new names
genotypePredictionTable <- rename(genotypePredictionTable, c("predictedGenotype"="MatchingGenotype.", "actualGenotype"="Genotype."))
genotypePredictionTable$predictedGenotype <- as.factor(genotypePredictionTable$predictedGenotype)
genotypePredictionTable$actualGenotype <- as.factor(genotypePredictionTable$actualGenotype)

summaryStatistics <- confusionMatrix(genotypePredictionTable$predictedGenotype, genotypePredictionTable$actualGenotype)