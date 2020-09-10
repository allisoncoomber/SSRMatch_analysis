library(dplyr)

#Create a list
confusionMatrices <- list()

# Read in the confusion matrices into a list of dataframes
for(i in 1:5) {
  confusionMatrices[[i]] <- read.table(paste("/Users/Allison/git/SSRmatch_analysis/K1_sample", i, "_matrix.txt", sep=""), header=TRUE, encoding="UTF-8", sep = ",")
}

#
bind_rows(confusionMatrices, .id = "column_label") %>%
  group_by(X) %>%
  summarise(across(c(BR_1.:US_17.), sum)) %>%
  write.csv(paste("/Users/Allison/git/SSRmatch_analysis/K1_summary_confusion_matrix.csv", sep=","), fileEncoding="UTF-8")
