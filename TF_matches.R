library(dplyr)

for(i in 1:5) {
  results <- read.table(paste("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/K1_sample", i, "_results.txt", sep = ""), header=TRUE, encoding="UTF-8",stringsAsFactors = FALSE, colClasses = c(rep("character",3)))
  confusion_matrix_for_all_genotypes <- table(results$MatchingGenotype, results$Genotype)
  write.csv(confusion_matrix_for_all_genotypes, paste(paste("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/K1_sample", i, "_matrix.txt", sep = "")))
  results$match <- NULL
  results[results$MatchingGenotype == results$Genotype,"match"] <- "TRUE"
  results[results$MatchingGenotype != results$Genotype,"match"] <- "FALSE"
  write.csv(results, paste("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/K1_sample", i, "_matches.txt", sep = ""))
}
  