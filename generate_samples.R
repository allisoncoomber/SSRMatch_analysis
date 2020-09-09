library (dplyr)

all_data <- read.table("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/all_12_defined.txt", header=TRUE, encoding="UTF-8",stringsAsFactors = FALSE, colClasses = c(rep("character",20)))

set.sample.size<-(nrow(all_data))/5
num.cycles<-5
all_data$id <- 1:nrow(all_data)
View(all_data)
for(i in 1:(num.cycles)) {
  nam <- paste("sample_", i, sep = "")
  assign(nam, all_data[sample(nrow(all_data), set.sample.size), ])
  write.csv(get(nam), paste("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/sample_", i, ".csv", sep=""))
  all_data<-all_data[!(all_data$id %in% get(nam)$id),]
  }

all_data <- read.table("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/all_12_defined.txt", header=TRUE, encoding="UTF-8",stringsAsFactors = FALSE, colClasses = c(rep("character",20)))
all_data$id <- 1:nrow(all_data)
for(i in 1:5) {
  write.csv(anti_join(all_data, get(paste("sample_", i, sep = "")), by = "Sample"), paste("/Users/Allison/Documents/NCSU/Ristaino Lab/Workshop/sample_", i, "_reference.csv", sep =""))
}
