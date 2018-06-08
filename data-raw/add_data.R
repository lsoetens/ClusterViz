example <- read.table("data_example.txt", stringsAsFactors = F, header=T, sep=",") 

devtools::use_data(example)