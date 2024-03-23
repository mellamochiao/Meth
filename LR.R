library(MEAL)
library(GEOquery)
library(minfi)
library(limma)
library(ggplot2)

gse30870 <- getGEO("GSE30870", GSEMatrix = TRUE, AnnotGPL = TRUE)


#gse40279_matrix <- gse40279[[1]]
class(gse30870)

gse30870_matrix <- gse30870[[1]]
data <- exprs(gse30870_matrix)
data <- as.data.frame(data)
#把DATA TRANSPOSE
tdata<-t(data)




summary(exprs(gse30870_matrix))


age <- pData(gse30870_matrix)$characteristics_ch1
ages_processed <- gsub("age:", "", age)
ages_processed <- gsub(" years", "", ages_processed)
ages_processed <- gsub("Newborn", "0", ages_processed)
ages_processed <- as.numeric(ages_processed)
agevector <- c(ages_processed)

class(age)
class(tdata)

agevector <- t(agevector)
agevector <- as.vector(agevector)
#先轉成dataframe然後合併年齡跟DATA
tdata <- as.data.frame(tdata)
tdata$age <- agevector

data<-t(tdata)

class(tdata)


require(datasets)  

model <- lm(age ~ cg25088412, data = tdata)
summary(model)
model <- lm(age ~ cg13485721, data = tdata)
summary(model)
model <- lm(age ~ cg15126211, data = tdata)
summary(model)

