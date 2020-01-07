
rm(list =ls())

setwd("C:/Users/ANeil/Google Drive/AA_PhD_TAMU/Hybrid Wheat Project/AE_Genomic Selection/GS_Analysis")

list.files()

geno <- read.csv("Geno_data.csv", header = FALSE, stringsAsFactors = TRUE)
pheno <- read.csv("AE_pheno.csv",header=TRUE, stringsAsFactors = TRUE)

GID_pheno <- pheno$GENO_ID2

# transpose the Geno file, markers in cols and Genotypes in row  
M1 <- t(geno)
colnames(M1) <- M1[1,]
M1 <- M1[-1,]

# convert the M1 into a dataframe
df <- data.frame(M1)

# Remove duplicate genotypes based on GID
library(data.table)
df2 <- unique(setDT(df), by = "GID")

# Extract marker data for lines in phenotypic dataset
df3 <- subset(df2, df2$GID %in% GID_pheno)

# write the final dataset to a csv file
write.csv(df3,"GS_AE_Final_Genotype_Data.csv")



#-----------------------------------------------------------------------------------------------------
#               RECODE GENOTYPE DATA TO MATCH BGLR REQUIREMENTS
#-----------------------------------------------------------------------------------------------------

rm(list =ls())
setwd("C:/Users/ANeil/Google Drive/AA_PhD_TAMU/Hybrid Wheat Project/AE_Genomic Selection/GS_Analysis")

Genotype_info=read.csv(file="GS_AE_Final_Genotype_Data.csv",
                       header=TRUE,na.strings="-/-",stringsAsFactors=FALSE)

#get rid of "/"
geno <- as.data.frame(sapply(Genotype_info[-1,-1], gsub, pattern = "/", replacement = "", fixed = TRUE))

t <- as.data.frame(t(geno))

colnames(t) <- Genotype_info$GID[2:606]
t<-t[-1,]

# create a new columns for counts of A T C G and missing
t$AA <- apply(t, 1, function(x) length(which(x=="AA")))
t$CC <- apply(t, 1, function(x) length(which(x=="CC")))
t$GG <- apply(t, 1, function(x) length(which(x=="GG")))
t$TT <- apply(t, 1, function(x) length(which(x=="TT")))
t$missing <- apply(t, 1, function(x) length(which(x=="NA")))


n <- ncol(t)

# find the major allele
for(i in 1:nrow(t)){
  t[i,(n+1)] <- colnames(t[which.max(t[i,])])
}

# find the minor allele
## sort the select columns and give the column name of the second highest value
for(i in 1:nrow(t)){
  t[i,(n+2)] <- colnames(
    sort(t[i,(n-4):(n-1)],TRUE)[2]
    )
}

# create an allele id based on major/minor 

for(i in 1:nrow(t)){
  t[i,(n+3)] <- paste(substring(t[i,n+1],2),
                      substring(t[i,n+2],2),
                      sep = "/")
}


# lets add the allele call column major/minor to the begining of the dataframe
library(tibble)
t <- add_column(t, alleles = t$V613, .after = 0)

write.csv(t,"tmp.csv")

tmp <- t
tmp <- add_column(tmp, Major = tmp$V611, .after = 1)
tmp <- add_column(tmp, Minor = tmp$V612, .after = 2)
tmp <- tmp[,1:609]

# add other columns as needed for hapmap format
# other columns and arbritary chromosome and position info were added in excel outside 
# file is saved as GS_hapmap.csv

#---------------------------------------------------------------------------------------
#           Recode to numeric format 
#---------------------------------------------------------------------------------------

tmp.mat <- as.matrix(tmp)
tmp.mat[is.na(tmp.mat)] <- 0 

for(i in 1:nrow(tmp.mat)){
  for(j in 4:ncol(tmp.mat)){ 
    if(tmp.mat[i,j] == tmp.mat[i,3]){
     tmp.mat[i,j] <- '-1'}
  }
}


tmp2 <-tmp.mat

# recode hetrozygote to 0
for(i in 1:nrow(tmp.mat)){
  for(j in 4:ncol(tmp.mat)){ 
    if(tmp.mat[i,j] != c(1,0,-1)){
      tmp.mat[i,j] <- '0'}
  }
}

View(tmp.mat[1:20,1:20])
write.csv(tmp.mat,"Geno_recoded.csv")

