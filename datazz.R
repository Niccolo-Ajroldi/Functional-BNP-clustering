
setwd("D:/Poli/Corsi/BAYESIAN/Proj/dati")

#### IMPORTAZIONE DATI ####
data <- read.table('dati.csv', head=TRUE, sep=',', row.names = "id")
data <- data[,-c(5,9,21,22,23,24,25,26)]
data[,2] <- ifelse(data[,2] == 1, "M", "F")
data[,3] <- ifelse(data[,3] == 1, "TR", "VS")
data[,4] <- ifelse(data[,4] == 1, "DX", "SX")
data[,7] <- ifelse(data[,7] == 'FRMN', "SI", "NO")
data[,11] <- ifelse(data[,11] == 'B', 2, 1)
data[,8] <- as.factor(toupper(data[,8]))
data[,9] <- as.factor(toupper(data[,9]))
data[,12] <- as.factor(toupper(data[,12]))
data[,13] <- as.factor(toupper(data[,13]))
data[data[,12]=='NO ',12] <- c('NO')
data[data[,13]=='NO ',13] <- c('NO')
data[,12] <- factor(data[,12])
data[,13] <- factor(data[,13])  

colnames(data) <- c("AGE", "SEX", "EZIOLOGIA", "LATE", "MESIC", "MESIR", "FRMN",
                    "SLDX", "SLSX", "SLSEP", "SLSEPINC", "MLDX", "MLSX", 
                    "MLSEP", "MLSEP2", "GOSE", "LCF", "DRS")


#' GOSE: 1: male
#'       2: buono
#' 
#' LCF: Level of Cognitive Functioning
#'      1: sfavorevole
#'      2: bene
#'      
#' DRS: Disability Rate Scale (menomazione)
#'      1: male
#'      2: bene
#'      

levels(as.factor(data$GOSE))
levels(as.factor(data$LCF))
levels(as.factor(data$DRS))


#### CLUSTER-GOSE COMPARISON ####

# GOSE
GOSE <- data$GOSE
LCF  <- data$LCF
DRS  <- data$DRS

# select a partition
partition <- minbinder.ext(psm,cls.draw = K, method="comp")[[1]]
partition

x11()
matplot(t(X), type="l", col=partition)

# partition by clusters + GOSE
x11()
par(mfrow = n2mfrow(length(unique(partition))))
count <- 1
for(j in levels(as.factor(partition)))
{
  indexes.j <- which(partition==j)
  if(length(indexes.j)==1)
    matplot(X[indexes.j,], type='l', col=GOSE[indexes.j], xlab="", main=paste0("Cluster ",count))
  else
    matplot(t(X[indexes.j,]), type='l', col=GOSE[indexes.j], xlab="", main=paste0("Cluster ",count))
  count <- count+1
}

# partition by clusters + LCF
x11()
par(mfrow = n2mfrow(length(unique(partition))))
count <- 1
for(j in levels(as.factor(partition)))
{
  indexes.j <- which(partition==j)
  if(length(indexes.j)==1)
    matplot(X[indexes.j,], type='l', col=LCF[indexes.j], xlab="", main=paste0("Cluster ",count))
  else
    matplot(t(X[indexes.j,]), type='l', col=LCF[indexes.j], xlab="", main=paste0("Cluster ",count))
  count <- count+1
}

# partition by clusters + DRS
x11()
par(mfrow = n2mfrow(length(unique(partition))))
count <- 1
for(j in levels(as.factor(partition)))
{
  indexes.j <- which(partition==j)
  if(length(indexes.j)==1)
    matplot(X[indexes.j,], type='l', col=DRS[indexes.j], xlab="", main=paste0("Cluster ",count))
  else
    matplot(t(X[indexes.j,]), type='l', col=DRS[indexes.j], xlab="", main=paste0("Cluster ",count))
  count <- count+1
}

