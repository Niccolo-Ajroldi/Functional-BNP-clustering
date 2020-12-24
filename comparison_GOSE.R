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


#### CLUSTER-GOSE COMPARISON ####

# GOSE
GOSE <- data$GOSE

# remove the observations that I removed from X
GOSE.rm <- GOSE[-c(12,13,19,24)]

GOSE.rm
partition

x11()
matplot(t(X), type='l', col=GOSE.rm, main="Curves clustered by GOSE")
legend("bottomright",
       legend=c("GOSE=1 (bene)", "GOSE=2 (male)"),
       col=c(1,2), 
       lty=c(1,1),
       lwd=c(2,2),
       cex=0.75)

x11()
matplot(t(X), type='l', col=partition, main="Curves clustered by FBNP")

# get the group with highest numerosity in our clustering 
library(plyr)
biggest.group <- which( count(partition)$freq == max(count(partition)) )
biggest.group

x11()
matplot(t(X), type='l', col=1+(partition!=biggest.group), main="Curves cluster by FBNP, minor clusters merged")
legend("bottomright",
       legend=c("Curves in largest cluster", "Curves in all other clusters"),
       col=c(1,2), 
       lty=c(1,1),
       lwd=c(2,2),
       cex=0.75)

# patients with GOSE==2 che riusciamo a separare dagli altri
which(partition!=biggest.group)[which(partition!=biggest.group) %in% which(GOSE.rm==2)]

# percentage of GOSE==2 patients identified
length(which(which(partition!=biggest.group) %in% which(GOSE.rm==2)))/length(which(GOSE.rm==2))


