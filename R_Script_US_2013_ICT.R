######################################################################################################################################################
# THIS SCRIPT CALCULATES THE ECONOMIC LOSSES ASSOCIATED WITH CYBER-ATTACS IN THE US ICT SECTOR for 180 days
# IT IS BASED ON THE WIOD DATABASE
# IT CAN BE REPLICATED USING DATA ON OTHER REGIONS
# @Rokhaya DIEYE, September 2019
######################################################################################################################################################
# Load libraries
require(foreign)
library(MASS)
require(utils)
require(gdata)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(gridExtra)
library(xtable)
######################################################################################################################################################
# Load the I/O matrix for the US in 2013
wiot_US13 <- read.dta("wiot_US_2013.dta")
A_star1 <- as.matrix(wiot_US13[,7:61]) 

######################################################################
# STEP 1: CALCULATE THE INOPERABILITY VECTOR FOR THE WHOLE TIME PERIOD
######################################################################

# Use the S-aggregation technique to aggregate the "Telecommunications" and the "Computer programming, consultancy and related activities; information service activities" sector to form the ICT sector
S <- matrix(0,55,55)
diag(S) <- 1
S[39,40] <- 1
S <- S[-40,]
A_star <- S%*%A_star1%*%t(S)
industr_name <- wiot_US13$IndustryDescription
industr_code <- wiot_US13$IndustryCode
industr_name <- industr_name[-40]
industr_code <- industr_code[-40]

# Rename the newly formed sector and record its position on the technical coefficients matrix A_star (39)
industr_name[39] <- "ICT"
g <- 39
dep <- A_star[g,]

# Define the vector RecTT of possible time horizon (in days)
RecTT <- c(10,30, 60, 90, 180, 365)

# Define the sector resiliebce coeffuciebt matrix KK
KK <-  matrix(0, nrow(A_star), ncol(A_star)) 

# Define the vector of initial inoperability levels
qq_0 <- seq(0.1,0.4,0.1)

# Start of the loop in order to assess the inoperability levels during the time span of the attack (RecTT)
for(t in 1:length(RecTT)){
  RecT <- RecTT[t]
  Q <- matrix(0,ncol(A_star),RecT)
  avg_q <- matrix(0,54,length(qq_0))
  med_q <- matrix(0,54,length(qq_0))
  
  for(k in 1:length(qq_0)){
    q_00 <- qq_0[k]
    Q_0 <- matrix(q_00/dep[g],ncol(A_star),1)
    K <- KK
    Q[RecT] <- 0.01
    
    for(j in 1:nrow(A_star)){
      diag(K)[j] <- log(q_00/Q[RecT])/(RecT*(1-diag(A_star)[j])) 
    }
    
    i <- 1
    Q[,1] <- Q_0*dep
    
    while(i<RecT-1){
      Q[,i+1] <- Q[,i] + K%*%(A_star%*%Q[,i] - Q[,i])
      i <- i+1
    }
    nam <- paste("Q",q_00,RecT,sep="_")
    assign(nam,Q)
    
    # Determining the average and median inoperability vector by sector over the period
    avg_q[,k] <- rowMeans(Q)
    med_q[,k] <- rowMedians(Q)
    
  } 
  
  nam1 <- paste("avg_q",RecT, sep=".")
  nam2 <- paste("med_q",RecT, sep=".")
  assign(nam1,avg_q)
  assign(nam2,med_q)
} 

# end of the loop

# Build thedataframe of the results for the corresponding year and industry codes
year <- matrix(2013,54,1)
output_inop_us13 <- data.frame(cbind(wiot_US13$IndustryCode[1:54],wiot_US13$IndustryDescription[1:54],year,med_q.180,med_q.60))
names(output_inop_us13) <- c("ISIC_SNA","IndustryName","year","inop1_180","inop2_180","inop3_180","inop4_180","inop1_60","inop2_60","inop3_60","inop4_60")

# Save the dataframe on a dta file 
write.dta(output_inop_us13,"inopUS2013.dta")


######################################################################
# STEP 2: CALCULATE THE ASSOCIATED ECONOMIC LOSSES
######################################################################

# MEdian inoperability level by industry for T = 180 days
MED_180 <- as.matrix(cbind(industr_name,med_q.180))

# For the 10 first sectors that are inoperable (10 most inoperable sectors) by initial inperability level
MED_180_class <- data.frame(0,10,1)

for(i in 2:5){
  tempp <- cbind(MED_180[,1], MED_180[,i])
  temp <- tempp[order(as.numeric(tempp[,2]), decreasing = T),]
  temp <- temp[1:10,]
  MED_180_class <- cbind(MED_180_class,temp[,1],temp[,2])
}

# Keep the matrix of the 10 most inoperable sectors per initial inoperability level
MED_180_class1 <- MED_180_class[,c(4,6,8,10)]

#########################################
n_days <- 180

# Determining economic losses
daily_output <- wiot_US13$Output/365
daily_output <- S%*%daily_output

  # For initial inoperability = 0.1
  econ_loss2 <- t(t(Q_0.1_180)*as.numeric(c(daily_output)))
  econ_loss2 <- t(apply(econ_loss2, 1, cumsum))
  econ_losss2 <- cbind(econ_loss2[,180],industr_name,as.character(industr_code))
  econ_loss_class2 <-  econ_losss2[order(as.numeric(as.character(econ_losss2[,1])), decreasing = T),]
  
  # For initial inoperability = 0.2
  econ_loss3 <- t(t(Q_0.2_180)*as.numeric(c(daily_output)))
  econ_loss3 <- t(apply(econ_loss3, 1, cumsum))
  econ_losss3 <- cbind(econ_loss3[,180],industr_name,as.character(industr_code))
  econ_loss_class3 <-  econ_losss3[order(as.numeric(as.character(econ_losss3[,1])), decreasing = T),]
  
  # For initial inoperability = 0.3
  econ_loss4 <- t(t(Q_0.3_180)*as.numeric(c(daily_output)))
  econ_loss4 <- t(apply(econ_loss4, 1, cumsum))
  econ_losss4 <- cbind(econ_loss4[,180],industr_name,as.character(industr_code))
  econ_loss_class4 <-  econ_losss4[order(as.numeric(as.character(econ_losss4[,1])), decreasing = T),]
  
  # For initial inoperability = 0.4
  econ_loss1 <- t(Q_0.4_180)*as.numeric(c(daily_output))
  econ_loss1 <- t(econ_loss1)
  econ_loss1 <- t(apply(econ_loss1, 1, cumsum))
  econ_loss <- econ_loss1[,180]
  econ_loss <- cbind(econ_loss,as.character(industr_name),as.character(industr_code))
  econ_loss_class <-  econ_loss[order(as.numeric(as.character(econ_loss[,1])), decreasing = T),]

# Ranking
rank <- cbind(as.character(AVG_180_class1[,4]),round(as.numeric(as.character(MED_180_class[,11])),3),econ_loss_class[1:10,2],round(as.numeric(econ_loss_class[1:10,1]),0))
xtable(cbind(rank))

rank_inop <- cbind(as.character(AVG_180_class1[,1]),round(as.numeric(as.character(MED_180_class[,5])),3),
                   as.character(AVG_180_class1[,2]),round(as.numeric(as.character(MED_180_class[,7])),3),
                   as.character(AVG_180_class1[,3]),round(as.numeric(as.character(MED_180_class[,9])),3))
xtable(cbind(rank_inop))

rank_loss <-  cbind(econ_loss_class2[1:10,2],round(as.numeric(econ_loss_class2[1:10,1]),0),
                    econ_loss_class3[1:10,2],round(as.numeric(econ_loss_class3[1:10,1]),0),
                    econ_loss_class4[1:10,2],round(as.numeric(econ_loss_class4[1:10,1]),0))
xtable(cbind(rank_loss))

# Outputting the results

year1 <- matrix(2013,44,1)
inop <- as.matrix(c(matrix(0.1,11,1),matrix(0.2,11,1),matrix(0.3,11,1),matrix(0.4,11,1)))
econ_loss_class_other <- c(sum(as.numeric(as.character(econ_loss_class[11:54,1]))),"Other sectors","Other")
econ_loss_class_other2 <- c(sum(as.numeric(as.character(econ_loss_class2[11:54,1]))),"Other sectors","Other")
econ_loss_class_other3 <- c(sum(as.numeric(as.character(econ_loss_class3[11:54,1]))),"Other sectors","Other")
econ_loss_class_other4 <- c(sum(as.numeric(as.character(econ_loss_class4[11:54,1]))),"Other sectors","Other")

econ_loss_class_class <- rbind(econ_loss_class2[1:10,],econ_loss_class_other2,
                               econ_loss_class3[1:10,],econ_loss_class_other3,
                               econ_loss_class4[1:10,],econ_loss_class_other4,
                               econ_loss_class[1:10,],econ_loss_class_other)

econ_loss_class_class <- data.frame(cbind(econ_loss_class_class,year1,inop))
names(econ_loss_class_class) <- c("Loss","IndustryName","IndustryCode","year","inop")
econ_loss_class_class$EconLoss <- as.numeric(as.character(econ_loss_class_class$Loss))
econ_loss_class_class$Inopty <- as.numeric(as.character(econ_loss_class_class$inop))
write.dta(econ_loss_class_class,"LossUS2013Ranked.dta")

# Plotting

data_180 <- as.data.frame(rbind(Q_0.4_180, seq(1,180,1)))

  # Inoperability
  par(mfrow = c(1, 1))
  plot(seq(1,180,1),Q_0.4_180[g,],type="l",col="black", xlab="Days", ylab="Inoperability",main="Inoperability for top-ten sectors" )
  lines(seq(1,180,1),Q_0.4_180[49,], col="red")
  lines(seq(1,180,1),Q_0.4_180[42,],col="brown")
  lines(seq(1,180,1),Q_0.4_180[44,],col="orange")
  lines(seq(1,180,1),Q_0.4_180[43,],col="blue")
  lines(seq(1,180,1),Q_0.4_180[39,],col="purple")
  lines(seq(1,180,1),Q_0.4_180[41,],col="green")
  lines(seq(1,180,1),Q_0.4_180[29,],col="yellow")
  lines(seq(1,180,1),Q_0.4_180[45,],col="grey")
  lines(seq(1,180,1),Q_0.4_180[50,],col="cyan")
  legend(140,0.4, c("FINANCE","ADMIN","AUXIL","LEGAL","REAL","ICT","INSURANCE","WHOLESALE","ARCHITECT","PUBLIC"),
         lwd=c(2.5,2.5),col=c("black","red", "brown", "orange","blue",
                              "purple","green","yellow","grey","cyan"),
         cex = 0.6)
  
  # Economic losses
  plot(seq(1,180,1),econ_loss1[g,],type="l",col="black", xlab="Days", ylab="Economic loss",main="Economic losses for top-ten sectors" )
  lines(seq(1,180,1),econ_loss1[42,], col="red")
  lines(seq(1,180,1),econ_loss1[49,],col="brown")
  lines(seq(1,180,1),econ_loss1[44,],col="orange")
  lines(seq(1,180,1),econ_loss1[39,],col="blue")
  lines(seq(1,180,1),econ_loss1[41,],col="purple")
  lines(seq(1,180,1),econ_loss1[43,],col="green")
  lines(seq(1,180,1),econ_loss1[29,],col="yellow")
  lines(seq(1,180,1),econ_loss1[45,],col="grey")
  lines(seq(1,180,1),econ_loss1[53,],col="cyan")
  legend(0,29000, c("FINANCE","AUXIL","ADMIN","LEGAL","ICT","INSURANCE","REAL","WHOLESALE","ARCHITECT","OTHER"),
         lwd=c(2.5,2.5),col=c("black","red", "brown", "orange","blue",
                              "purple","green","yellow","grey","cyan"),
         cex = 0.6)

######################################################################################################################################################
