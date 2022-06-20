library(gridExtra)

#Plot True Gammas
sim.set <- 1; sim.rep<-1
source("GenSim.R")
gamma.true.sim1 <- gamma.true
gammaTrueplot1 <- ggplot(mapping = aes(x = t, y=gamma.true.sim1)) + geom_line() + theme_bw() + ggtitle("Simulation 1") + theme(text = element_text(size=27)) +ylab("gamma") +xlab("Location") + theme(axis.text.x = element_text(size = 24, hjust = 1))
sim.set <- 2; sim.rep<-1
source("GenSim.R")
gamma.true.sim2 <- gamma.true
gammaTrueplot2 <- ggplot(mapping = aes(x = t, y=gamma.true.sim2)) + geom_line() + theme_bw() + ggtitle("Simulation 2") + theme(text = element_text(size=27)) +ylab("gamma") +xlab("Location") + theme(axis.text.x = element_text(size = 24, hjust = 1))
sim.set <- 3; sim.rep<-1
source("GenSim.R")
gamma.true.sim3 <- gamma.true
gammaTrueplot3 <- ggplot(mapping = aes(x = t, y=gamma.true.sim3)) + geom_line() + theme_bw() + ggtitle("Simulation 3") + theme(text = element_text(size=27)) +ylab("gamma") +xlab("Location") + theme(axis.text.x = element_text(size = 24, hjust = 1))
sim.set <- 4; sim.rep<-1
M <- matrix(0.95,nrow=5000,ncol=5000)
diag(M) <- rep(1,5000)
cholM <- t(chol(M))
source("GenSim.R")
gamma.true.sim4 <- gamma.true
gammaTrueplot4 <- ggplot(mapping = aes(x = t, y=gamma.true.sim4)) + geom_line() + theme_bw() + ggtitle("Simulation 4") + theme(text = element_text(size=27)) +ylab("gamma") +xlab("Location") + theme(axis.text.x = element_text(size = 24, hjust = 1))
pdf("NonStatTrueGammas.pdf",width=18, height=16)
grid.arrange(gammaTrueplot1,gammaTrueplot2,gammaTrueplot3,gammaTrueplot4,nrow=2,ncol=2)
dev.off()


#Plot True Mean Functions
sim.set <- 1; sim.rep <- 1
source("GenSim.R")
mu.true1.sim1 <- mu.true1
mu.true0.sim1 <- mu.true0
DF1 <-cbind(seq(1,5000,1),c(mu.true1,mu.true0),rep(c(1,0),rep(5000,2)))
DF1 <- rbind(DF1, cbind(1:50,mu.true1[1:50],rep(2,50)))
DF1 <- rbind(DF1, cbind(301:700,mu.true1[301:700],rep(3,length(301:700))))
DF1 <- rbind(DF1, cbind(951:1250,mu.true1[951:1250],rep(4,length(951:1250))))
DF1 <- rbind(DF1, cbind(1501:1600,mu.true1[1501:1600],rep(5,length(1501:1600))))
DF1 <- rbind(DF1, cbind(1851:2000,mu.true1[1851:2000],rep(6,length(1851:2000))))
DF1 <- rbind(DF1, cbind(2251:3010,mu.true1[2251:3010],rep(7,length(2251:3010))))
DF1 <- rbind(DF1, cbind(3261:3800,mu.true1[3261:3800],rep(8,length(3261:3800))))
DF1 <- rbind(DF1, cbind(4051:4400,mu.true1[4051:4400],rep(9,length(4051:4400))))
DF1 <- rbind(DF1, cbind(4651:5000,mu.true1[4651:5000],rep(10,length(4651:5000))))
DF1 <- as.data.frame(DF1)
names(DF1) <- c("Location","mu","Group")
DF1$Group <- as.factor(DF1$Group)
muplots1 <-ggplot(DF1, aes(Location, mu, colour=as.factor(Group) )) + geom_line(aes(group=as.factor(Group) )) + scale_alpha_manual(values=rep(0.01,11)) + scale_colour_manual(values=c("green","red",rep("blue",9)))+theme_bw() + theme(legend.position="none") + ggtitle("Simulation 1") + theme(strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1)) 


sim.set <- 2; sim.rep <- 1
source("GenSim.R")
mu.true1.sim2 <- mu.true1
mu.true0.sim2 <- mu.true0
DF2 <-cbind(seq(1,5000,1),c(mu.true1,mu.true0),rep(c(1,0),rep(5000,2)))
DF2 <- rbind(DF2, cbind(1:1600,mu.true1[1:1600],rep(2,length(1:1600))))
DF2 <- rbind(DF2, cbind(1726:3200,mu.true1[1726:3200],rep(3,length(1726:3200))))
DF2 <- rbind(DF2, cbind(3376:5000,mu.true1[3376:5000],rep(4,length(3376:5000))))
DF2 <- as.data.frame(DF2)
names(DF2) <- c("Location","mu","Group")
DF2$Group <- as.factor(DF2$Group)
muplots2 <-ggplot(DF2, aes(Location, mu, colour=as.factor(Group) )) + geom_line(aes(group=as.factor(Group))) + scale_colour_manual(values=c("green","red",rep("blue",3)))+theme_bw() + theme(legend.position="none") + ggtitle("Simulation 2") + theme(strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))


sim.set <- 3; sim.rep <- 1
source("GenSim.R")
mu.true1.sim3 <- mu.true1
mu.true0.sim3 <- mu.true0
DF3 <-cbind(seq(1,5000,1),c(mu.true1,mu.true0),rep(c(1,0),rep(5000,2)))
DF3 <- rbind(DF3, cbind(1:1200,mu.true1[1:1200],rep(2,length(1:1200))))
DF3 <- rbind(DF3, cbind(1451:2200,mu.true1[1451:2200],rep(3,length(1451:2200))))
DF3 <- rbind(DF3, cbind(2451:5000,mu.true1[2451:5000],rep(4,length(2451:5000))))
DF3 <- as.data.frame(DF3)
names(DF3) <- c("Location","mu","Group")
DF3$Group <- as.factor(DF3$Group)
muplots3 <-ggplot(DF3, aes(Location, mu, colour=as.factor(Group) )) + geom_line(aes(group=as.factor(Group))) + scale_colour_manual(values=c("green","red",rep("blue",3)))+theme_bw() + theme(legend.position="none") + ggtitle("Simulation 3") + theme(strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))



sim.set <- 4; sim.rep <- 1
M <- matrix(0.95,nrow=5000,ncol=5000)
diag(M) <- rep(1,5000)
cholM <- t(chol(M))
source("GenSim.R")
mu.true1.sim4 <- mu.true1
mu.true0.sim4 <- mu.true0
DF4 <-cbind(seq(1,5000,1),c(mu.true1,mu.true0),rep(c(1,0),rep(5000,2)))
DF4 <- rbind(DF4, cbind(1:1400,mu.true1[1:1400],rep(2,length(1:1400))))
DF4 <- rbind(DF4, cbind(1701:3400,mu.true1[1701:3400],rep(3,length(1701:3400))))
DF4 <- rbind(DF4, cbind(3701:5000,mu.true1[3701:5000],rep(4,length(3701:5000))))
DF4 <- as.data.frame(DF4)
names(DF4) <- c("Location","mu","Group")
DF4$Group <- as.factor(DF4$Group)
muplots4 <-ggplot(DF4, aes(Location, mu, colour=as.factor(Group) )) + geom_line(aes(group=as.factor(Group))) + scale_colour_manual(values=c("green","red",rep("blue",3)))+theme_bw() + theme(legend.position="none") + ggtitle("Simulation 4") + theme(strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))


pdf("NonStatTrueMus.pdf",width=18, height=16)
grid.arrange(muplots1,muplots2,muplots3,muplots4,nrow=2,ncol=2)
dev.off()

#Plot True Noise Variance Functions
sim.set <- 1; sim.rep <- 1
source("GenSim.R")
sigma.eps.true1.sim1 <- sigma.eps.true1
sigma.eps.true0.sim1 <- sigma.eps.true0
dataframeEps1 <- data.frame(c(sigma.eps.true1.sim1,sigma.eps.true0.sim1), factor(c(rep("Positive",length(sigma.eps.true1.sim1)), rep("Negative",length(sigma.eps.true0.sim1)))), rep(1:5000,2) )
names(dataframeEps1) <- c("value","status","Timepoint")
dataframeEps1[,2] <- factor(dataframeEps1[,2],levels=c("Positive","Negative"))
Epsplots1 <- ggplot(dataframeEps1, aes(x = Timepoint, y = value, col=status)) + geom_line(aes(color = status), size = 0.5, alpha=0.4) +scale_color_manual(values = c("red", "blue")) + theme_bw() + ggtitle("Simulation 1") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))

sim.set <- 2; sim.rep <- 1
source("GenSim.R")
sigma.eps.true1.sim2 <- sigma.eps.true1
sigma.eps.true0.sim2 <- sigma.eps.true0
dataframeEps2 <- data.frame(c(sigma.eps.true1.sim2,sigma.eps.true0.sim2), factor(c(rep("Positive",length(sigma.eps.true1.sim1)), rep("Negative",length(sigma.eps.true0.sim1)))), rep(1:5000,2) )
names(dataframeEps2) <- c("value","status","Timepoint")
dataframeEps2[,2] <- factor(dataframeEps2[,2],levels=c("Positive","Negative"))
Epsplots2 <- ggplot(dataframeEps2, aes(x = Timepoint, y = value, col=status)) + geom_line(aes(color = status), size = 0.5, alpha=0.4) +scale_color_manual(values = c("red", "blue")) + theme_bw() + ggtitle("Simulation 2") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))

sim.set <- 3; sim.rep <- 1
source("GenSim.R")
sigma.eps.true1.sim3 <- sigma.eps.true1
sigma.eps.true0.sim3 <- sigma.eps.true0
dataframeEps3 <- data.frame(c(sigma.eps.true1.sim3,sigma.eps.true0.sim3), factor(c(rep("Positive",length(sigma.eps.true1.sim1)), rep("Negative",length(sigma.eps.true0.sim1)))), rep(1:5000,2) )
names(dataframeEps3) <- c("value","status","Timepoint")
dataframeEps3[,2] <- factor(dataframeEps3[,2],levels=c("Positive","Negative"))
Epsplots3 <- ggplot(dataframeEps3, aes(x = Timepoint, y = value, col=status)) + geom_line(aes(color = status), size = 0.5, alpha=0.4) +scale_color_manual(values = c("red", "blue")) + theme_bw() + ggtitle("Simulation 3") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))

sim.set <- 4; sim.rep <- 1
M <- matrix(0.7,nrow=5000,ncol=5000)
diag(M) <- rep(1,5000)
cholM <- t(chol(M))
source("GenSim.R")
sigma.eps.true1.sim4 <- rep(0,n.varbs)
sigma.eps.true0.sim4 <- rep(0,n.varbs)
dataframeEps4 <- data.frame(c(sigma.eps.true1.sim4,sigma.eps.true0.sim4), factor(c(rep("Positive",length(sigma.eps.true1.sim4)), rep("Negative",length(sigma.eps.true0.sim4)))), rep(1:5000,2) )
names(dataframeEps4) <- c("value","status","Timepoint")
dataframeEps4[,2] <- factor(dataframeEps4[,2],levels=c("Positive","Negative"))
Epsplots4 <- ggplot(dataframeEps4, aes(x = Timepoint, y = value, col=status)) + geom_line(aes(color = status), size = 0.5, alpha=0.4) +scale_color_manual(values = c("red", "blue")) + theme_bw() + ggtitle("Simulation 4") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))
pdf("NonStatTrueSigmaEps.pdf",width=18, height=16)
grid.arrange(Epsplots1,Epsplots2,Epsplots3,Epsplots4,nrow=2,ncol=2)
dev.off()

#Plot true S
sim.set <- 1; sim.rep<-1
source("GenSim.R")
dataframeS1 <- data.frame(S.true,1:4999)
names(dataframeS1) <- c("S","Location")
Splots1 <- ggplot(dataframeS1, aes(x = Location, y = S)) + geom_line(size = 0.5)  + theme_bw() + ggtitle("Simulation 1") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))


sim.set <- 2; sim.rep<-1
source("GenSim.R")
dataframeS2 <- data.frame(S.true,1:4999)
names(dataframeS2) <- c("S","Location")
Splots2 <- ggplot(dataframeS2, aes(x = Location, y = S)) + geom_line(size = 0.5)  + theme_bw() + ggtitle("Simulation 2") + theme(legend.position = "bottom", strip.text = element_text(size = 24)) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, hjust = 1))
pdf("NonStatTrueS.pdf",width=14, height=9)
grid.arrange(Splots1,Splots2,nrow=1,ncol=2)
dev.off()


#Plot Classification Errors
load("NonStatSimulation1.RData")
Accmatrix1 <- 100 - cbind(Errors1,Errors2,Errors3, Errors4,Errors5, Errors6, Errors7 , Errors8)/500*100
colnames(Accmatrix1) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1")
datAcc1 = reshape2::melt(Accmatrix1,varnames=c("Iteration","Method"),value.name="ClassificationAccuracy")
names(datAcc1) <- c("Iteration","Classifier","ClassificationAccuracy")
Acc_plot1 <- ggplot(datAcc1, aes(x = Classifier, y = ClassificationAccuracy,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("Classification accuracy (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",7))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Simulation 1")
load("NonStatSimulation2.RData")
Accmatrix2 <- 100 - cbind(Errors1,Errors2,Errors3, Errors4,Errors5, Errors6, Errors7 , Errors8)/500*100
colnames(Accmatrix2) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1")
datAcc2 = reshape2::melt(Accmatrix2,varnames=c("Iteration","Method"),value.name="ClassificationAccuracy")
names(datAcc2) <- c("Iteration","Classifier","ClassificationAccuracy")
Acc_plot2 <- ggplot(datAcc2, aes(x = Classifier, y = ClassificationAccuracy,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("Classification accuracy (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",7))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Simulation 2")
load("NonStatSimulation3.RData")
Accmatrix3 <- 100 - cbind(Errors1,Errors2,Errors3, Errors4,Errors5, Errors6, Errors7 , Errors8)/500*100
colnames(Accmatrix3) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1")
datAcc3 = reshape2::melt(Accmatrix3,varnames=c("Iteration","Method"),value.name="ClassificationAccuracy")
names(datAcc3) <- c("Iteration","Classifier","ClassificationAccuracy")
Acc_plot3 <- ggplot(datAcc3, aes(x = Classifier, y = ClassificationAccuracy,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("Classification accuracy (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",7))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Simulation 3")
load("NonStatSimulation4.RData")
Accmatrix4 <- 100 - cbind(Errors1,Errors2,Errors3, Errors4,Errors5, Errors6, Errors7 , Errors8)/500*100
colnames(Accmatrix4) <- c("GPDA","VNPDA","penLDA-FL","RandomForest","VLDA","SparseLDA","SVM-L2","SVM-L1")
datAcc4 = reshape2::melt(Accmatrix4,varnames=c("Iteration","Method"),value.name="ClassificationAccuracy")
names(datAcc4) <- c("Iteration","Classifier","ClassificationAccuracy")
Acc_plot4 <- ggplot(datAcc4, aes(x = Classifier, y = ClassificationAccuracy,fill = Classifier)) + geom_boxplot()+ xlab(" ") + ylab("Classification accuracy (%)") + theme_bw() + scale_fill_manual(values=c("yellow",rep("gray",7))) + theme(text = element_text(size=27)) + theme(axis.text.x = element_text(size = 24, angle = 90, hjust = 1)) + theme(axis.text.y = element_text(size = 24))  +theme(legend.position = "none", strip.text = element_text(size = 24))+theme(panel.border = element_rect(color = "black", fill = NA, size = 1.5), strip.background = element_rect(color = "black", size = 1.5))+ggtitle("Simulation 4")
pdf("NonStatSimulationClassification.pdf",width=18, height=16)
grid.arrange(Acc_plot1,Acc_plot2,Acc_plot3,Acc_plot4,nrow=2,ncol=2)
dev.off()