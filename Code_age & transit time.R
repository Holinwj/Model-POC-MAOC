library(SoilR)

################ using mean kp and km from the paper (Zhou et al., 2024, NC); kp=1/23, km=1/129 (0-100 cm) ############
############################################### 1. without input to MAOC ########################################### 
rm(list = ls())

############################### 2pseries1
k2s1=c(1/23,1/129)
A2s1=diag(-k2s1)
A2s1[2,1]=k2s1[1]*0.017685
u2s1=matrix(c(100,0),ncol=1)

############################### 3pseries1
k3s1=c(0.023578,1/23,1/129)
A3s1=diag(-k3s1)
A3s1[2,1]=k3s1[1]*0.061840
A3s1[3,2]=k3s1[2]*0.017685
u3s1=matrix(c(100,0,0),ncol=1)

############################### 3pseries2
k3s2=c(0.070734,1/23,1/129)
A3s2=diag(-k3s2)
A3s2[2,1]=k3s2[1]*0.06184
A3s2[3,2]=k3s2[2]*0.017685
u3s2=matrix(c(100,0,0),ncol=1)

############################### 3pseries3
k3s3=c(0.070734,1/23,1/129)
A3s3=diag(-k3s3)
A3s3[2,1]=k3s3[1]*0.185520
A3s3[3,2]=k3s3[2]*0.017685
u3s3=matrix(c(100,0,0),ncol=1)

############################### 3pseries4
k3s4=c(0.070734,1/23,1/129)
A3s4=diag(-k3s4)
A3s4[2,1]=k3s4[1]*0.185520
A3s4[3,2]=k3s4[2]*0.053055
u3s4=matrix(c(100,0,0),ncol=1)

##################################### Age and transit time
SA2s1=systemAge(A=A2s1, u=u2s1)
SA3s1=systemAge(A=A3s1, u=u3s1)
SA3s2=systemAge(A=A3s2, u=u3s2)
SA3s3=systemAge(A=A3s3, u=u3s3)
SA3s4=systemAge(A=A3s4, u=u3s4)

TT2s1=transitTime(A=A2s1, u=u2s1)
TT3s1=transitTime(A=A3s1, u=u3s1)
TT3s2=transitTime(A=A3s2, u=u3s2)
TT3s3=transitTime(A=A3s3, u=u3s3)
TT3s4=transitTime(A=A3s4, u=u3s4)

SA2s1$meanSystemAge
SA3s1$meanSystemAge
SA3s2$meanSystemAge
SA3s3$meanSystemAge
SA3s4$meanSystemAge

TT2s1$meanTransitTime
TT3s1$meanTransitTime
TT3s2$meanTransitTime
TT3s3$meanTransitTime
TT3s4$meanTransitTime

########################################## plot age
plot(xlim=c(0, 100),ylim=c(-0.0005, 0.07),SA2s1$systemAgeDensity, type="l", col='black', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA2s1$meanSystemAge, lwd=2,lty=2, col='black')
lines(SA3s1$systemAgeDensity, type="l", col='#f1c6a0', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s1$meanSystemAge, lwd=2,lty=2, col='#f1c6a0')
lines(SA3s2$systemAgeDensity, type="l", col='#e89c8f', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s2$meanSystemAge, lwd=2,lty=2, col='#e89c8f')
lines(SA3s3$systemAgeDensity, type="l", col='#d67f7c', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s3$meanSystemAge, lwd=2,lty=2, col='#d67f7c')
lines(SA3s4$systemAgeDensity, type="l", col='#c25e6c', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s4$meanSystemAge, lwd=2,lty=2, col='#c25e6c')

########################################## plot transit time
plot(xlim=c(0, 100),ylim=c(-0.0005, 0.07),TT2s1$transitTimeDensity, type="l", col='black', xlab="Transit time", ylab="Density function",lwd=2)
abline(v=TT2s1$meanTransitTime, lwd=2,lty=2, col='black')
lines(TT3s1$transitTimeDensity, type="l", col='#a1d3e0', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s1$meanTransitTime, lwd=2,lty=2, col='#a1d3e0')
lines(TT3s2$transitTimeDensity, type="l", col='#74b1c1', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s2$meanTransitTime, lwd=2,lty=2, col='#74b1c1')
lines(TT3s3$transitTimeDensity, type="l", col='#4a8ea3', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s3$meanTransitTime, lwd=2,lty=2, col='#4a8ea3')
lines(TT3s4$transitTimeDensity, type="l", col='#2f6d8c', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s4$meanTransitTime, lwd=2,lty=2, col='#2f6d8c')

############################################### 2. with input to MAOC ########################################### 
rm(list = ls())

############################### 2pseries1
k2s1=c(1/23,1/129)
A2s1=diag(-k2s1)
A2s1[2,1]=k2s1[1]*0.017685
u2s1=matrix(c(90,10),ncol=1)

############################### 2pseries2
k2s2=c(1/23,1/129)
A2s2=diag(-k2s2)
A2s2[2,1]=k2s2[1]*0.017685
u2s2=matrix(c(70,30),ncol=1)

############################### 3pseries1
k3s1=c(0.023578,1/23,1/129)
A3s1=diag(-k3s1)
A3s1[2,1]=k3s1[1]*0.061840
A3s1[3,2]=k3s1[2]*0.017685
A3s1[3,1]=k3s1[1]*0.086723
u3s1=matrix(c(100,0,0),ncol=1)

############################### 3pseries2
k3s2=c(0.070734,1/23,1/129)
A3s2=diag(-k3s2)
A3s2[2,1]=k3s2[1]*0.061840
A3s2[3,2]=k3s2[2]*0.017685
A3s2[3,1]=k3s2[1]*0.086723
u3s2=matrix(c(100,0,0),ncol=1)

############################### 3pseries3
k3s3=c(0.070734,1/23,1/129)
A3s3=diag(-k3s3)
A3s3[2,1]=k3s3[1]*0.185520
A3s3[3,2]=k3s3[2]*0.017685
A3s3[3,1]=k3s3[1]*0.086723
u3s3=matrix(c(100,0,0),ncol=1)

############################### 3pseries4
k3s4=c(0.070734,1/23,1/129)
A3s4=diag(-k3s4)
A3s4[2,1]=k3s4[1]*0.185520
A3s4[3,2]=k3s4[2]*0.053055
A3s4[3,1]=k3s4[1]*0.086723
u3s4=matrix(c(100,0,0),ncol=1)

############################### 3pseries5
k3s5=c(0.070734,1/23,1/129)
A3s5=diag(-k3s5)
A3s5[2,1]=k3s5[1]*0.185520
A3s5[3,2]=k3s5[2]*0.053055
A3s5[3,1]=k3s5[1]*0.260169
u3s5=matrix(c(100,0,0),ncol=1)

##################################### Age and transit time
SA2s1=systemAge(a = seq(0, 230),A=A2s1, u=u2s1)
SA2s2=systemAge(a = seq(0, 230),A=A2s2, u=u2s2)
SA3s1=systemAge(a = seq(0, 230),A=A3s1, u=u3s1)
SA3s2=systemAge(a = seq(0, 230),A=A3s2, u=u3s2)
SA3s3=systemAge(a = seq(0, 230),A=A3s3, u=u3s3)
SA3s4=systemAge(a = seq(0, 230),A=A3s4, u=u3s4)
SA3s5=systemAge(a = seq(0, 230),A=A3s5, u=u3s5)

TT2s1=transitTime(A=A2s1, u=u2s1)
TT2s2=transitTime(A=A2s2, u=u2s2)
TT3s1=transitTime(A=A3s1, u=u3s1)
TT3s2=transitTime(A=A3s2, u=u3s2)
TT3s3=transitTime(A=A3s3, u=u3s3)
TT3s4=transitTime(A=A3s4, u=u3s4)
TT3s5=transitTime(A=A3s5, u=u3s5)

SA2s1$meanSystemAge
SA2s2$meanSystemAge
SA3s1$meanSystemAge
SA3s2$meanSystemAge
SA3s3$meanSystemAge
SA3s4$meanSystemAge
SA3s5$meanSystemAge

TT2s1$meanTransitTime
TT2s2$meanTransitTime
TT3s1$meanTransitTime
TT3s2$meanTransitTime
TT3s3$meanTransitTime
TT3s4$meanTransitTime
TT3s5$meanTransitTime

########################################## plot age
plot(xlim=c(0, 230),ylim=c(-0.0005, 0.04),SA2s1$systemAgeDensity, type="l", col='#c0c0c0', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA2s1$meanSystemAge, lwd=2,lty=2, col='#c0c0c0')
lines(SA2s2$systemAgeDensity, type="l", col='#787878', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA2s2$meanSystemAge, lwd=2,lty=2, col='#787878')
lines(SA3s1$systemAgeDensity, type="l", col='#f1c6a0', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s1$meanSystemAge, lwd=2,lty=2, col='#f1c6a0')
lines(SA3s2$systemAgeDensity, type="l", col='#e89c8f', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s2$meanSystemAge, lwd=2,lty=2, col='#e89c8f')
lines(SA3s3$systemAgeDensity, type="l", col='#d67f7c', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s3$meanSystemAge, lwd=2,lty=2, col='#d67f7c')
lines(SA3s4$systemAgeDensity, type="l", col='#c25e6c', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s4$meanSystemAge, lwd=2,lty=2, col='#c25e6c')
lines(SA3s5$systemAgeDensity, type="l", col='#a94d5c', xlab="Age", ylab="Density function",lwd=2)
abline(v=SA3s5$meanSystemAge, lwd=2,lty=2, col='#a94d5c')

########################################## plot transit time
plot(xlim=c(0, 100),ylim=c(-0.0005, 0.04),TT2s1$transitTimeDensity, type="l", col='#c0c0c0', xlab="Transit time", ylab="Density function",lwd=2)
abline(v=TT2s1$meanTransitTime, lwd=2,lty=2, col='#c0c0c0')
lines(TT2s2$transitTimeDensity, type="l", col='#787878', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT2s2$meanTransitTime, lwd=2,lty=2, col='#787878')
lines(TT3s1$transitTimeDensity, type="l", col='#a1d3e0', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s1$meanTransitTime, lwd=2,lty=2, col='#a1d3e0')
lines(TT3s2$transitTimeDensity, type="l", col='#74b1c1', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s2$meanTransitTime, lwd=2,lty=2, col='#74b1c1')
lines(TT3s3$transitTimeDensity, type="l", col='#4a8ea3', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s3$meanTransitTime, lwd=2,lty=2, col='#4a8ea3')
lines(TT3s4$transitTimeDensity, type="l", col='#2f6d8c', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s4$meanTransitTime, lwd=2,lty=2, col='#2f6d8c')
lines(TT3s5$transitTimeDensity, type="l", col='#1c4b65', xlab="Age", ylab="Density function",lwd=2)
abline(v=TT3s5$meanTransitTime, lwd=2,lty=2, col='#1c4b65')


















