library(SoilR)
library(FME)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)

# 2-pool models 

data = read_excel("./2_pool_model/data.xlsx")

twopsFit=function(timeSeries, initial, inipars){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  tt=seq(from=0, to=unlist(tail(complete[,1],1)), length.out = 500)
  
  Func2ps=function(pars){
    mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], 
                               C0=c(initial[1,2]*initial[3,2], initial[2,2]*initial[3,2]), In=0)
    Ct=SoilR::getC(mod)
    rt = SoilR::getAccumulatedRelease(mod)
    return(data.frame(time=tt, rt=rowSums(rt), Ct=rowSums(Ct), pom = Ct[, 1], maom = Ct[, 2]))
  }
  
  costFunc=function(pars){
    output=Func2ps(pars)
    cost1 = modCost(model = output, obs = as.data.frame(rt), x="time", err = "sd" )
    cost2 = modCost(model = output, obs = as.data.frame(pom), x="time", err = "sd", cost = cost1)
    return(modCost(model=output, obs=as.data.frame(maom), x="time", err = "sd", cost = cost2)) 
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Nelder-Mead", lower=c(0,0,0), upper=c(Inf, Inf, 1))
  bestMod=Func2ps(pars=Fit$par) 

  AIC=((n+2*(length(Fit$par)+1))/(n))+log(Fit$ms) 
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3],  
                                    C0=c(initial[1,2]*initial[3,2], initial[2,2]*initial[3,2]), In=0)
  A=matrix(c(-Fit$par[1],Fit$par[3]*Fit$par[1],0,-Fit$par[2]),ncol=2)
  u=matrix(c(initial[1,2]*initial[3,2], initial[2,2]*initial[3,2]),ncol=1)
  TT=transitTime(A,u)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TT=TT, AIC=AIC))
}

#0-5 cm depth control samples Fitted to C contents

rt=na.omit(data[-c(1:3),c(1,2,6)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./2_pool_model/data_C.xlsx")
pom = data_C[,c(1:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = as.data.frame(data[c(1:3),c(1,2)])

# fit the two pool series model with known initial C values 
M1=twopsFit(rt, initial, inipars=c(0.005, 0.002, 0.01))  # fit C pool
t1=M1$SoilRmodel@times                   
C1=getC(M1$SoilRmodel)                
M1$FMEmodel$par                     
M1$FMEmodel$ms

t1 = as.data.frame(t1)
t1$time = seq(from = 1, to = 182, length.out = 500)
C1 = as.data.frame(C1)
C1$time = seq(from = 1, to = 182, length.out=500)
C1$total <- C1$V1 + C1$V2

g1 = ggplot(C1, aes(x = time)) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_point(data = data.frame(time = maom$time, C = maom$maom), 
             aes(x = time, y = C), size=2, color = "#0e2f4e", shape = 19) +
  geom_point(data = data.frame(time = pom$time, C = pom$pom), 
             aes(x = time, y = C), size=2, color = "#8a3c4c", shape = 19) +
  geom_errorbar(data = pom, aes(x = time, ymin = pom - sd, ymax = pom + sd), 
                width = 0.1, color = "#8a3c4c") +
  geom_errorbar(data = maom, aes(x = time, ymin = maom - sd, ymax = maom + sd), 
                width = 0.1, color = "#0e2f4e") +
  scale_color_manual(values = c("C total" = "black", "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "POC", "MAOC")) +
  labs(y = expression("C content (mg C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") +
  scale_y_continuous(limits = c(0, 95))+
  theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g1

rt1=getAccumulatedRelease(M1$SoilRmodel)
rt1=as.data.frame(rt1)
rt1$time = seq(from = 1, to = 182, length.out=500)
rt1$total <- rt1$V1 + rt1$V2

g2= ggplot(rt1, aes(x = time)) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_point(data = rt, aes(x = time, y = rt), size = 2, color = "black", shape = 19) +
  geom_errorbar(data = rt, aes(x = time, ymin = rt - sd, ymax = rt + sd), 
                width = 0.1, color = "black") +
  scale_color_manual(values = c("C total" = "black", 
                                "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/C_C05_2ps.rds")
saveRDS(g2, "./outputs/Resp_C05_2ps.rds")
saveRDS(M1, "./2_pool_model/C05_2ps_2.rds")

#################################################################################################

# 0-5 cm depth control samples Fitted to C respiration

rt=na.omit(data[-c(1:3),c(1,2,6)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./2_pool_model/data_C.xlsx")
pom = data_C[,c(1:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = as.data.frame(data[c(1:3),c(1,2)])

# fit the two pool series model with known initial C values 
M1=twopsFit(rt, initial, inipars=c(0.01, 0.2, 0.1)) # fit resp
t1=M1$SoilRmodel@times                    
C1=getC(M1$SoilRmodel)                
M1$FMEmodel$par                     
M1$FMEmodel$ms

t1 = as.data.frame(t1)
t1$time = seq(from = 1, to = 182, length.out = 500)
C1 = as.data.frame(C1)
C1$time = seq(from = 1, to = 182, length.out=500)
C1$total <- C1$V1 + C1$V2

g1 = ggplot(C1, aes(x = time)) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_point(data = data.frame(time = maom$time, C = maom$maom), 
             aes(x = time, y = C), size=2, color = "#0e2f4e", shape = 19) +
  geom_point(data = data.frame(time = pom$time, C = pom$pom), 
             aes(x = time, y = C), size=2, color = "#8a3c4c", shape = 19) +
  geom_errorbar(data = pom, aes(x = time, ymin = pom - sd, ymax = pom + sd), 
                width = 0.1, color = "#8a3c4c") +
  geom_errorbar(data = maom, aes(x = time, ymin = maom - sd, ymax = maom + sd), 
                width = 0.1, color = "#0e2f4e") +
  scale_color_manual(values = c("C total" = "black", "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "POC", "MAOC")) +
  labs(y = expression("C content (mg C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") +
  scale_y_continuous(limits = c(0, 95))+
  theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g1

rt1=getAccumulatedRelease(M1$SoilRmodel)
rt1=as.data.frame(rt1)
rt1$time = seq(from = 1, to = 182, length.out=500)
rt1$total <- rt1$V1 + rt1$V2

g2= ggplot(rt1, aes(x = time)) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_point(data = rt, aes(x = time, y = rt), size = 2, color = "black", shape = 19) +
  geom_errorbar(data = rt, aes(x = time, ymin = rt - sd, ymax = rt + sd), 
                width = 0.1, color = "black") +
  scale_color_manual(values = c("C total" = "black", 
                                "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/C_C05_2ps_resp.rds")
saveRDS(g2, "./outputs/Resp_C05_2ps_resp.rds")
saveRDS(M1, "./2_pool_model/C05_2ps_resp.rds")

#################################################################################################

#0-5 cm depth litter treatment samples Fit C pool

rt=na.omit(data[-c(1:3),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./2_pool_model/data_C.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = as.data.frame(data[c(1:3),c(1,4)])

# fit the two pool series model with known initial C values 
M1=twopsFit(rt, initial, inipars=c(0.001, 0.0001, 0.0001)) # C pool
t1=M1$SoilRmodel@times    
C1=getC(M1$SoilRmodel)             
M1$FMEmodel$par
M1$FMEmodel$ms

t1 = as.data.frame(t1)
t1$time = seq(from = 1, to = 182, length.out = 500)
C1 = as.data.frame(C1)
C1$time = seq(from = 1, to = 182, length.out=500)
C1$total <- C1$V1 + C1$V2


g1 = ggplot(C1, aes(x = time)) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_point(data = data.frame(time = maom$time, C = maom$maom), 
             aes(x = time, y = C), size=2, color = "#0e2f4e", shape = 19) +
  geom_point(data = data.frame(time = pom$time, C = pom$pom), 
             aes(x = time, y = C), size=2, color = "#8a3c4c", shape = 19) +
  geom_errorbar(data = pom, aes(x = time, ymin = pom - sd, ymax = pom + sd), 
                width = 0.1, color = "#8a3c4c") +
  geom_errorbar(data = maom, aes(x = time, ymin = maom - sd, ymax = maom + sd), 
                width = 0.1, color = "#0e2f4e") +
  scale_color_manual(values = c("C total" = "black", "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "POC", "MAOC")) +
  labs(y = expression("C content (mg C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") +
  theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g1

rt1=getAccumulatedRelease(M1$SoilRmodel)
rt1=as.data.frame(rt1)
rt1$time = seq(from = 1, to = 182, length.out=500)
rt1$total <- rt1$V1 + rt1$V2

g2= ggplot(rt1, aes(x = time)) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_point(data = rt, aes(x = time, y = rt), size = 2, color = "black", shape = 19) +
  geom_errorbar(data = rt, aes(x = time, ymin = rt - sd, ymax = rt + sd), 
                width = 0.1, color = "black") +
  scale_color_manual(values = c("C total" = "black", 
                                "POM" = "#8a3c4c", "MAOM" = "#0e2f4e")) +
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/C_L05_2ps.rds")
saveRDS(g2, "./outputs/Resp_L05_2ps.rds")
saveRDS(M1, "./2_pool_model/L05_2ps_2.rds")

#################################################################################################

#0-5 cm depth litter treatment samples Fit C respiration

rt=na.omit(data[-c(1:3),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./2_pool_model/data_C.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = as.data.frame(data[c(1:3),c(1,4)])

# fit the two pool series model with known initial C values 
M1=twopsFit(rt, initial, inipars=c(0.03, 0.025, 0.07)) #Fit resp 
t1=M1$SoilRmodel@times    
C1=getC(M1$SoilRmodel)              
M1$FMEmodel$par
M1$FMEmodel$ms

t1 = as.data.frame(t1)
t1$time = seq(from = 1, to = 182, length.out = 500)
C1 = as.data.frame(C1)
C1$time = seq(from = 1, to = 182, length.out=500)
C1$total <- C1$V1 + C1$V2


g1 = ggplot(C1, aes(x = time)) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_point(data = data.frame(time = maom$time, C = maom$maom), 
             aes(x = time, y = C), size=2, color = "#0e2f4e", shape = 19) +
  geom_point(data = data.frame(time = pom$time, C = pom$pom), 
             aes(x = time, y = C), size=2, color = "#8a3c4c", shape = 19) +
  geom_errorbar(data = pom, aes(x = time, ymin = pom - sd, ymax = pom + sd), 
                width = 0.1, color = "#8a3c4c") +
  geom_errorbar(data = maom, aes(x = time, ymin = maom - sd, ymax = maom + sd), 
                width = 0.1, color = "#0e2f4e") +
  scale_color_manual(values = c("C total" = "black", "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "POC", "MAOC")) +
  labs(y = expression("C content (mg C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") +
  theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g1

rt1=getAccumulatedRelease(M1$SoilRmodel)
rt1=as.data.frame(rt1)
rt1$time = seq(from = 1, to = 182, length.out=500)
rt1$total <- rt1$V1 + rt1$V2

g2= ggplot(rt1, aes(x = time)) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_line(aes(y = V1, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "MAOM"), linewidth = 0.8) +
  geom_point(data = rt, aes(x = time, y = rt), size = 2, color = "black", shape = 19) +
  geom_errorbar(data = rt, aes(x = time, ymin = rt - sd, ymax = rt + sd), 
                width = 0.1, color = "black") +
  scale_color_manual(values = c("C total" = "black", 
                                "POM" = "#8a3c4c", "MAOM" = "#0e2f4e")) +
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


saveRDS(g1, "./outputs/C_L05_2ps.rds")
saveRDS(g2, "./outputs/Resp_L05_2ps.rds")
saveRDS(M1, "./2_pool_model/L05_2ps.rds")

#################################################################################################


library(writexl)
C05_resp = readRDS("./2_pool_model/C05_2ps_resp.rds") # C respiration
C05_pools = readRDS("./2_pool_model/C05_2ps_2.rds") # C contents
L05_resp = readRDS("./2_pool_model/L05_2ps_resp.rds") # C respiration
L05_pools = readRDS("./2_pool_model/L05_2ps_2.rds") # C contents

outputs <- data.frame(
  Object = c("C05_resp", "C05_pools", "L05_resp", "L05_pools"),
  
  k1 = c(
    C05_resp$FMEmodel$par[1],
    C05_pools$FMEmodel$par[1],
    L05_resp$FMEmodel$par[1],
    L05_pools$FMEmodel$par[1]
  ),
  
  k2 = c(
    C05_resp$FMEmodel$par[2],
    C05_pools$FMEmodel$par[2],
    L05_resp$FMEmodel$par[2],
    L05_pools$FMEmodel$par[2]
  ),
  
  a21 = c(
    C05_resp$FMEmodel$par[3],
    C05_pools$FMEmodel$par[3],
    L05_resp$FMEmodel$par[3],
    L05_pools$FMEmodel$par[3]
  ),
  
  ms = c(
    C05_resp$FMEmodel$ms,
    C05_pools$FMEmodel$ms,
    L05_resp$FMEmodel$ms,
    L05_pools$FMEmodel$ms
  ),
  
  AIC = c(
    C05_resp[["AIC"]],
    C05_pools[["AIC"]],
    L05_resp[["AIC"]],
    L05_pools[["AIC"]]
  ),
  
  meanTT = c(
    C05_resp$TT$meanTransitTime,
    C05_pools$TT$meanTransitTime,
    L05_resp$TT$meanTransitTime,
    L05_pools$TT$meanTransitTime
  ),
  
  medianTT = c(
    C05_resp$TT$quantiles[2],
    C05_pools$TT$quantiles[2],
    L05_resp$TT$quantiles[2],
    L05_pools$TT$quantiles[2]
  )
)

write_xlsx(outputs, "./2_pool_model/outputs2ps.xlsx")
