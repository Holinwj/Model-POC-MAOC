
library(SoilR)
library(FME)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)

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
  pred_rt   <- approx(x = bestMod$time, y = bestMod$rt,   xout = rt$time)$y
  pred_pom  <- approx(x = bestMod$time, y = bestMod$pom,  xout = pom$time)$y
  pred_maom <- approx(x = bestMod$time, y = bestMod$maom, xout = maom$time)$y
  rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
  rmse_rt   <- rmse(rt$rt, pred_rt)
  rmse_pom  <- rmse(pom$pom, pred_pom)
  rmse_maom <- rmse(maom$maom, pred_maom)
  AIC=((n+2*(length(Fit$par)+1))/(n))+log(Fit$ms) # formula from Shumway & Stoffer
  AICc=log(Fit$ms)+((n+length(Fit$par)+1)/max(0,n-length(Fit$par)-2)) # for small sample size
  if(AICc == Inf) stop("Time series is too short. No degrees of freedom")
  SoilRmodel=SoilR::TwopSeriesModel(t=tt,ks=Fit$par[1:2], a21=Fit$par[1]*Fit$par[3],  
                                    C0=c(initial[1,2]*initial[3,2], initial[2,2]*initial[3,2]), In=0)
  A=matrix(c(-Fit$par[1],Fit$par[3]*Fit$par[1],0,-Fit$par[2]),ncol=2)
  u=matrix(c(initial[1,2]*initial[3,2], initial[2,2]*initial[3,2]),ncol=1)
  TT=transitTime(A,u)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TT=TT, AIC=AIC, AICc=AICc,
              rmse_rt=rmse_rt, rmse_pom=rmse_pom, rmse_maom=rmse_maom))
}

# Nelder-Mead    /    Marq

#0-5 cm depth control samples 

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
M1=twopsFit(rt, initial, inipars=c(0.01, 0.2, 0.1))   #inipars=c(0.01, 0.2, 0.1)) Good fit resp
t1=M1$SoilRmodel@times             #inipars=c(0.005, 0.002, 0.01)) Good fit C pool               
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
saveRDS(g1, "./outputs/C_C05_2ps_resp.rds")
saveRDS(g2, "./outputs/Resp_C05_2ps_resp.rds")

saveRDS(M1, "./2_pool_model/C05_2ps_resp.rds")
ggsave("./Outputs/C05_2ps_resp.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

saveRDS(M1, "./2_pool_model/C05_2ps_2.rds")
ggsave("./Outputs/C05_2ps.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

model_results <- data.frame(
  AIC = M1$AIC,
  AICc = M1$AICc,
  RMSE_rt = M1$rmse_rt,
  RMSE_pom = M1$rmse_pom,
  RMSE_maom = M1$rmse_maom)

#################################################################################################

#15-30 cm depth control samples 

rt=na.omit(data[-c(1:3),c(1,3,7)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./2_pool_model/data_C.xlsx")
pom = data_C[,c(1,4:5)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,12:13)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = as.data.frame(data[c(1:3),c(1,3)])

# fit the two pool series model with known initial C values 
M1=twopsFit(rt, initial, inipars=c(0.001, 0.005, 0.01))
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


saveRDS(g1, "./outputs/C_C1530_2ps.rds")
saveRDS(g2, "./outputs/Resp_C1530_2ps.rds")

saveRDS(M1, "./2_pool_model/C1530_2ps.rds")
ggsave("./Outputs/C1530_2ps.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

model_results <- data.frame(
  AIC = M1$AIC,
  AICc = M1$AICc,
  RMSE_rt = M1$rmse_rt,
  RMSE_pom = M1$rmse_pom,
  RMSE_maom = M1$rmse_maom)
#################################################################################################

#0-5 cm depth litter treatment samples 

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
M1=twopsFit(rt, initial, inipars=c(0.001, 0.0001, 0.0001))
t1=M1$SoilRmodel@times    #inipars=c(0.001, 0.0001, 0.0001)) Fit C pool
C1=getC(M1$SoilRmodel)        #inipars=c(0.03, 0.025, 0.07)) Fit resp      
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
saveRDS(g1, "./outputs/C_L05_2ps_resp.rds")
saveRDS(g2, "./outputs/Resp_L05_2ps_resp.rds")

saveRDS(M1, "./2_pool_model/L05_2ps_resp.rds")
ggsave("./Outputs/L05_2ps_resp.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

saveRDS(M1, "./2_pool_model/L05_2ps.rds")
ggsave("./Outputs/L05_2ps.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

model_results <- data.frame(
  AIC = M1$AIC,
  AICc = M1$AICc,
  RMSE_rt = M1$rmse_rt,
  RMSE_pom = M1$rmse_pom,
  RMSE_maom = M1$rmse_maom)
#################################################################################################

#15-30 cm depth litter treatment samples 

rt=na.omit(data[-c(1:3),c(1,5,9)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./2_pool_model/data_C.xlsx")
pom = data_C[,c(1,8:9)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,16:17)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = as.data.frame(data[c(1:3),c(1,5)])

# fit the two pool series model with known initial C values 
M1=twopsFit(rt, initial, inipars=c(0.03, 0.01, 0.01)) 
t1=M1$SoilRmodel@times    #inipars=c(0.03, 0.01, 0.01)) Fit resp
C1=getC(M1$SoilRmodel)               # inipars=c(0.001, 0.0001, 0.001))  Fit C pools
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

saveRDS(g1, "./outputs/C_L1530_2ps.rds")
saveRDS(g2, "./outputs/Resp_L1530_2ps.rds")
saveRDS(g1, "./outputs/C_L1530_2ps_resp.rds")
saveRDS(g2, "./outputs/Resp_L1530_2ps_resp.rds")

saveRDS(M1, "./2_pool_model/L1530_2ps.rds")
ggsave("./Outputs/L1530_2ps.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

saveRDS(M1, "./2_pool_model/L1530_2ps_resp.rds")
ggsave("./Outputs/L1530_2ps_resp.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

model_results <- data.frame(
  AIC = M1$AIC,
  AICc = M1$AICc,
  RMSE_rt = M1$rmse_rt,
  RMSE_pom = M1$rmse_pom,
  RMSE_maom = M1$rmse_maom)
#########################################################################################

library(writexl)
C05 = readRDS("./2_pool_model/C05_2ps.rds")
C05_2 = readRDS("./2_pool_model/C05_2ps_2.rds")
C1530 = readRDS("./2_pool_model/C1530_2ps.rds")
C1530_2 = readRDS("./2_pool_model/C1530_2ps_2.rds")
L05 = readRDS("./2_pool_model/L05_2ps.rds")
L1530 = readRDS("./2_pool_model/L1530_2ps.rds")

outputs <- data.frame(
  Object = c("C05", "C05_2", "C1530", "C1530_2", "L05", "L1530"),
  k1 = c(C05$FMEmodel$par[1], C05_2$FMEmodel$par[1], C1530$FMEmodel$par[1], C1530_2$FMEmodel$par[1], L05$FMEmodel$par[1], L1530$FMEmodel$par[1]),
  k2 = c(C05$FMEmodel$par[2], C05_2$FMEmodel$par[2], C1530$FMEmodel$par[2], C1530_2$FMEmodel$par[2], L05$FMEmodel$par[2], L1530$FMEmodel$par[2]),
  a21 = c(C05$FMEmodel$par[3], C05_2$FMEmodel$par[3], C1530$FMEmodel$par[3], C1530_2$FMEmodel$par[3], L05$FMEmodel$par[3], L1530$FMEmodel$par[3]),
  ms = c(C05$FMEmodel$ms, C05_2$FMEmodel$ms, C1530$FMEmodel$ms, C1530_2$FMEmodel$ms, L05$FMEmodel$ms, L1530$FMEmodel$ms),
  AIC =(c(C05[["AIC"]], C05_2[["AIC"]], C1530[["AIC"]],C1530_2[["AIC"]], L05[["AIC"]], L1530[["AIC"]])),
  AICc =(c(C05[["AICc"]], C05_2[["AICc"]], C1530[["AICc"]],C1530_2[["AICc"]], L05[["AICc"]], L1530[["AICc"]])),
  meanTT =(c(C05$TT$meanTransitTime,C05_2$TT$meanTransitTime,C1530$TT$meanTransitTime,C1530_2$TT$meanTransitTime,L05$TT$meanTransitTime,L1530$TT$meanTransitTime)),
  medianTT = (c(C05$TT$quantiles[2], C05_2$TT$quantiles[2], C1530$TT$quantiles[2], C1530_2$TT$quantiles[2], L05$TT$quantiles[2], L1530$TT$quantiles[2]))
  )

write_xlsx(outputs, "./2_pool_model/outputs2ps.xlsx")
