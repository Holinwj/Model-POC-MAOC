
# Three pools series model with a31# 

library(SoilR)
library(FME)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)

data = read_excel("./3_pool_model/data.xlsx")

ThreepSeriesLitterModel = function (t, ks, a21, a32, a31, C0, In, xi = 1, 
                                    solver = deSolve.lsoda.wrapper, pass = FALSE) 
{
  t_start = min(t)
  t_end = max(t)
  if (length(ks) != 3) 
    stop("ks must be of length = 3")
  if (length(C0) != 3) 
    stop("the vector with initial conditions must be of length = 3")
  if (length(In) == 1) {
    inputFluxes = BoundInFluxes(function(t) {
      matrix(nrow = 3, ncol = 1, c(In, 0, 0))
    }, t_start, t_end)
  }
  if (inherits(In, "data.frame")) {
    x = In[, 1]
    y = In[, 2]
    inputFlux = splinefun(x, y)
    inputFluxes = BoundInFluxes(function(t) {
      matrix(nrow = 3, ncol = 1, c(inputFlux(t), 0, 0))
    }, min(x), max(x))
  }
  A = -1 * abs(diag(ks))
  A[2, 1] = a21
  A[3, 2] = a32
  A[3, 1] = a31
  if (length(xi) == 1) 
    fX = function(t) {
      xi
    }
  if (inherits(xi, "data.frame")) {
    X = xi[, 1]
    Y = xi[, 2]
    fX = splinefun(X, Y)
  }
  Af = BoundLinDecompOp(function(t) {
    fX(t) * A
  }, t_start, t_end)
  Mod = GeneralModel(t = t, A = Af, ivList = C0, inputFluxes = inputFluxes, 
                     pass = pass)
  return(Mod)
}

threepFit=function(timeSeries, initial, inipars){
  complete=data.frame(time=timeSeries[complete.cases(timeSeries),1],Ct=timeSeries[complete.cases(timeSeries),2])
  n=nrow(complete)
  tt=seq(from=0, to=unlist(tail(complete[,1],1)), length.out = 500)
  
  Func=function(pars){
  mod=ThreepSeriesLitterModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5],
  a31=pars[1]*pars[6], C0=c(initial[1,2]*initial[4,2], initial[2,2]*initial[4,2],initial[3,2]*initial[4,2]), In=0)
  Ct=SoilR::getC(mod)
  rt = SoilR::getAccumulatedRelease(mod)
  return(data.frame(time=tt, rt=rowSums(rt), Ct=rowSums(Ct),litter=Ct[,1], pom = Ct[, 2], maom = Ct[, 3]))
}
  
  costFunc=function(pars){
    output=Func(pars)
    cost1 = modCost(model = output, obs = as.data.frame(rt), x="time", err = "sd" )
    cost2 = modCost(model = output, obs = as.data.frame(pom), x="time", err = "sd", cost = cost1)
    return(modCost(model=output, obs=as.data.frame(maom), x="time", err = "sd", cost = cost2)) 
  }
  
  Fit=modFit(f=costFunc, p=inipars, method="Nelder-Mead", lower=c(0,0,0,0,0,0), upper=c(Inf, Inf, Inf,1,1,1))
  bestMod=Func(pars=Fit$par)
  #RMSE
  rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
  rmse_rt   <- rmse(rt$rt, approx(x = bestMod$time, y = bestMod$rt,   xout = rt$time)$y)
  rmse_pom  <- rmse(pom$pom, approx(x = bestMod$time, y = bestMod$pom,  xout = pom$time)$y)
  rmse_maom <- rmse(maom$maom, approx(x = bestMod$time, y = bestMod$maom, xout = maom$time)$y)
  
  #AIC 
  AIC=((n+2*(length(Fit$par)+1))/(n))+log(Fit$ms) # formula from Shumway & Stoffer
  AICc=log(Fit$ms)+((n+length(Fit$par)+1)/max(0,n-length(Fit$par)-2)) # for small sample size
  if(AICc == Inf) stop("Time series is too short. No degrees of freedom")
  SoilRmodel=ThreepSeriesLitterModel(t=tt,ks=Fit$par[1:3], a21=Fit$par[1]*Fit$par[4], 
  a32=Fit$par[2]*Fit$par[5],a31=Fit$par[1]*Fit$par[6], C0=c(initial[1,2]*initial[4,2], 
  initial[2,2]*initial[4,2],initial[3,2]*initial[4,2]), In=0)
  A=matrix(c(-Fit$par[1],Fit$par[4]*Fit$par[1],Fit$par[1]*Fit$par[6],0,-Fit$par[2],Fit$par[2]*Fit$par[5],0,0,-Fit$par[3]),ncol=3)
  u=matrix(c(initial[1,2]*initial[4,2], initial[2,2]*initial[4,2], initial[3,2]*initial[4,2]),ncol=1)
  TT=transitTime(A,u)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TT=TT, AIC=AIC, AICc=AICc,
            rmse_rt=rmse_rt, rmse_pom=rmse_pom, rmse_maom=rmse_maom))
}
# Nelder-Mead & Marq


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_model/data_C.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])

# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.3, 0.008, 0.08, 0.03, 0.02, 0.03))
t1=M1$SoilRmodel@times     
C1=getC(M1$SoilRmodel)               
M1$FMEmodel$par         
M1$FMEmodel$ms         

rt1=getAccumulatedRelease(M1$SoilRmodel)

t1 = as.data.frame(t1)
t1$time = seq(from = 1, to = 182, length.out = 500)
C1 = as.data.frame(C1)
C1$time = seq(from = 1, to = 182, length.out=500)
C1$total <- C1$V1 + C1$V2 + C1$V3

g1 = ggplot(C1, aes(x = time)) +
  geom_line(aes(y = V1, color = "litter"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V3, color = "MAOM"), linewidth = 0.8) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_point(data = data.frame(time = maom$time, C = maom$maom), 
             aes(x = time, y = C), size=2, color = "#0e2f4e", shape = 19) +
  geom_point(data = data.frame(time = pom$time, C = pom$pom), 
             aes(x = time, y = C), size=2, color = "#8a3c4c", shape = 19) +
  geom_point(data = data.frame(time = litter$time, C = litter$litter), 
             aes(x = time, y = C), size=2,color = "#e89c8f", shape = 19) +
  geom_errorbar(data = pom, aes(x = time, ymin = pom - sd, ymax = pom + sd), 
                width = 0.1, color = "#8a3c4c") +
  geom_errorbar(data = maom, aes(x = time, ymin = maom - sd, ymax = maom + sd), 
                width = 0.1, color = "#0e2f4e") +
  scale_color_manual(values = c("C total" = "black", "litter" = "#e89c8f", "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "Litter C", "POC", "MAOC")) +
  labs(y = expression("C content (mg C g soil "^"-1"*")"), 
                      color = "", x = "Time (days)") +
  scale_y_continuous(limits = c(0, 105))+
  theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g1

rt1=as.data.frame(rt1)
rt1$time = seq(from = 1, to = 182, length.out=500)
rt1$total <- rt1$V1 + rt1$V2 + rt1$V3

g2= ggplot(rt1, aes(x = time)) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_line(aes(y = V1, color = "litter"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V3, color = "MAOM"), linewidth = 0.8) +
  geom_point(data = rt, aes(x = time, y = rt), size = 2, color = "black", shape = 19) +
  geom_errorbar(data = rt, aes(x = time, ymin = rt - sd, ymax = rt + sd), 
                width = 0.1, color = "black") +
  scale_color_manual(values = c("C total" = "black", "litter" = "#e89c8f", 
                                "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/C_L05_3ps.rds")
saveRDS(g2, "./outputs/Resp_L05_3ps.rds")
saveRDS(M1, "./3_pool_model/L05_3ps.rds")
ggsave("./outputs/L05_3ps.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

model_results <- data.frame(
  AIC = M1$AIC,
  AICc = M1$AICc,
  RMSE_rt = M1$rmse_rt,
  RMSE_pom = M1$rmse_pom,
  RMSE_maom = M1$rmse_maom)
#################################################################################################

#15-30 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,5,9)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_model/data_C.xlsx")
pom = data_C[,c(1,8:9)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,16:17)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,24:25)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = as.data.frame(data[c(1:4),c(1,5)])

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.3, 0.08, 0.08, 0.04, 0.002, 0.001))
t1=M1$SoilRmodel@times    #inipars=c(0.3, 0.08, 0.08, 0.04, 0.002, 0.001)) Good Ct and "good" rt Nelder
C1=getC(M1$SoilRmodel)    #inipars=c(0.01, 0.01, 0.001, 0.01, 0.1, 0.01)) Bad ct good rt Nelder 
M1$FMEmodel$par           # inipars=c(0.2, 0.005, 0.003, 0.08, 0.05, 0.01)) Bad Ct good rt
M1$FMEmodel$ms

rt1=getAccumulatedRelease(M1$SoilRmodel)

t1 = as.data.frame(t1)
t1$time = seq(from = 1, to = 182, length.out = 500)
C1 = as.data.frame(C1)
C1$time = seq(from = 1, to = 182, length.out=500)
C1$total <- C1$V1 + C1$V2 + C1$V3

g1 = ggplot(C1, aes(x = time)) +
  geom_line(aes(y = V1, color = "litter"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V3, color = "MAOM"), linewidth = 0.8) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_point(data = data.frame(time = maom$time, C = maom$maom), 
             aes(x = time, y = C), size=2, color = "#0e2f4e", shape = 19) +
  geom_point(data = data.frame(time = pom$time, C = pom$pom), 
             aes(x = time, y = C), size=2, color = "#8a3c4c", shape = 19) +
  geom_point(data = data.frame(time = litter$time, C = litter$litter), 
             aes(x = time, y = C), size=2,color = "#e89c8f", shape = 19) +
  geom_errorbar(data = pom, aes(x = time, ymin = pom - sd, ymax = pom + sd), 
                width = 0.1, color = "#8a3c4c") +
  geom_errorbar(data = maom, aes(x = time, ymin = maom - sd, ymax = maom + sd), 
                width = 0.1, color = "#0e2f4e") +
  scale_color_manual(values = c("C total" = "black", "litter" = "#e89c8f", "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "Litter C", "POC", "MAOC")) +
  labs(y = expression("C content (mg C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") +
  scale_y_continuous(limits = c(0, 48))+
  theme_classic() + theme(
    axis.text.x = element_text(size = 10, color="black"),
    axis.text.y = element_text(size = 10, color="black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g1


rt1=as.data.frame(rt1)
rt1$time = seq(from = 1, to = 182, length.out=500)
rt1$total <- rt1$V1 + rt1$V2 + rt1$V3

g2= ggplot(rt1, aes(x = time)) +
  geom_line(aes(y = total, color = "C total"), linewidth = 0.8) +
  geom_line(aes(y = V1, color = "litter"), linewidth = 0.8) +
  geom_line(aes(y = V2, color = "POM"), linewidth = 0.8) +
  geom_line(aes(y = V3, color = "MAOM"), linewidth = 0.8) +
  geom_point(data = rt, aes(x = time, y = rt), size = 2, color = "black", shape = 19) +
  geom_errorbar(data = rt, aes(x = time, ymin = rt - sd, ymax = rt + sd), 
                width = 0.1, color = "black") +
  scale_color_manual(values = c("C total" = "black", "litter" = "#e89c8f", 
                                "POM" = "#8a3c4c", "MAOM" = "#0e2f4e"),
                     labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 8))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/C_L1530_3ps.rds")
saveRDS(g2, "./outputs/Resp_L1530_3ps.rds")
saveRDS(M1, "./3_pool_model/L1530_3ps.rds")
ggsave("./outputs/L1530_3ps.jpg", plot = graph, device = "jpeg", width = 6, height = 6, dpi = 300)

model_results <- data.frame(
  AIC = M1$AIC,
  AICc = M1$AICc,
  RMSE_rt = M1$rmse_rt,
  RMSE_pom = M1$rmse_pom,
  RMSE_maom = M1$rmse_maom)
#########################################################################################################


library(writexl)
C05 = readRDS("./3_pool_model/C05_3ps.rds")
C1530 = readRDS("./3_pool_model/C1530_3ps.rds")
L05 = readRDS("./3_pool_model/L05_3ps.rds")
L1530 = readRDS("./3_pool_model/L1530_3ps.rds")

outputs <- data.frame(
  Object = c("C05", "C1530", "L05", "L1530"),
  k1 = c(C05$FMEmodel$par[1], C1530$FMEmodel$par[1], L05$FMEmodel$par[1], L1530$FMEmodel$par[1]),
  k2 = c(C05$FMEmodel$par[2], C1530$FMEmodel$par[2], L05$FMEmodel$par[2], L1530$FMEmodel$par[2]),
  k3 = c(C05$FMEmodel$par[3], C1530$FMEmodel$par[3], L05$FMEmodel$par[3], L1530$FMEmodel$par[3]),
  a21 = c(C05$FMEmodel$par[4], C1530$FMEmodel$par[4], L05$FMEmodel$par[4], L1530$FMEmodel$par[4]),
  a32 = c(C05$FMEmodel$par[5], C1530$FMEmodel$par[5], L05$FMEmodel$par[5], L1530$FMEmodel$par[5]),
  a31 = c(C05$FMEmodel$par[6], C1530$FMEmodel$par[6], L05$FMEmodel$par[6], L1530$FMEmodel$par[6]),
  ms = c(C05$FMEmodel$ms, C1530$FMEmodel$ms, L05$FMEmodel$ms, L1530$FMEmodel$ms),
  AIC =(c(C05[["AIC"]], C1530[["AIC"]], L05[["AIC"]], L1530[["AIC"]])),
  AICc =(c(C05[["AICc"]], C1530[["AICc"]], L05[["AICc"]], L1530[["AICc"]])),
  meanTT =(c(C05$TT$meanTransitTime,C1530$TT$meanTransitTime,L05$TT$meanTransitTime,L1530$TT$meanTransitTime)),
  medianTT = (c(C05$TT$quantiles[2], C1530$TT$quantiles[2], L05$TT$quantiles[2], L1530$TT$quantiles[2]  ))
  )

write_xlsx(outputs, "./3_pool_model/outputs3ps.xlsx")
