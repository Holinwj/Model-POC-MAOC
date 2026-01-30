library(SoilR)
library(FME)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(writexl)

# 3-pools model -10% of POC as litter C

data = read_excel("./3_pool_model_10litter/data.xlsx")

ThreepSeriesLitterModel = function (t, ks, a21, a32, a31, C0, In, xi = 1, solver = deSolve.lsoda.wrapper, 
                                    pass = FALSE) 
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
  
  AIC=((n+2*(length(Fit$par)+1))/(n))+log(Fit$ms) 
  SoilRmodel=ThreepSeriesLitterModel(t=tt,ks=Fit$par[1:3], a21=Fit$par[1]*Fit$par[4], 
                                     a32=Fit$par[2]*Fit$par[5],a31=Fit$par[1]*Fit$par[6], C0=c(initial[1,2]*initial[4,2], 
                                                                                               initial[2,2]*initial[4,2],initial[3,2]*initial[4,2]), In=0)
  A=matrix(c(-Fit$par[1],Fit$par[4]*Fit$par[1],Fit$par[1]*Fit$par[6],0,-Fit$par[2],Fit$par[2]*Fit$par[5],0,0,-Fit$par[3]),ncol=3)
  u=matrix(c(initial[1,2]*initial[4,2], initial[2,2]*initial[4,2], initial[3,2]*initial[4,2]),ncol=1)
  TT=transitTime(A,u)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TT=TT, AIC=AIC))
}

# Control 0-5 cm 

rt=na.omit(data[-c(1:4),c(1,2,6)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_model_10litter/data_C.xlsx")
pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.12, 0.04, 0.05, 0.05, 0.5, 0.3))
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
  scale_y_continuous(limits = c(0, 95))+
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
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/C_C05_10litter.rds")
saveRDS(g2, "./outputs/Resp_C05_10litter.rds")
saveRDS(M1, "./3_pool_model_10litter/C05_3ps_10litter.rds")

####################################################################

library(writexl)

C05 = readRDS("./3_pool_model_10litter/C05_3ps_10litter.rds")

outputs <- data.frame(
  Object   = "C05",
  k1       = C05$FMEmodel$par[1],
  k2       = C05$FMEmodel$par[2],
  k3       = C05$FMEmodel$par[3],
  a21      = C05$FMEmodel$par[4],
  a32      = C05$FMEmodel$par[5],
  a31      = C05$FMEmodel$par[6],
  ms       = C05$FMEmodel$ms,
  AIC      = C05[["AIC"]],
  meanTT   = C05$TT$meanTransitTime,
  medianTT = C05$TT$quantiles[2]
)

write_xlsx(outputs, "./3_pool_model_10litter/outputs3ps_10litter.xlsx")
