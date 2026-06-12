library(SoilR)
library(FME)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(writexl)
library(furrr)
library(purrr)



# 3-pools model -15% of POC as litter C

data = read_excel("./3_pool_assumptions/data_15.xlsx")

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
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TT=TT, AIC=AIC,
              costFunc   = costFunc))
}

# Control 0-5 cm 

rt=na.omit(data[-c(1:4),c(1,2,6)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_15.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(1234)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.92809341, 0.42115686, 0.30142211, 0.69913124, 0.70244682, 0.03197319))
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 12))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


saveRDS(g1, "./outputs/model_assumptions/C_C05_15litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_15litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_15litter.rds")

saveRDS(g1, "./outputs/model_assumptions/g1_controldata_withlitterpars.rds")
saveRDS(g2, "./outputs/model_assumptions/g2_controldata_withlitterpars.rds")

####################################################################


########################################################################################

# 15% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_15.xlsx")

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


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_15.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.11930840, 0.16870946, 0.60112933, 0.43024695, 0.71055360, 0.01699834)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_L05_15litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_15litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_15litter.rds")

######################################################################################


# 3-pools model -5% of POC as litter C

data = read_excel("./3_pool_assumptions/data_5.xlsx")

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
data_C = read_excel("./3_pool_assumptions/data_C_5.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.6292380, 0.9158102, 0.2683922, 0.0264171, 0.4133907, 0.8874679))
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_C05_5litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_5litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_5litter.rds")

######################################################################################


# 3-pools model -10% of POC as litter C

data = read_excel("./3_pool_assumptions/data_10.xlsx")

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
data_C = read_excel("./3_pool_assumptions/data_C_10.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.06516906, 0.64272948, 0.21793446, 0.02871518, 0.65517999, 0.03304440))
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


M1= readRDS("./3_pool_assumptions/C05_3ps_10litter.rds")
M1$FMEmodel$ms # 0.3990198
saveRDS(g1, "./outputs/model_assumptions/C_C05_10litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_10litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_10litter.rds")

####################################################################

######################################################################################



# 3-pools model -20% of POC as litter C

data = read_excel("./3_pool_assumptions/data_20.xlsx")

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
data_C = read_excel("./3_pool_assumptions/data_C_20.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.74436459, 0.31378990, 0.59009145, 0.01792305, 0.26332320, 0.45614771))
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_C05_20litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_20litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_20litter.rds")

####################################################################


####################################################################

# 3-pools model -25% of POC as litter C

data = read_excel("./3_pool_assumptions/data_25.xlsx")

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
data_C = read_excel("./3_pool_assumptions/data_C_25.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.04436767, 0.65142417, 0.18200370, 0.01874155, 0.84480997, 0.43366116))
t1=M1$SoilRmodel@times    
C1=getC(M1$SoilRmodel)               
M1$FMEmodel$par          
M1$FMEmodel$ms # 0.4179912

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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_C05_25litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_25litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_25litter.rds")

####################################################################

# 3-pools model -30% of POC as litter C

data = read_excel("./3_pool_assumptions/data_30.xlsx")

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
data_C = read_excel("./3_pool_assumptions/data_C_30.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:24]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.10933215, 0.20560047, 0.01100966, 0.62632208, 0.86121363, 0.02338490))
t1=M1$SoilRmodel@times    
C1=getC(M1$SoilRmodel)               
M1$FMEmodel$par          
M1$FMEmodel$ms # 0.3028235

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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 5.5))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_C05_30litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_30litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_30litter.rds")

####################################################################


# 3-pools model -40% of POC as litter C

data = read_excel("./3_pool_assumptions/data_40.xlsx")

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
data_C = read_excel("./3_pool_assumptions/data_C_40.xlsx")

pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,18)]
colnames(litter)=c("time", "litter")

# get initial values
initial = data.frame(data[c(1:4),c(1,2)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models

# fit the two pool series model with known initial C values 
M1=threepFit(rt, initial, inipars=c(0.04436767, 0.65142417, 0.18200370, 0.01874155, 0.84480997, 0.43366116))
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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

saveRDS(g1, "./outputs/model_assumptions/C_C05_40litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_C05_40litter.rds")
saveRDS(M1, "./3_pool_assumptions/C05_3ps_40litter.rds")

############################################################################################################
############################################################################################################

library(purrr)
library(dplyr)
library(writexl)

litters <- c(1, 5, 10, 15, 20, 30, 40)

outputs <- map_dfr(litters, function(lit) {
  
  file_path <- paste0("./3_pool_assumptions/C05_3ps_", lit, "litter.rds")
  
  mod <- readRDS(file_path)
  
  data.frame(
    Object   = paste0("L05_", lit),
    k1       = mod$FMEmodel$par[1],
    k2       = mod$FMEmodel$par[2],
    k3       = mod$FMEmodel$par[3],
    a21      = mod$FMEmodel$par[1] * mod$FMEmodel$par[4],
    a32      = mod$FMEmodel$par[2] * mod$FMEmodel$par[5],
    a31      = mod$FMEmodel$par[1] * mod$FMEmodel$par[6],
    ms       = mod$FMEmodel$ms,
    AIC      = mod$AIC,
    meanTT   = mod$TT$meanTransitTime,
    medianTT = mod$TT$quantiles[2]
  )
})

write_xlsx(outputs, "./3_pool_assumptions/outputs3ps_C_litter.xlsx")
M1$FMEmodel$par[4]*M1$FMEmodel$par[1]*365

########################################################################################

# 0% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data.xlsx")

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


rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])

niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.06516906, 0.64272948, 0.21793446, 0.02871518, 0.65517999, 0.03304440)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


saveRDS(g1, "./outputs/model_assumptions/C_L05_0litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_0L05_litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_0litter.rds")

########################################################################################

# 5% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_5.xlsx")

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


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_5.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.06516906, 0.64272948, 0.21793446, 0.02871518, 0.65517999, 0.03304440)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_L05_5litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_5litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_5litter.rds")

########################################################################################

# 10% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_10.xlsx")

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


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_10.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.11930840, 0.16870946, 0.60112933, 0.43024695, 0.71055360, 0.01699834)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_L05_10litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_10litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_10litter.rds")


# 20% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_20.xlsx")

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


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_20.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.11930840, 0.16870946, 0.60112933, 0.43024695, 0.71055360, 0.01699834)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)

saveRDS(g1, "./outputs/model_assumptions/C_L05_20litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_20litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_20litter.rds")

########################################################################################

# 25% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_25.xlsx")

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


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_25.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.06516906, 0.64272948, 0.21793446, 0.02871518, 0.65517999, 0.03304440)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


saveRDS(g1, "./outputs/model_assumptions/C_L05_25litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_25litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_25litter.rds")

########################################################################################
# 30% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_30.xlsx")

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


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_30.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.05660153, 0.51178771, 0.48631939, 0.30130161, 0.31423121, 0.14053317)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


saveRDS(g1, "./outputs/model_assumptions/C_L05_30litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_30litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_30litter.rds")


########################################################################################

########################################################################################
# 40% POM as litter for litter-addition soils

data = read_excel("./3_pool_assumptions/data_40.xlsx")

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
  u=matrix(c(initia=l[1,2]*initial[4,2], initial[2,2]*initial[4,2], initial[3,2]*initial[4,2]),ncol=1)
  TT=transitTime(A,u)
  return(list(FMEmodel=Fit, SoilRmodel=SoilRmodel, TT=TT, AIC=AIC))
}


#0-5 cm depth litter treatment samples 

rt=na.omit(data[-c(1:4),c(1,4,8)])
rt$time=as.numeric(rt$time)
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

# get pom and maom values for cost function
data_C = read_excel("./3_pool_assumptions/data_C_40.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

# get initial values
initial = data.frame(data[c(1:4),c(1,4)])


niter <- 200

plan(multisession, workers = parallel::detectCores())

set.seed(123)

ranges <- list(
  k1  = c(0, 1),
  k2  = c(0, 1),
  k3  = c(0, 1),
  a21 = c(0,1),
  a32 = c(0, 1),
  a31 = c(0, 1)
)

random_inipars <- function(ranges) {
  c(
    runif(1, ranges$k1[1],  ranges$k1[2]),
    runif(1, ranges$k2[1],  ranges$k2[2]),
    runif(1, ranges$k3[1],  ranges$k3[2]),
    runif(1, ranges$a21[1], ranges$a21[2]),
    runif(1, ranges$a32[1], ranges$a32[2]),
    runif(1, ranges$a31[1], ranges$a31[2])
  )
}

results <- future_map(
  1:niter,
  function(i) {
    inipars_try <- random_inipars(ranges)
    
    out <- try({
      fit_i <- threepFit(rt, initial, inipars_try)
      list(
        fit     = fit_i,
        ms      = fit_i$FMEmodel$ms,
        AIC     = fit_i$AIC,
        inipars = inipars_try
      )
    }, silent = TRUE)
    
    if (inherits(out, "try-error")) return(NULL)
    return(out)
  },
  .options = furrr::furrr_options(seed = TRUE)   
)

results <- compact(results)

best_models <- {
  AIC_values <- map_dbl(results, "AIC")
  top50 <- results[order(AIC_values)][1:2]
  
  map(top50, ~ list(
    ms      = .x$ms,
    inipars = .x$inipars,
    pars    = .x$fit$FMEmodel$par
  ))
}

best_models


# fit the two pool series model with known initial C values  
M1=threepFit(rt, initial, inipars=c(0.11372618, 0.35837425, 0.06984211, 0.26851796, 0.45169216, 0.23250124)) 
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
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
  scale_color_manual(
    values = c("C total" = "black", "litter"  = "#e89c8f",
               "POM"     = "#8a3c4c", "MAOM"    = "#0e2f4e"),
    breaks = c("C total", "litter", "POM", "MAOM"),
    labels = c("Total C", "Litter C", "POC", "MAOC")) +
  scale_y_continuous(limits = c(0, 13.2))+
  labs(y = expression("Respired C (mg CO2-C g soil "^"-1"*")"), 
       color = "", x = "Time (days)") + theme_classic() + theme(
         axis.text.x = element_text(size = 10, color="black"),
         axis.text.y = element_text(size = 10, color="black"),
         panel.border = element_rect(color = "black", fill = NA, size = 0.5))
g2

graph = grid.arrange(g1, g2)


saveRDS(g1, "./outputs/model_assumptions/C_L05_40litter.rds")
saveRDS(g2, "./outputs/model_assumptions/Resp_L05_40litter.rds")
saveRDS(M1, "./3_pool_assumptions/L05_3ps_40litter.rds")


########################################################################################

library(purrr)
library(dplyr)
library(writexl)

litters <- c(0, 1, 5, 10, 15, 20, 30, 40)

outputs <- map_dfr(litters, function(lit) {
  
  file_path <- paste0("./3_pool_assumptions/L05_3ps_", lit, "litter.rds")
  
  mod <- readRDS(file_path)
  
  data.frame(
    Object   = paste0("L05_", lit),
    k1       = mod$FMEmodel$par[1],
    k2       = mod$FMEmodel$par[2],
    k3       = mod$FMEmodel$par[3],
    a21      = mod$FMEmodel$par[1] * mod$FMEmodel$par[4],
    a32      = mod$FMEmodel$par[2] * mod$FMEmodel$par[5],
    a31      = mod$FMEmodel$par[1] * mod$FMEmodel$par[6],
    ms       = mod$FMEmodel$ms,
    AIC      = mod$AIC,
    meanTT   = mod$TT$meanTransitTime,
    medianTT = mod$TT$quantiles[2]
  )
})

write_xlsx(outputs, "./3_pool_assumptions/outputs3ps_L_litter.xlsx")
