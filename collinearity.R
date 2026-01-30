library(SoilR)
library(FME)
library(dplyr)
library(readxl)

# 2-pool models

data_col = read_excel("./2_pool_model/collinearity/data_col.xlsx")

# Control soils

rt=(data_col[,c(1,2:3)]) # keep only time and respiration data
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

data_C = read_excel("./2_pool_model/collinearity/data_C.xlsx")
pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,10:11)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = data.frame(c(pom[1,2], maom[1,2]))
colnames(initial) = c("pom", "maom")

costFunc=function(pars){
  output=Func(pars)
  cost1 = modCost(model = output, obs = as.data.frame(rt), x="time", err = "sd" )
  cost2 = modCost(model = output, obs = as.data.frame(pom), x="time", err = "sd", cost = cost1)
  return(modCost(model=output, obs=as.data.frame(maom), x="time", err = "sd", cost = cost2)) 
}

tt=seq(from=0, to=unlist(tail(rt[,1],1)), length.out = 500)

inipars=c(0.01, 0.2, 0.1)

Func=function(pars){
  mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=c(unlist(initial[1,1]), unlist(initial[1,2])), In=0)
  Ct=SoilR::getC(mod)
  rt = SoilR::getAccumulatedRelease(mod)
  return(data.frame(time=tt, rt=rowSums(rt), Ct=rowSums(Ct), 
                    pom = Ct[, 1], maom = Ct[, 2]))
}

Sfun <- sensFun(costFunc, inipars)  
twops=collin(Sfun)

twops$collinearity

saveRDS(twops, "./2_pool_model/collinearity_2ps_control.rds")

######################################################################################

# Litter-addition treatment

rt=(data_col[,c(1,6:7)]) # keep only time and respiration data
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

data_C = read_excel("./2_pool_model/collinearity/data_C.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")

# get initial values
initial = data.frame(c(pom[1,2], maom[1,2]))
colnames(initial) = c("pom", "maom")

costFunc=function(pars){
  output=Func(pars)
  cost1 = modCost(model = output, obs = as.data.frame(rt), x="time", err = "sd" )
  cost2 = modCost(model = output, obs = as.data.frame(pom), x="time", err = "sd", cost = cost1)
  return(modCost(model=output, obs=as.data.frame(maom), x="time", err = "sd", cost = cost2)) 
}

tt=seq(from=0, to=unlist(tail(rt[,1],1)), length.out = 500)

inipars=c(0.03, 0.025, 0.07)

Func=function(pars){
  mod=SoilR::TwopSeriesModel(t=tt,ks=pars[1:2], a21=pars[1]*pars[3], C0=c(unlist(initial[1,1]), unlist(initial[1,2])), In=0)
  Ct=SoilR::getC(mod)
  rt = SoilR::getAccumulatedRelease(mod)
  return(data.frame(time=tt, rt=rowSums(rt), Ct=rowSums(Ct), 
                    pom = Ct[, 1], maom = Ct[, 2]))
}

Sfun <- sensFun(costFunc, inipars)  
twops=collin(Sfun)

twops$collinearity

saveRDS(twops, "./2_pool_model/collinearity_2ps_litter.rds")

######################################################################################

# 3-pool models 

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

data_col = read_excel("./3_pool_model/collinearity_3p/data_resp.xlsx")

rt=(data_col[,c(1,6:7)]) # keep only time and respiration data
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

data_C = read_excel("./3_pool_model/data_C.xlsx")
pom = data_C[,c(1,6:7)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,14:15)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,22:23)]
colnames(litter)=c("time", "litter", "sd")

initial = data.frame(c(litter[1,2], pom[1,2], maom[1,2]))
colnames(initial) = c("litter", "pom", "maom")

costFunc=function(pars){
  output=Func(pars)
  cost1 = modCost(model = output, obs = as.data.frame(rt), x="time", err = "sd" )
  cost2 = modCost(model = output, obs = as.data.frame(pom), x="time", err = "sd", cost = cost1)
  return(modCost(model=output, obs=as.data.frame(maom), x="time", err = "sd", cost = cost2)) 
}

tt=seq(from=0, to=unlist(tail(rt[,1],1)), length.out = 500)

inipars=c(0.3, 0.005, 0.08, 0.03, 0.02, 0.03)

Func=function(pars){
  mod=ThreepSeriesLitterModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5],
                              a31=pars[1]*pars[6], C0=c(unlist(initial[1,1]), unlist(initial[1,2]), unlist(initial[1,3])), In=0)
  Ct=SoilR::getC(mod)
  rt = SoilR::getAccumulatedRelease(mod)
  return(data.frame(time=tt, rt=rowSums(rt), Ct=rowSums(Ct), 
                    litter=Ct[,1], pom = Ct[, 2], maom = Ct[, 3]))
}

Sfun <- sensFun(costFunc, inipars)
threeps=collin(Sfun)

threeps$collinearity

saveRDS(threeps, "./3_pool_model/collinearity_3ps.rds")

##########################################################################3

# 3-pool model 10% of POC as litter C

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

data_col = read_excel("./3_pool_model_10litter/data_resp.xlsx")

rt=(data_col[,c(1,2:3)]) # keep only time and respiration data
colnames(rt)=c("time", "rt", "sd")
rt$sd[rt$sd <= 0] = 1e-9

data_C = read_excel("./3_pool_model_10litter/data_C_col.xlsx")
pom = data_C[,c(1,2:3)]
colnames(pom)=c("time", "pom", "sd")
maom = data_C[,c(1,4:5)] 
colnames(maom)=c("time", "maom", "sd")
litter = data_C[,c(1,6:7)]
colnames(litter)=c("time", "litter", "sd")

initial = data.frame(c(litter[1,2], pom[1,2], maom[1,2]))
colnames(initial) = c("litter", "pom", "maom")

costFunc=function(pars){
  output=Func(pars)
  cost1 = modCost(model = output, obs = as.data.frame(rt), x="time", err = "sd" )
  cost2 = modCost(model = output, obs = as.data.frame(pom), x="time", err = "sd", cost = cost1)
  return(modCost(model=output, obs=as.data.frame(maom), x="time", err = "sd", cost = cost2)) 
}

tt=seq(from=0, to=unlist(tail(rt[,1],1)), length.out = 500)

inipars=c(0.4, 0.05, 0.04, 0.05, 0.5, 0.3)

Func=function(pars){
  mod=ThreepSeriesLitterModel(t=tt,ks=pars[1:3], a21=pars[1]*pars[4], a32=pars[2]*pars[5],
                              a31=pars[1]*pars[6], C0=c(unlist(initial[1,1]), unlist(initial[1,2]), unlist(initial[1,3])), In=0)
  Ct=SoilR::getC(mod)
  rt = SoilR::getAccumulatedRelease(mod)
  return(data.frame(time=tt, rt=rowSums(rt), Ct=rowSums(Ct),litter=Ct[,1], pom = Ct[, 2], maom = Ct[, 3]))
}

Sfun <- sensFun(costFunc, inipars)
threeps=collin(Sfun)

threeps$collinearity

saveRDS(threeps, "./3_pool_model/collinearity_3ps_10litter.rds")

##########################################################################
# end