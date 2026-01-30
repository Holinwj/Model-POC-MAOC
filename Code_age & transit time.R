library(SoilR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


############################################### 1. Age and TT for part 1 ########################################### 
rm(list = ls())

###############################control 2-pool_fit respiration
k2s1=c(0.0139*365,0.00011*365)
A2s1=diag(-k2s1)
A2s1[2,1]=k2s1[1]*0.8318
u2s1=matrix(c(100,0),ncol=1)

###############################control 3-pool_10% litter
k3s1=c(0.0252*365,0.00033*365,0.00027*365)
A3s1=diag(-k3s1)
A3s1[2,1]=k3s1[1]*0.1734
A3s1[3,2]=k3s1[2]*0.9678
A3s1[3,1]=k3s1[1]*0.0310
u3s1=matrix(c(100,0,0),ncol=1)

############################### +litter 2-pool_fit respiration
k2s2=c(0.0168*365,0.0000296*365)
A2s2=diag(-k2s2)
A2s2[2,1]=k2s2[1]*0.6726
u2s2=matrix(c(100,0),ncol=1)

############################### +litter 3-pool
k3s2=c(0.01878*365,0.00027*365,0.000044*365)
A3s2=diag(-k3s2)
A3s2[2,1]=k3s2[1]*0.100003
A3s2[3,2]=k3s2[2]*0.02162
A3s2[3,1]=k3s2[1]*0.034561
u3s2=matrix(c(100,0,0),ncol=1)

##################################### Age and transit time
SA2s1=systemAge(a = seq(0, 20000),A=A2s1, u=u2s1)
SA2s2=systemAge(a = seq(0, 20000),A=A2s2, u=u2s2)
SA3s1=systemAge(a = seq(0, 20000),A=A3s1, u=u3s1)
SA3s2=systemAge(a = seq(0, 20000),A=A3s2, u=u3s2)

TT2s1=transitTime(a = seq(0, 20000),A=A2s1, u=u2s1)
TT2s2=transitTime(a = seq(0, 20000),A=A2s2, u=u2s2)
TT3s1=transitTime(a = seq(0, 20000),A=A3s1, u=u3s1)
TT3s2=transitTime(a = seq(0, 20000),A=A3s2, u=u3s2)

SA2s1$meanSystemAge
SA2s2$meanSystemAge
SA3s1$meanSystemAge
SA3s2$meanSystemAge

TT2s1$meanTransitTime
TT2s2$meanTransitTime
TT3s1$meanTransitTime
TT3s2$meanTransitTime

saveRDS(SA2s1, file = "./C_stability/SA2s1.rds")
saveRDS(SA2s2, file = "./C_stability/SA2s2.rds")
saveRDS(SA3s1, file = "./C_stability/SA3s1.rds")
saveRDS(SA3s2, file = "./C_stability/SA3s2.rds")
saveRDS(TT2s1, file = "./C_stability/TT2s1.rds")
saveRDS(TT2s2, file = "./C_stability/TT2s2.rds")
saveRDS(TT3s1, file = "./C_stability/TT3s1.rds")
saveRDS(TT3s2, file = "./C_stability/TT3s2.rds")

SA2s1 = readRDS("./C_stability/SA2s1.rds")
SA2s2 = readRDS("./C_stability/SA2s2.rds")
SA3s1 = readRDS("./C_stability/SA3s1.rds")
SA3s2 = readRDS("./C_stability/SA3s2.rds")

TT2s1 = readRDS("./C_stability/TT2s1.rds")
TT2s2 = readRDS("./C_stability/TT2s2.rds")
TT3s1 = readRDS("./C_stability/TT3s1.rds")
TT3s2 = readRDS("./C_stability/TT3s2.rds")


cage <- bind_rows(
  data.frame(
    age = seq_along(SA2s1$systemAgeDensity),
    density = SA2s1$systemAgeDensity,
    meanAge = SA2s1$meanSystemAge,
    model = "Control 2-pools"
  ),
  data.frame(
    age = seq_along(SA3s1$systemAgeDensity),
    density = SA3s1$systemAgeDensity,
    meanAge = SA3s1$meanSystemAge,
    model = "Control 3-pools"
  ),
  data.frame(
    age = seq_along(SA2s2$systemAgeDensity),
    density = SA2s2$systemAgeDensity,
    meanAge = SA2s2$meanSystemAge,
    model = "Litter 2-pools"
  ),
  data.frame(
    age = seq_along(SA3s2$systemAgeDensity),
    density = SA3s2$systemAgeDensity,
    meanAge = SA3s2$meanSystemAge,
    model = "Litter 3-pools"
  )
)

model_labels <- cage %>% 
  distinct(model, meanAge) %>%
  mutate(label = paste0(model, " = ", round(meanAge, 0), " yr")) %>%
  select(model, label) %>%     
  tibble::deframe()

colors = c("#17154FFF", "#BF3729FF", "#E69B00FF", "#355828FF")

C_age <- ggplot(cage, aes(x = age, y = density, color = model)) +
  geom_line(linewidth = 0.5) +
  geom_vline(
    aes(xintercept = meanAge, color = model),
    linetype = 5,
    linewidth = 0.5,
    show.legend = FALSE
  ) + 
  scale_color_manual(
    name = NULL,
    values = colors,
    labels = model_labels,
    guide = guide_legend(override.aes = list(linewidth = 1))
  ) +
  coord_cartesian(xlim = c(15, 350), ylim = c(-0.0005, 0.01)) +
  labs(
    x = "C age (years)",
    y = "Density function"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = c(0.715, 0.75),   
    legend.background = element_rect(
      fill = NA,
      color = NA
    ),
    legend.key = element_blank()
  )

C_age
#



get_age_tt <- function(TT) {
  if (!is.null(TT$age) && length(TT$age) > 0) {
    TT$age
  } else {
    seq_along(TT$transitTimeDensity)
  }
}

df_TT <- bind_rows(
  data.frame(
    age = get_age_tt(TT2s1),
    density = TT2s1$transitTimeDensity,
    meanTT = TT2s1$meanTransitTime,
    model = "Control 2-pools"
  ),
  data.frame(
    age = get_age_tt(TT3s1),
    density = TT3s1$transitTimeDensity,
    meanTT = TT3s1$meanTransitTime,
    model = "Control 3-pools"
  ),
  data.frame(
    age = get_age_tt(TT2s2),
    density = TT2s2$transitTimeDensity,
    meanTT = TT2s2$meanTransitTime,
    model = "Litter 2-pools"
  ),
  data.frame(
    age = get_age_tt(TT3s2),
    density = TT3s2$transitTimeDensity,
    meanTT = TT3s2$meanTransitTime,
    model = "Litter 3-pools"
  )
)

df_TT$model <- factor(
  cage$model,
  levels = c(
    "Control 2-pools",
    "Control 3-pools",
    "Litter 2-pools",
    "Litter 3-pools"
  )
)


model_labels_TT <- df_TT %>% 
  distinct(model, meanTT) %>%
  mutate(label = paste0(model, " = ", round(meanTT, 0), " yr")) %>%
  select(model, label) %>%     
  tibble::deframe()

TT_plot <- ggplot(df_TT, aes(x = age, y = density, color = model)) +
  geom_line(linewidth = 0.5) +
  geom_vline(
    aes(xintercept = meanTT, color = model),
    linetype = 5,
    linewidth = 0.5,
    show.legend = FALSE
  ) + 
  scale_color_manual(
    name = NULL,
    values = colors,
    labels = model_labels_TT,
    guide = guide_legend(override.aes = list(linewidth = 1))
  ) +
  coord_cartesian(xlim = c(0, 530), ylim = c(0, 0.001)) +
  labs(
    x = "Transit time (years)",
    y = "Density function"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = c(0.715, 0.75),   
    legend.background = element_rect(
      fill = NA,
      color = NA
    ),
    legend.key = element_blank()
  )

TT_plot

inc_CAge_TT = C_age + TT_plot
inc_CAge_TT

ggsave("./manuscript_figures/incubation_Cage_TT.png", inc_CAge_TT, width = 8.25, height = 3, dpi = 300)
ggsave("./manuscript_figures/incubation_Cage_TT.pdf", inc_CAge_TT, width = 8.5, height = 3, dpi = 300)
#####################################################################################

#####################################################################################

#####################################################################################

############# Model comparison: effects of structures and parameters changes #########################################
rm(list = ls())

############################### 2p1 # 2-pool baseline model
k2s1=c(0.00027099*365,0.00004431*365)
A2s1=diag(-k2s1)
A2s1[2,1]=k2s1[1]*0.02162174
u2s1=matrix(c(100,0),ncol=1)

############################### 2p2 # Model 1: 10% of C inputs enter into MAOC
k2s2=c(0.00027099*365,0.00004431*365)
A2s2=diag(-k2s2)
A2s2[2,1]=k2s2[1]*0.02162174
u2s2=matrix(c(90,10),ncol=1)

############################### 3p1 # 3-pool baseline model
k3s1=c(0.01878309*365, 0.00027099*365,0.00004431*365)
A3s1=diag(-k3s1)
A3s1[2,1]=k3s1[1]*0.10003170
A3s1[3,2]=k3s1[2]*0.02162174
u3s1=matrix(c(100,0,0),ncol=1)

############################### 3p2 # Model 2: direct C pathway from litter C to MAOC
k3s2=c(0.01878309*365, 0.00027099*365,0.00004431*365)
A3s2=diag(-k3s2)
A3s2[2,1]=k3s2[1]*0.10003170
A3s2[3,2]=k3s2[2]*0.02162174
A3s2[3,1]=k3s2[1]*0.03456093
u3s2=matrix(c(100,0,0),ncol=1)

############################### 3p3 # Model 3: triples k1
k3s3=c(0.01878309*365*3, 0.00027099*365,0.00004431*365)
A3s3=diag(-k3s3)
A3s3[2,1]=k3s3[1]*0.10003170
A3s3[3,2]=k3s3[2]*0.02162174
u3s3=matrix(c(100,0,0),ncol=1)

############################### 3p4 # Model 4: triples k1 and a21
k3s4=c(0.01878309*365*3, 0.00027099*365,0.00004431*365)
A3s4=diag(-k3s4)
A3s4[2,1]=k3s4[1]*0.10003170*3
A3s4[3,2]=k3s4[2]*0.02162174
u3s4=matrix(c(100,0,0),ncol=1)

############################### 3p5 # Model 5: triples k1, a21 and a32
k3s5=c(0.01878309*365*3, 0.00027099*365,0.00004431*365)
A3s5=diag(-k3s5)
A3s5[2,1]=k3s5[1]*0.10003170*3
A3s5[3,2]=k3s5[2]*0.02162174*3
u3s5=matrix(c(100,0,0),ncol=1)

###### 3p6 # Model 6: Adds a pathway from Litter C to MAOC and triples k1, a21 and a32
k3s6=c(0.01878309*365*3, 0.00027099*365,0.00004431*365)
A3s6=diag(-k3s6)
A3s6[2,1]=k3s6[1]*0.10003170*3
A3s6[3,2]=k3s6[2]*0.02162174*3
A3s6[3,1]=k3s6[1]*0.03456093
u3s6=matrix(c(100,0,0),ncol=1)

### 3p7 # Model 7: Adds a pathway from Litter C to MAOC and triples k1, a21, a32 and a31
k3s7=c(0.01878309*365*3, 0.00027099*365,0.00004431*365)
A3s7=diag(-k3s7)
A3s7[2,1]=k3s7[1]*0.10003170*3
A3s7[3,2]=k3s7[2]*0.02162174*3
A3s7[3,1]=k3s7[1]*0.03456093*3
u3s7=matrix(c(100,0,0),ncol=1)

##################################### Age and transit time
SA2s1=systemAge(a = seq(0, 2000),A=A2s1, u=u2s1)
SA2s2=systemAge(a = seq(0, 2000),A=A2s2, u=u2s2)
SA3s1=systemAge(a = seq(0, 2000),A=A3s1, u=u3s1)
SA3s2=systemAge(a = seq(0, 2000),A=A3s2, u=u3s2)
SA3s3=systemAge(a = seq(0, 2000),A=A3s3, u=u3s3)
SA3s4=systemAge(a = seq(0, 2000),A=A3s4, u=u3s4)
SA3s5=systemAge(a = seq(0, 2000),A=A3s5, u=u3s5)
SA3s6=systemAge(a = seq(0, 2000),A=A3s6, u=u3s6)
SA3s7=systemAge(a = seq(0, 2000),A=A3s7, u=u3s7)

TT2s1=transitTime(a = seq(0, 2000),A=A2s1, u=u2s1)
TT2s2=transitTime(a = seq(0, 2000),A=A2s2, u=u2s2)
TT3s1=transitTime(a = seq(0, 2000),A=A3s1, u=u3s1)
TT3s2=transitTime(a = seq(0, 2000),A=A3s2, u=u3s2)
TT3s3=transitTime(a = seq(0, 2000),A=A3s3, u=u3s3)
TT3s4=transitTime(a = seq(0, 2000),A=A3s4, u=u3s4)
TT3s5=transitTime(a = seq(0, 2000),A=A3s5, u=u3s5)
TT3s6=transitTime(a = seq(0, 2000),A=A3s6, u=u3s6)
TT3s7=transitTime(a = seq(0, 2000),A=A3s7, u=u3s7)

SA2s1$meanSystemAge
SA2s2$meanSystemAge
SA3s1$meanSystemAge
SA3s2$meanSystemAge
SA3s3$meanSystemAge
SA3s4$meanSystemAge
SA3s5$meanSystemAge
SA3s6$meanSystemAge
SA3s7$meanSystemAge

TT2s1$meanTransitTime
TT2s2$meanTransitTime
TT3s1$meanTransitTime
TT3s2$meanTransitTime
TT3s3$meanTransitTime
TT3s4$meanTransitTime
TT3s5$meanTransitTime
TT3s6$meanTransitTime
TT3s7$meanTransitTime

########################################## plot age

cage <- bind_rows(
  data.frame(age = 1:length(SA2s1$systemAgeDensity),
             density = SA2s1$systemAgeDensity,
             model = "2-Base",
             meanAge = SA2s1$meanSystemAge),
  data.frame(age = 1:length(SA2s2$systemAgeDensity),
             density = SA2s2$systemAgeDensity,
             model = "2.A",
             meanAge = SA2s2$meanSystemAge),
  data.frame(age = 1:length(SA3s1$systemAgeDensity),
             density = SA3s1$systemAgeDensity,
             model = "3-Base",
             meanAge = SA3s1$meanSystemAge),
  data.frame(age = 1:length(SA3s2$systemAgeDensity),
             density = SA3s2$systemAgeDensity,
             model = "3.A",
             meanAge = SA3s2$meanSystemAge),
  data.frame(age = 1:length(SA3s3$systemAgeDensity),
             density = SA3s3$systemAgeDensity,
             model = "3.B",
             meanAge = SA3s3$meanSystemAge),
  data.frame(age = 1:length(SA3s4$systemAgeDensity),
             density = SA3s4$systemAgeDensity,
             model = "3.C",
             meanAge = SA3s4$meanSystemAge),
  data.frame(age = 1:length(SA3s5$systemAgeDensity),
             density = SA3s5$systemAgeDensity,
             model = "3.D",
             meanAge = SA3s5$meanSystemAge),
  data.frame(age = 1:length(SA3s6$systemAgeDensity),
             density = SA3s6$systemAgeDensity,
             model = "3.E",
             meanAge = SA3s6$meanSystemAge),
  data.frame(age = 1:length(SA3s7$systemAgeDensity),
             density = SA3s7$systemAgeDensity,
             model = "3.F",
             meanAge = SA3s7$meanSystemAge)
)

cage$model <- factor(cage$model,
                     levels = c("2-Base", "2.A", "3-Base","3.A","3.B","3.C","3.D","3.E","3.F"))

model_labels <- cage %>% 
  distinct(model, meanAge) %>%
  mutate(label = paste0(model, " = ", round(meanAge, 1), " yr")) %>%
  select(model, label) %>%     
  tibble::deframe()

colors = c("#17154FFF", "#6C5D9EFF", "#B0799AFF", 
           "#F6B3B0FF", "#BF3729FF", "#E69B00FF", "#FFD380", 
           "#948619", "#355828FF")

C_age <- ggplot(cage, aes(x = age, y = density, color = model)) +
  geom_line(linewidth = 0.5) +
  geom_vline(aes(xintercept = meanAge, color = model),
             linetype = 5, linewidth = 1,
             show.legend = FALSE) + 
  scale_color_manual(name = NULL, values = colors,
                     labels = model_labels, 
                     guide = guide_legend(override.aes = list(linewidth = 1.5))) +
  coord_cartesian(xlim = c(15, 160), ylim = c(-0.0005, 0.01)) +
  labs(x = "C Age (years)", y = "Density function") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = c(0.99, 1.01),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(0.5, "cm")
  )

C_age 

ggsave("./manuscript_figures/C_age.png", C_age, width = 7, height = 3, dpi = 300)
ggsave("./manuscript_figures/C_age.pdf", C_age, width = 7, height = 3, dpi = 300)


# plot transit time 
tt <- bind_rows(
  data.frame(
    age     = seq(0, 2000),
    density = TT2s1$transitTimeDensity,
    model   = "2-Base",
    meanAge = TT2s1$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT2s2$transitTimeDensity,
    model   = "2.A",
    meanAge = TT2s2$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s1$transitTimeDensity,
    model   = "3-Base",
    meanAge = TT3s1$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s2$transitTimeDensity,
    model   = "3.A",
    meanAge = TT3s2$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s3$transitTimeDensity,
    model   = "3.B",
    meanAge = TT3s3$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s4$transitTimeDensity,
    model   = "3.C",
    meanAge = TT3s4$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s5$transitTimeDensity,
    model   = "3.D",
    meanAge = TT3s5$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s6$transitTimeDensity,
    model   = "3.E",
    meanAge = TT3s6$meanTransitTime
  ),
  
  data.frame(
    age     = seq(0, 2000),
    density = TT3s7$transitTimeDensity,
    model   = "3.F",
    meanAge = TT3s7$meanTransitTime
  )
)

tt$model <- factor(tt$model,
                   levels = c("2-Base", "2.A", "3-Base","3.A","3.B","3.C","3.D","3.E","3.F"))

tt_model_labels <- tt %>% 
  distinct(model, meanAge) %>%
  mutate(label = paste0(model, " = ", round(meanAge, 1), " yr")) %>%
  select(model, label) %>%     
  tibble::deframe()

colors = c("#17154FFF", "#6C5D9EFF", "#B0799AFF", 
           "#F6B3B0FF", "#BF3729FF", "#E69B00FF", "#FFD380", 
           "#948619", "#355828FF")

Tt = ggplot(tt, aes(x = age, y = density, color = model)) +
  geom_line(linewidth = 0.5) +
  geom_vline(aes(xintercept = meanAge, color = model),
             linetype = 5, linewidth = 1,
             show.legend = FALSE) + 
  scale_color_manual(
    name = NULL,
    values = colors,
    labels = tt_model_labels, 
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  coord_cartesian(xlim = c(0, 75), ylim = c(-0.0005, 0.0105)) +
  labs(x = "Transit time (years)", y = "Density function") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = c(0.99, 1.01),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.height = unit(0.5, "cm")
  )
Tt

ggsave("./manuscript_figures/transittime.png", Tt, width = 7, height = 3, dpi = 300)
ggsave("./manuscript_figures/transittime.pdf", Tt, width = 7, height = 3, dpi = 300)
#


library(patchwork)

final_plot <- C_age / Tt
final_plot

ggsave("./manuscript_figures/final_plot.png", final_plot, width = 7, height = 6, dpi = 300)
ggsave("./manuscript_figures/final_plot.pdf", final_plot, width = 7, height = 6, dpi = 300)

