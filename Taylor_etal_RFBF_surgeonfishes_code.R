#### Code for Taylor et al - Comparative demography of surgeonfishes from the tropical Western Pacific ############
### 2023 ###

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(investr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(nlstools)
library(scales)
library(ellipse)
library(car)

NUdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Naso unicornis", col_types = c("text", 
                      "text", "text", "date", "text", "numeric", 
                      "numeric", "numeric", "text", 
                      "numeric", "numeric", "numeric"))
ANdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Acanthurus nigricauda", col_types = c("text", 
                     "text", "text", "date", "text", "numeric", 
                      "numeric", "numeric", "text", "numeric", 
                      "numeric", "numeric"))
ATdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Acanthurus triostegus", col_types = c("text", 
                     "text", "text", "date", "text", "numeric", 
                     "numeric", "numeric", "text", "numeric", 
                     "numeric", "numeric"))
ABdat <- read_excel("RFBF_Species_Data.xlsx", 
                     sheet = "Acanthurus blochii", col_types = c("text", 
                                                                 "text", "text", "date", "text", "numeric", 
                                                                 "numeric", "text", "numeric", 
                                                                 "numeric", "numeric"))
AGdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Acanthurus guttatus", col_types = c("text", 
                                                                "text", "text", "date", "text", "numeric", 
                                                                "numeric", "text", "numeric", 
                                                                "numeric", "numeric"))
ALdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Acanthurus lineatus", col_types = c("text", 
                                                                "text", "text", "date", "text", "numeric", 
                                                                "numeric", "text", "numeric", 
                                                                "numeric", "numeric"))
AOdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Acanthurus olivaceus", col_types = c("text", 
                                                                 "text", "text", "date", "text", "numeric", 
                                                                 "numeric", "text", "numeric", 
                                                                 "numeric", "numeric"))
AXdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Acanthurus xanthopterus", col_types = c("text", 
                                                                  "text", "text", "date", "text", "numeric", 
                                                                  "numeric", "text", "numeric", 
                                                                  "numeric", "numeric"))
CStdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Ctenochaetus striatus", col_types = c("text", 
                                                                     "text", "text", "date", "text", "numeric", 
                                                                     "numeric", "text", "numeric", 
                                                                     "numeric", "numeric"))
NTdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Naso tonganus", col_types = c("text", 
                                                                     "text", "text", "date", "text", "numeric", 
                                                                     "numeric", "text", "numeric", 
                                                                     "numeric", "numeric"))
NVdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Naso vlamingii", col_types = c("text", 
                                                                     "text", "text", "date", "text", "numeric", 
                                                                     "numeric", "text", "numeric", 
                                                                     "numeric", "numeric"))
ZVdat <- read_excel("RFBF_Species_Data.xlsx", 
                    sheet = "Zebrasoma velifer", col_types = c("text", 
                                                                     "text", "text", "date", "text", "numeric", 
                                                                     "numeric", "text", "numeric", 
                                                                     "numeric", "numeric"))


##### Female maturity ###########
#################################
matfunc <- maturity ~ (1 + exp(-log(19)*(length - L50)/(L95 - L50)))^-1
agematfunc <- maturity ~ (1 + exp(-log(19)*(age - t50)/(t95 - t50)))^-1

#Subset the data in new data frames to deal with NAs
ABmatdf <- subset(ABdat, finalsex=="female")
ABmatdf <- ABmatdf[,c(6,9)]
ABmatdf <- na.omit(ABmatdf)

AGmatdf <- subset(AGdat, finalsex=="female")
AGmatdf <- AGmatdf[,c(6,9)]
AGmatdf <- na.omit(AGmatdf)

ALmatdf <- subset(ALdat, finalsex=="female")
ALmatdf <- ALmatdf[,c(6,9)]
ALmatdf <- na.omit(ALmatdf)

AOmatdf <- subset(AOdat, finalsex=="female")
AOmatdf <- AOmatdf[,c(6,9)]
AOmatdf <- na.omit(AOmatdf)

AXmatdf <- subset(AXdat, finalsex=="female")
AXmatdf <- AXmatdf[,c(6,9)]
AXmatdf <- na.omit(AXmatdf)

CStmatdf <- subset(CStdat, finalsex=="female")
CStmatdf <- CStmatdf[,c(6,9)]
CStmatdf <- na.omit(CStmatdf)

NTmatdf <- subset(NTdat, finalsex=="female")
NTmatdf <- NTmatdf[,c(6,9)]
NTmatdf <- na.omit(NTmatdf)

NVmatdf <- subset(NVdat, finalsex=="female")
NVmatdf <- NVmatdf[,c(6,9)]
NVmatdf <- na.omit(NVmatdf)

ZVmatdf <- subset(ZVdat, finalsex=="female")
ZVmatdf <- ZVmatdf[,c(6,9)]
ZVmatdf <- na.omit(ZVmatdf)

NUmatdf <- subset(NUdat, finalsex=="female")
NUmatdf <- NUmatdf[,c(6,10)]
NUmatdf <- na.omit(NUmatdf)

ANmatdf <- subset(ANdat, finalsex=="female")
ANmatdf <- ANmatdf[,c(6,10)]
ANmatdf <- na.omit(ANmatdf)

ATmatdf <- subset(ATdat, finalsex=="female")
ATmatdf <- ATmatdf[,c(6,10)]
ATmatdf <- na.omit(ATmatdf)




###Model fitting for logistic models
ABmat <- nls(matfunc, data=ABmatdf, start = c(L50 = 250, L95 = 300), algorithm = "port")
AGmat <- nls(matfunc, data=AGmatdf, start = c(L50 = 150, L95 = 200), algorithm = "port")
ALmat <- nls(matfunc, data=ALmatdf, start = c(L50 = 180, L95 = 200), algorithm = "port")
AOmat <- nls(matfunc, data=AOmatdf, start = c(L50 = 180, L95 = 250), algorithm = "port")
AXmat <- nls(matfunc, data=AXmatdf, start = c(L50 = 250, L95 = 300), algorithm = "port")
CStmat <- nls(matfunc, data=CStmatdf, start = c(L50 = 250, L95 = 300), algorithm = "port")
NTmat <- nls(matfunc, data=NTmatdf, start = c(L50 = 250, L95 = 300), algorithm = "port")
NVmat <- nls(matfunc, NVmatdf, start = c(L50 = 150, L95 = 200), algorithm = "port") 
NUmat <- nls(matfunc, data=NUmatdf, start = c(L50 = 300, L95 = 400), algorithm = "port")
ANmat <- nls(matfunc, data=ANmatdf, start = c(L50 = 200, L95 = 250), algorithm = "port")
ATmat <- nls(matfunc, ATmatdf, start = c(L50 = 150, L95 = 170), algorithm = "port") 
ZVmat <- nls(matfunc, ZVmatdf, start = c(L50 = 150, L95 = 200), algorithm = "port") 



##bootstrapping estimates

new.data <- data.frame(length=seq(50, 600, by = 1))


NUboot <- nlsBoot(NUmat, niter = 999)
NUint <- as_tibble(nlsBootPredict(NUboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(NUint) <- c("fit", "lower", "upper", "length")

ANboot <- nlsBoot(ANmat, niter = 999)
ANint <- as_tibble(nlsBootPredict(ANboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(ANint) <- c("fit", "lower", "upper", "length")

ATboot <- nlsBoot(ATmat, niter = 999)
ATint <- as_tibble(nlsBootPredict(ATboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(ATint) <- c("fit", "lower", "upper", "length")

ABboot <- nlsBoot(ABmat, niter = 999)
ABint <- as_tibble(nlsBootPredict(ABboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(ABint) <- c("fit", "lower", "upper", "length")

AGboot <- nlsBoot(AGmat, niter = 999)
AGint <- as_tibble(nlsBootPredict(AGboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(AGint) <- c("fit", "lower", "upper", "length")

ALboot <- nlsBoot(ALmat, niter = 999)
ALint <- as_tibble(nlsBootPredict(ALboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(ALint) <- c("fit", "lower", "upper", "length")

AOboot <- nlsBoot(AOmat, niter = 999)
AOint <- as_tibble(nlsBootPredict(AOboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(AOint) <- c("fit", "lower", "upper", "length")

AXboot <- nlsBoot(AXmat, niter = 999)
AXint <- as_tibble(nlsBootPredict(AXboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(AXint) <- c("fit", "lower", "upper", "length")

CStboot <- nlsBoot(CStmat, niter = 999)
CStint <- as_tibble(nlsBootPredict(CStboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(CStint) <- c("fit", "lower", "upper", "length")

NTboot <- nlsBoot(NTmat, niter = 999)
NTint <- as_tibble(nlsBootPredict(NTboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(NTint) <- c("fit", "lower", "upper", "length")

NVboot <- nlsBoot(NVmat, niter = 999)
NVint <- as_tibble(nlsBootPredict(NVboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(NVint) <- c("fit", "lower", "upper", "length")

ZVboot <- nlsBoot(ZVmat, niter = 999)
ZVint <- as_tibble(nlsBootPredict(ZVboot, newdata = new.data, interval = "confidence")) %>% 
  mutate(length = new.data$length)
colnames(ZVint) <- c("fit", "lower", "upper", "length")



#plotting maturity

pNTmat <- ggplot(NTmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=NTint, aes(x = length, y = fit ))+ xlim(200, 600) +
  geom_ribbon(data=NTint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Naso tonganus") + theme_light()+
  geom_segment(aes(x = 200, xend = environment(NTmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(NTmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(NTmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pNUmat <- ggplot(NUmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=NUint, aes(x = length, y = fit))+ xlim(50, 600) +
  geom_ribbon(data=NUint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Naso unicornis")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(NUmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(NUmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(NUmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pANmat <- ggplot(ANmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=ANint, aes(x = length, y = fit))+ xlim(120, 250) +
  geom_ribbon(data=ANint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus nigricauda")+ theme_light()+
  geom_segment(aes(x = 120, xend = environment(ANmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(ANmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(ANmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pATmat <- ggplot(ATmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=ATint, aes(x = length, y = fit))+ xlim(50, 200) +
  geom_ribbon(data=ATint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus triostegus")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(ATmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(ATmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(ATmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pABmat <- ggplot(ABmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=ABint, aes(x = length, y = fit))+ xlim(100, 350) +
  geom_ribbon(data=ABint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus blochii")+ theme_light()+
  geom_segment(aes(x = 100, xend = environment(ABmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(ABmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(ABmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pAGmat <- ggplot(AGmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=AGint, aes(x = length, y = fit))+ xlim(50, 250) +
  geom_ribbon(data=AGint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus guttatus")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(AGmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(AGmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(AGmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pALmat <- ggplot(ALmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=ALint, aes(x = length, y = fit))+ xlim(50, 250) +
  geom_ribbon(data=ALint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus lineatus")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(ALmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(ALmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(ALmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pAOmat <- ggplot(AOmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=AOint, aes(x = length, y = fit))+ xlim(50, 250) +
  geom_ribbon(data=AOint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus olivaceus")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(AOmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(AOmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(AOmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pAXmat <- ggplot(AXmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=AXint, aes(x = length, y = fit))+ xlim(100, 450) +
  geom_ribbon(data=AXint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Acanthurus xanthopterus")+ theme_light()+
  geom_segment(aes(x = 100, xend = environment(AXmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(AXmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(AXmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pCStmat <- ggplot(CStmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=CStint, aes(x = length, y = fit))+ xlim(50, 250) +
  geom_ribbon(data=CStint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Ctenochaetus striatus")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(CStmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(CStmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(CStmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pNVmat <- ggplot(NVmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=NVint, aes(x = length, y = fit))+ xlim(100, 400) +
  geom_ribbon(data=NVint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Naso vlamingii")+ theme_light()+
  geom_segment(aes(x = 100, xend = environment(NVmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(NVmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(NVmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")

pZVmat <- ggplot(ZVmatdf) +  geom_point(aes(x=length, y=maturity),size=2, colour="black", alpha = 0.25) + xlab("Fork length (mm)") + ylab("Proportion mature")  +
  geom_line(data=ZVint, aes(x = length, y = fit))+ xlim(50, 250) +
  geom_ribbon(data=ZVint, aes(x=length, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey")+ ggtitle("Zebrasoma velifer")+ theme_light()+
  geom_segment(aes(x = 50, xend = environment(ZVmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0.5, yend = .5), color ="red", size =.5, linetype = "dashed")+
  geom_segment(aes(x = environment(ZVmat[["m"]][["formula"]])[["internalPars"]][1], xend = environment(ZVmat[["m"]][["formula"]])[["internalPars"]][1] , y = 0, yend = .5), color ="red", size =.5, linetype = "dashed")



grid.arrange(pABmat, pAGmat, pALmat, pANmat, pAOmat, pATmat, pAXmat, pCStmat, pNTmat, pNUmat, pNVmat, pZVmat, nrow = 4)




#################################################
#################GROWTH##########################
#################################################

#Subset data

ABfgrodat <- subset(ABdat, finalsex=="female")
ABfgrodat <- ABfgrodat[,c(6,11)]
ABfgrodat <- na.omit(ABfgrodat)

ABmgrodat <- subset(ABdat, finalsex=="male")
ABmgrodat <- ABmgrodat[,c(6,11)]
ABmgrodat <- na.omit(ABmgrodat)

AGfgrodat <- subset(AGdat, finalsex=="female")
AGfgrodat <- AGfgrodat[,c(6,11)]
AGfgrodat <- na.omit(AGfgrodat)

AGmgrodat <- subset(AGdat, finalsex=="male")
AGmgrodat <- AGmgrodat[,c(6,11)]
AGmgrodat <- na.omit(AGmgrodat)

ALfgrodat <- subset(ALdat, finalsex=="female")
ALfgrodat <- ALfgrodat[,c(6,11)]
ALfgrodat <- na.omit(ALfgrodat)

ALmgrodat <- subset(ALdat, finalsex=="male")
ALmgrodat <- ALmgrodat[,c(6,11)]
ALmgrodat <- na.omit(ALmgrodat)

AOfgrodat <- subset(AOdat, finalsex=="female")
AOfgrodat <- AOfgrodat[,c(6,11)]
AOfgrodat <- na.omit(AOfgrodat)

AOmgrodat <- subset(AOdat, finalsex=="male")
AOmgrodat <- AOmgrodat[,c(6,11)]
AOmgrodat <- na.omit(AOmgrodat)

AXfgrodat <- subset(AXdat, finalsex=="female")
AXfgrodat <- AXfgrodat[,c(6,11)]
AXfgrodat <- na.omit(AXfgrodat)

AXmgrodat <- subset(AXdat, finalsex=="male")
AXmgrodat <- AXmgrodat[,c(6,11)]
AXmgrodat <- na.omit(AXmgrodat)

CStfgrodat <- subset(CStdat, finalsex=="female")
CStfgrodat <- CStfgrodat[,c(6,11)]
CStfgrodat <- na.omit(CStfgrodat)

CStmgrodat <- subset(CStdat, finalsex=="male")
CStmgrodat <- CStmgrodat[,c(6,11)]
CStmgrodat <- na.omit(CStmgrodat)

NTfgrodat <- subset(NTdat, finalsex=="female")
NTfgrodat <- NTfgrodat[,c(6,11)]
NTfgrodat <- na.omit(NTfgrodat)

NTmgrodat <- subset(NTdat, finalsex=="male")
NTmgrodat <- NTmgrodat[,c(6,11)]
NTmgrodat <- na.omit(NTmgrodat)

NUfgrodat <- subset(NUdat, finalsex=="female")
NUfgrodat <- NUfgrodat[,c(6,12)]
NUfgrodat <- na.omit(NUfgrodat)

NUmgrodat <- subset(NUdat, finalsex=="male")
NUmgrodat <- NUmgrodat[,c(6,12)]
NUmgrodat <- na.omit(NUmgrodat)

NVfgrodat <- subset(NVdat, finalsex=="female")
NVfgrodat <- NVfgrodat[,c(6,11)]
NVfgrodat <- na.omit(NVfgrodat)

NVmgrodat <- subset(NVdat, finalsex=="male")
NVmgrodat <- NVmgrodat[,c(6,11)]
NVmgrodat <- na.omit(NVmgrodat)

ANfgrodat <- subset(ANdat, finalsex=="female")
ANfgrodat <- ANfgrodat[,c(6,12)]
ANfgrodat <- na.omit(ANfgrodat)

ANmgrodat <- subset(ANdat, finalsex=="male")
ANmgrodat <- ANmgrodat[,c(6,12)]
ANmgrodat <- na.omit(ANmgrodat)

ATfgrodat <- subset(ATdat, finalsex=="female")
ATfgrodat <- ATfgrodat[,c(6,12)]
ATfgrodat <- na.omit(ATfgrodat)

ATmgrodat <- subset(ATdat, finalsex=="male")
ATmgrodat <- ATmgrodat[,c(6,12)]
ATmgrodat <- na.omit(ATmgrodat)

ZVfgrodat <- subset(ZVdat, finalsex=="female")
ZVfgrodat <- ZVfgrodat[,c(6,11)]
ZVfgrodat <- na.omit(ZVfgrodat)

ZVmgrodat <- subset(ZVdat, finalsex=="male")
ZVmgrodat <- ZVmgrodat[,c(6,11)]
ZVmgrodat <- na.omit(ZVmgrodat)


AB <- ABdat[,c(6,8,11)]
AB <- na.omit(AB)
AG <- AGdat[,c(6,8,11)]
AG <- na.omit(AG)
AL <- ALdat[,c(6,8,11)]
AL <- na.omit(AL)
AO <- AOdat[,c(6,8,11)]
AO <- na.omit(AO)
AX <- AXdat[,c(6,8,11)]
AX <- na.omit(AX)
CSt <- CStdat[,c(6,8,11)]
CSt <- na.omit(CSt)
NT <- NTdat[,c(6,8,11)]
NT <- na.omit(NT)
NV <- NVdat[,c(6,8,11)]
NV <- na.omit(NV)
ZV <- ZVdat[,c(6,8,11)]
ZV <- na.omit(ZV)
NU <- NUdat[,c(6,9,12)]
NU <- na.omit(NU)
NU <- NU[-c(600:610),] #Removing NAs from character string
AN <- ANdat[,c(6,9,12)]
AN <- na.omit(AN)
AN <- AN[-c(455:559),] #Removing NAs from character string
AT <- ATdat[,c(6,9,12)]
AT <- na.omit(AT)


####Fitting growth models using non-linear least squares
NUfGro <- nls(length ~ Linf-(Linf-50)*exp(-K*age), data=NUfgrodat, start = c(Linf = 400, K = 0.5), algorithm = "port")
NUmGro <- nls(length ~ Linf-(Linf-50)*exp(-K*age), data=NUmgrodat, start = c(Linf = 400, K = 0.5), algorithm = "port")
NUGro <- nls(length ~Linf-(Linf-50)*exp(-K*age), data = NU, start = c(Linf = 400, K = 0.5), algorithm = "port")
NUfGro
NUmGro
NUGro


#Naso unicornis growth
new.grodatNU <- data.frame(age=seq(0, 22, by = .25))

NUfGroboot <- nlsBoot(NUfGro, niter = 999)
NUfGroint <- as_tibble(nlsBootPredict(NUfGroboot, newdata = new.grodatNU, interval = "confidence")) %>% 
  mutate(age = new.grodatNU$age)
colnames(NUfGroint) <- c("fit", "lower", "upper", "age")

NUmGroboot <- nlsBoot(NUmGro, niter = 999)
NUmGroint <- as_tibble(nlsBootPredict(NUmGroboot, newdata = new.grodatNU, interval = "confidence")) %>% 
  mutate(age = new.grodatNU$age)
colnames(NUmGroint) <- c("fit", "lower", "upper", "age")

NUGroboot <- nlsBoot(NUGro, niter = 999)
NUGroint <- as_tibble(nlsBootPredict(NUGroboot, newdata = new.grodatNU, interval = "confidence")) %>% 
  mutate(age = new.grodatNU$age)
colnames(NUGroint) <- c("fit", "lower", "upper", "age")

estsNUm <- NUmGroboot$coefboot
estsNUf <- NUfGroboot$coefboot

plot(Linf~K, data = estsNUm, main = "Males")
plot(Linf~K, data = estsNUf, main = "Females")
estsNUfdf <- as.data.frame(estsNUf)
estsNUmdf <- as.data.frame(estsNUm)
NUfplot <- ggplot(estsNUfdf, aes(K, Linf)) + geom_point()
NUbiv <- NUfplot + geom_point(data = estsNUmdf, aes(K, Linf), color = "grey") + theme_bw()
NUbiv

vbNU <- ggplot(NU) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=NUmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=NUfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=NUGroint, aes(x = age, y = fit ))+ 
  xlim(0, 25) + ylim(0, 600)+
  geom_ribbon(data=NUGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Naso unicornis") + theme_classic()
vbNU


#Acanthurus nigricauda growth
ANfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ANfgrodat, start = c(Linf = 200, K = 0.5), algorithm = "port")
ANmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ANmgrodat, start = c(Linf = 200, K = 0.5), algorithm = "port")
ANGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AN, start = c(Linf = 200, K = 0.5), algorithm = "port")
ANfGro
ANmGro
ANGro


new.grodatAN <- data.frame(age=seq(0, 17, by = .25))

ANfGroboot <- nlsBoot(ANfGro, niter = 999)
ANfGroint <- as_tibble(nlsBootPredict(ANfGroboot, newdata = new.grodatAN, interval = "confidence")) %>% 
  mutate(age = new.grodatAN$age)
colnames(ANfGroint) <- c("fit", "lower", "upper", "age")

ANmGroboot <- nlsBoot(ANmGro, niter = 999)
ANmGroint <- as_tibble(nlsBootPredict(ANmGroboot, newdata = new.grodatAN, interval = "confidence")) %>% 
  mutate(age = new.grodatAN$age)
colnames(ANmGroint) <- c("fit", "lower", "upper", "age")

ANGroboot <- nlsBoot(ANGro, niter = 999)
ANGroint <- as_tibble(nlsBootPredict(ANGroboot, newdata = new.grodatAN, interval = "confidence")) %>% 
  mutate(age = new.grodatAN$age)
colnames(ANGroint) <- c("fit", "lower", "upper", "age")

estsANm <- ANmGroboot$coefboot
estsANf <- ANfGroboot$coefboot

plot(Linf~K, data = estsANm, main = "Males")
plot(Linf~K, data = estsANf, main = "Females")
estsANfdf <- as.data.frame(estsANf)
estsANmdf <- as.data.frame(estsANm)
ANfplot <- ggplot(estsANfdf, aes(K, Linf)) + geom_point()
ANbiv <- ANfplot + geom_point(data = estsANmdf, aes(K, Linf), color = "grey") + theme_bw()
ANbiv

vbAN <- ggplot(AN) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=ANmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=ANfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=ANGroint, aes(x = age, y = fit ))+ 
  xlim(0, 20) + ylim(0, 250)+
  geom_ribbon(data=ANGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus nigricauda") + theme_classic()
vbAN


#Acanthurus triostegus growth
ATfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ATfgrodat, start = c(Linf = 150, K = 0.5), algorithm = "port")
ATmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ATmgrodat, start = c(Linf = 150, K = 0.5), algorithm = "port")
ATGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AT, start = c(Linf = 150, K = 0.5), algorithm = "port")
ATfGro
ATmGro
ATGro


new.grodatAT <- data.frame(age=seq(0, 12, by = .25))

ATfGroboot <- nlsBoot(ATfGro, niter = 999)
ATfGroint <- as_tibble(nlsBootPredict(ATfGroboot, newdata = new.grodatAT, interval = "confidence")) %>% 
  mutate(age = new.grodatAT$age)
colnames(ATfGroint) <- c("fit", "lower", "upper", "age")

ATmGroboot <- nlsBoot(ATmGro, niter = 999)
ATmGroint <- as_tibble(nlsBootPredict(ATmGroboot, newdata = new.grodatAT, interval = "confidence")) %>% 
  mutate(age = new.grodatAT$age)
colnames(ATmGroint) <- c("fit", "lower", "upper", "age")

ATGroboot <- nlsBoot(ATGro, niter = 999)
ATGroint <- as_tibble(nlsBootPredict(ATGroboot, newdata = new.grodatAT, interval = "confidence")) %>% 
  mutate(age = new.grodatAT$age)
colnames(ATGroint) <- c("fit", "lower", "upper", "age")

estsATm <- ATmGroboot$coefboot
estsATf <- ATfGroboot$coefboot

plot(Linf~K, data = estsATm, main = "Males")
plot(Linf~K, data = estsATf, main = "Females")
estsATfdf <- as.data.frame(estsATf)
estsATmdf <- as.data.frame(estsATm)
ATfplot <- ggplot(estsATfdf, aes(K, Linf)) + geom_point()
ATbiv <- ATfplot + geom_point(data = estsATmdf, aes(K, Linf), color = "grey") + theme_bw()
ATbiv

vbAT <- ggplot(AT) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=ATmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=ATfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=ATGroint, aes(x = age, y = fit ))+ 
  xlim(0, 15) + ylim(0, 200)+
  geom_ribbon(data=ATGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus triostegus") + theme_classic()
vbAT


#Naso tonganus growth
NTfGro <- nls(length ~ Linf-(Linf-50)*exp(-K*age), data=NTfgrodat, start = c(Linf = 400, K = 0.5), algorithm = "port")
NTmGro <- nls(length ~ Linf-(Linf-50)*exp(-K*age), data=NTmgrodat, start = c(Linf = 400, K = 0.5), algorithm = "port")
NTGro <- nls(length ~Linf-(Linf-50)*exp(-K*age), data = NT, start = c(Linf = 400, K = 0.5), algorithm = "port")
NTfGro
NTmGro
NTGro


new.grodatNT <- data.frame(age=seq(0, 23, by = .25))

NTfGroboot <- nlsBoot(NTfGro, niter = 999)
NTfGroint <- as_tibble(nlsBootPredict(NTfGroboot, newdata = new.grodatNT, interval = "confidence")) %>% 
  mutate(age = new.grodatNT$age)
colnames(NTfGroint) <- c("fit", "lower", "upper", "age")

NTmGroboot <- nlsBoot(NTmGro, niter = 999)
NTmGroint <- as_tibble(nlsBootPredict(NTmGroboot, newdata = new.grodatNT, interval = "confidence")) %>% 
  mutate(age = new.grodatNT$age)
colnames(NTmGroint) <- c("fit", "lower", "upper", "age")

NTGroboot <- nlsBoot(NTGro, niter = 999)
NTGroint <- as_tibble(nlsBootPredict(NTGroboot, newdata = new.grodatNT, interval = "confidence")) %>% 
  mutate(age = new.grodatNT$age)
colnames(NTGroint) <- c("fit", "lower", "upper", "age")

estsNTm <- NTmGroboot$coefboot
estsNTf <- NTfGroboot$coefboot

plot(Linf~K, data = estsNTm, main = "Males")
plot(Linf~K, data = estsNTf, main = "Females")
estsNTfdf <- as.data.frame(estsNTf)
estsNTmdf <- as.data.frame(estsNTm)
NTfplot <- ggplot(estsNTfdf, aes(K, Linf)) + geom_point()
NTbiv <- NTfplot + geom_point(data = estsNTmdf, aes(K, Linf), color = "grey") + theme_bw()
NTbiv

vbNT <- ggplot(NT) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=NTmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=NTfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=NTGroint, aes(x = age, y = fit ))+ 
  xlim(0, 21) + ylim(0, 620)+
  geom_ribbon(data=NTGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Naso tonganus") + theme_classic()
vbNT



#Naso vlamingii growth
NVfGro <- nls(length ~ Linf-(Linf-50)*exp(-K*age), data=NVfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
NVmGro <- nls(length ~ Linf-(Linf-50)*exp(-K*age), data=NVmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
NVGro <- nls(length ~Linf-(Linf-50)*exp(-K*age), data = NV, start = c(Linf = 300, K = 0.5), algorithm = "port")
NVfGro
NVmGro
NVGro



new.grodatNV <- data.frame(age=seq(0, 9, by = .25))

NVfGroboot <- nlsBoot(NVfGro, niter = 999)
NVfGroint <- as_tibble(nlsBootPredict(NVfGroboot, newdata = new.grodatNV, interval = "confidence")) %>% 
  mutate(age = new.grodatNV$age)
colnames(NVfGroint) <- c("fit", "lower", "upper", "age")

NVmGroboot <- nlsBoot(NVmGro, niter = 999)
NVmGroint <- as_tibble(nlsBootPredict(NVmGroboot, newdata = new.grodatNV, interval = "confidence")) %>% 
  mutate(age = new.grodatNV$age)
colnames(NVmGroint) <- c("fit", "lower", "upper", "age")

NVGroboot <- nlsBoot(NVGro, niter = 999)
NVGroint <- as_tibble(nlsBootPredict(NVGroboot, newdata = new.grodatNV, interval = "confidence")) %>% 
  mutate(age = new.grodatNV$age)
colnames(NVGroint) <- c("fit", "lower", "upper", "age")

estsNVm <- NVmGroboot$coefboot
estsNVf <- NVfGroboot$coefboot

plot(Linf~K, data = estsNVm, main = "Males")
plot(Linf~K, data = estsNVf, main = "Females")
estsNVfdf <- as.data.frame(estsNVf)
estsNVmdf <- as.data.frame(estsNVm)
NVfplot <- ggplot(estsNVfdf, aes(K, Linf)) + geom_point()
NVbiv <- NVfplot + geom_point(data = estsNVmdf, aes(K, Linf), color = "grey") + theme_bw()
NVbiv

vbNV <- ggplot(NV) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=NVmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=NVfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=NVGroint, aes(x = age, y = fit ))+ 
  xlim(0, 10) + ylim(0, 400)+
  geom_ribbon(data=NVGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Naso vlamingii") + theme_classic()
vbNV



#Zebrasoma velifer growth
ZVfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ZVfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
ZVmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ZVmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
ZVGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = ZV, start = c(Linf = 300, K = 0.5), algorithm = "port")
ZVfGro
ZVmGro
ZVGro



new.grodatZV <- data.frame(age=seq(0, 22, by = .25))

ZVfGroboot <- nlsBoot(ZVfGro, niter = 999)
ZVfGroint <- as_tibble(nlsBootPredict(ZVfGroboot, newdata = new.grodatZV, interval = "confidence")) %>% 
  mutate(age = new.grodatZV$age)
colnames(ZVfGroint) <- c("fit", "lower", "upper", "age")

ZVmGroboot <- nlsBoot(ZVmGro, niter = 999)
ZVmGroint <- as_tibble(nlsBootPredict(ZVmGroboot, newdata = new.grodatZV, interval = "confidence")) %>% 
  mutate(age = new.grodatZV$age)
colnames(ZVmGroint) <- c("fit", "lower", "upper", "age")

ZVGroboot <- nlsBoot(ZVGro, niter = 999)
ZVGroint <- as_tibble(nlsBootPredict(ZVGroboot, newdata = new.grodatZV, interval = "confidence")) %>% 
  mutate(age = new.grodatZV$age)
colnames(ZVGroint) <- c("fit", "lower", "upper", "age")

estsZVm <- ZVmGroboot$coefboot
estsZVf <- ZVfGroboot$coefboot

plot(Linf~K, data = estsZVm, main = "Males")
plot(Linf~K, data = estsZVf, main = "Females")
estsZVfdf <- as.data.frame(estsZVf)
estsZVmdf <- as.data.frame(estsZVm)
ZVfplot <- ggplot(estsZVfdf, aes(K, Linf)) + geom_point()
ZVbiv <- ZVfplot + geom_point(data = estsZVmdf, aes(K, Linf), color = "grey") + theme_bw()
ZVbiv

vbZV <- ggplot(ZV) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=ZVmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=ZVfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=ZVGroint, aes(x = age, y = fit ))+ 
  xlim(0, 22) + ylim(0, 260)+
  geom_ribbon(data=ZVGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Zebrasoma velifer") + theme_classic()
vbZV



#Ctenochaetus striatus growth
CStfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=CStfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
CStmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=CStmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
CStGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = CSt, start = c(Linf = 300, K = 0.5), algorithm = "port")
CStfGro
CStmGro
CStGro



new.grodatCSt <- data.frame(age=seq(0, 20, by = .25))

CStfGroboot <- nlsBoot(CStfGro, niter = 999)
CStfGroint <- as_tibble(nlsBootPredict(CStfGroboot, newdata = new.grodatCSt, interval = "confidence")) %>% 
  mutate(age = new.grodatCSt$age)
colnames(CStfGroint) <- c("fit", "lower", "upper", "age")

CStmGroboot <- nlsBoot(CStmGro, niter = 999)
CStmGroint <- as_tibble(nlsBootPredict(CStmGroboot, newdata = new.grodatCSt, interval = "confidence")) %>% 
  mutate(age = new.grodatCSt$age)
colnames(CStmGroint) <- c("fit", "lower", "upper", "age")

CStGroboot <- nlsBoot(CStGro, niter = 999)
CStGroint <- as_tibble(nlsBootPredict(CStGroboot, newdata = new.grodatCSt, interval = "confidence")) %>% 
  mutate(age = new.grodatCSt$age)
colnames(CStGroint) <- c("fit", "lower", "upper", "age")

estsCStm <- CStmGroboot$coefboot
estsCStf <- CStfGroboot$coefboot

plot(Linf~K, data = estsCStm, main = "Males")
plot(Linf~K, data = estsCStf, main = "Females")
estsCStfdf <- as.data.frame(estsCStf)
estsCStmdf <- as.data.frame(estsCStm)
CStfplot <- ggplot(estsCStfdf, aes(K, Linf)) + geom_point()
CStbiv <- CStfplot + geom_point(data = estsCStmdf, aes(K, Linf), color = "grey") + theme_bw()
CStbiv

vbCSt <- ggplot(CSt) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=CStmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=CStfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=CStGroint, aes(x = age, y = fit ))+ 
  xlim(0, 17) + ylim(0, 250)+
  geom_ribbon(data=CStGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Ctenochaetus striatus") + theme_classic()
vbCSt



#Acanthurus blochii growth
ABfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ABfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
ABmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ABmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
ABGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AB, start = c(Linf = 300, K = 0.5), algorithm = "port")
ABfGro
ABmGro
ABGro



new.grodatAB <- data.frame(age=seq(0, 20, by = .25))

ABfGroboot <- nlsBoot(ABfGro, niter = 999)
ABfGroint <- as_tibble(nlsBootPredict(ABfGroboot, newdata = new.grodatAB, interval = "confidence")) %>% 
  mutate(age = new.grodatAB$age)
colnames(ABfGroint) <- c("fit", "lower", "upper", "age")

ABmGroboot <- nlsBoot(ABmGro, niter = 999)
ABmGroint <- as_tibble(nlsBootPredict(ABmGroboot, newdata = new.grodatAB, interval = "confidence")) %>% 
  mutate(age = new.grodatAB$age)
colnames(ABmGroint) <- c("fit", "lower", "upper", "age")

ABGroboot <- nlsBoot(ABGro, niter = 999)
ABGroint <- as_tibble(nlsBootPredict(ABGroboot, newdata = new.grodatAB, interval = "confidence")) %>% 
  mutate(age = new.grodatAB$age)
colnames(ABGroint) <- c("fit", "lower", "upper", "age")

estsABm <- ABmGroboot$coefboot
estsABf <- ABfGroboot$coefboot

plot(Linf~K, data = estsABm, main = "Males")
plot(Linf~K, data = estsABf, main = "Females")
estsABfdf <- as.data.frame(estsABf)
estsABmdf <- as.data.frame(estsABm)
ABfplot <- ggplot(estsABfdf, aes(K, Linf)) + geom_point()
ABbiv <- ABfplot + geom_point(data = estsABmdf, aes(K, Linf), color = "grey") + theme_bw()
ABbiv

vbAB <- ggplot(AB) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=ABmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=ABfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=ABGroint, aes(x = age, y = fit ))+ 
  xlim(0, 10) + ylim(0, 400)+
  geom_ribbon(data=ABGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus blochii") + theme_classic()
vbAB



#Acanthurus guttatus growth
AGfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=AGfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
AGmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=AGmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
AGGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AG, start = c(Linf = 300, K = 0.5), algorithm = "port")
AGfGro
AGmGro
AGGro



new.grodatAG <- data.frame(age=seq(0, 20, by = .25))

AGfGroboot <- nlsBoot(AGfGro, niter = 999)
AGfGroint <- as_tibble(nlsBootPredict(AGfGroboot, newdata = new.grodatAG, interval = "confidence")) %>% 
  mutate(age = new.grodatAG$age)
colnames(AGfGroint) <- c("fit", "lower", "upper", "age")

AGmGroboot <- nlsBoot(AGmGro, niter = 999)
AGmGroint <- as_tibble(nlsBootPredict(AGmGroboot, newdata = new.grodatAG, interval = "confidence")) %>% 
  mutate(age = new.grodatAG$age)
colnames(AGmGroint) <- c("fit", "lower", "upper", "age")

AGGroboot <- nlsBoot(AGGro, niter = 999)
AGGroint <- as_tibble(nlsBootPredict(AGGroboot, newdata = new.grodatAG, interval = "confidence")) %>% 
  mutate(age = new.grodatAG$age)
colnames(AGGroint) <- c("fit", "lower", "upper", "age")

estsAGm <- AGmGroboot$coefboot
estsAGf <- AGfGroboot$coefboot

plot(Linf~K, data = estsAGm, main = "Males")
plot(Linf~K, data = estsAGf, main = "Females")
estsAGfdf <- as.data.frame(estsAGf)
estsAGmdf <- as.data.frame(estsAGm)
AGfplot <- ggplot(estsAGfdf, aes(K, Linf)) + geom_point()
AGbiv <- AGfplot + geom_point(data = estsAGmdf, aes(K, Linf), color = "grey") + theme_bw()
AGbiv

vbAG <- ggplot(AG) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=AGmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=AGfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=AGGroint, aes(x = age, y = fit ))+ 
  xlim(0, 21) + ylim(0, 250)+
  geom_ribbon(data=AGGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus guttatus") + theme_classic()
vbAG



#Acanthurus lineatus growth
ALfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ALfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
ALmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=ALmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
ALGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AL, start = c(Linf = 300, K = 0.5), algorithm = "port")
ALfGro
ALmGro
ALGro



new.grodatAL <- data.frame(age=seq(0, 20, by = .25))

ALfGroboot <- nlsBoot(ALfGro, niter = 999)
ALfGroint <- as_tibble(nlsBootPredict(ALfGroboot, newdata = new.grodatAL, interval = "confidence")) %>% 
  mutate(age = new.grodatAL$age)
colnames(ALfGroint) <- c("fit", "lower", "upper", "age")

ALmGroboot <- nlsBoot(ALmGro, niter = 999)
ALmGroint <- as_tibble(nlsBootPredict(ALmGroboot, newdata = new.grodatAL, interval = "confidence")) %>% 
  mutate(age = new.grodatAL$age)
colnames(ALmGroint) <- c("fit", "lower", "upper", "age")

ALGroboot <- nlsBoot(ALGro, niter = 999)
ALGroint <- as_tibble(nlsBootPredict(ALGroboot, newdata = new.grodatAL, interval = "confidence")) %>% 
  mutate(age = new.grodatAL$age)
colnames(ALGroint) <- c("fit", "lower", "upper", "age")

estsALm <- ALmGroboot$coefboot
estsALf <- ALfGroboot$coefboot

plot(Linf~K, data = estsALm, main = "Males")
plot(Linf~K, data = estsALf, main = "Females")
estsALfdf <- as.data.frame(estsALf)
estsALmdf <- as.data.frame(estsALm)
ALfplot <- ggplot(estsALfdf, aes(K, Linf)) + geom_point()
ALbiv <- ALfplot + geom_point(data = estsALmdf, aes(K, Linf), color = "grey") + theme_bw()
ALbiv

vbAL <- ggplot(AL) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=ALmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=ALfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=ALGroint, aes(x = age, y = fit ))+ 
  xlim(0, 20) + ylim(0, 250)+
  geom_ribbon(data=ALGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus lineatus") + theme_classic()
vbAL


#Acanthurus olivaceus growth
AOfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=AOfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
AOmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=AOmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
AOGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AO, start = c(Linf = 300, K = 0.5), algorithm = "port")
AOfGro
AOmGro
AOGro



new.grodatAO <- data.frame(age=seq(0, 20, by = .25))

AOfGroboot <- nlsBoot(AOfGro, niter = 999)
AOfGroint <- as_tibble(nlsBootPredict(AOfGroboot, newdata = new.grodatAO, interval = "confidence")) %>% 
  mutate(age = new.grodatAO$age)
colnames(AOfGroint) <- c("fit", "lower", "upper", "age")

AOmGroboot <- nlsBoot(AOmGro, niter = 999)
AOmGroint <- as_tibble(nlsBootPredict(AOmGroboot, newdata = new.grodatAO, interval = "confidence")) %>% 
  mutate(age = new.grodatAO$age)
colnames(AOmGroint) <- c("fit", "lower", "upper", "age")

AOGroboot <- nlsBoot(AOGro, niter = 999)
AOGroint <- as_tibble(nlsBootPredict(AOGroboot, newdata = new.grodatAO, interval = "confidence")) %>% 
  mutate(age = new.grodatAO$age)
colnames(AOGroint) <- c("fit", "lower", "upper", "age")

estsAOm <- AOmGroboot$coefboot
estsAOf <- AOfGroboot$coefboot

plot(Linf~K, data = estsAOm, main = "Males")
plot(Linf~K, data = estsAOf, main = "Females")
estsAOfdf <- as.data.frame(estsAOf)
estsAOmdf <- as.data.frame(estsAOm)
AOfplot <- ggplot(estsAOfdf, aes(K, Linf)) + geom_point()
AObiv <- AOfplot + geom_point(data = estsAOmdf, aes(K, Linf), color = "grey") + theme_bw()
AObiv

vbAO <- ggplot(AO) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=AOmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=AOfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=AOGroint, aes(x = age, y = fit ))+ 
  xlim(0, 15) + ylim(0, 250)+
  geom_ribbon(data=AOGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus olivaceus") + theme_classic()
vbAO



#Acanthurus xanthopterus growth
AXfGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=AXfgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
AXmGro <- nls(length ~ Linf-(Linf-35)*exp(-K*age), data=AXmgrodat, start = c(Linf = 300, K = 0.5), algorithm = "port")
AXGro <- nls(length ~Linf-(Linf-35)*exp(-K*age), data = AX, start = c(Linf = 300, K = 0.5), algorithm = "port")
AXfGro
AXmGro
AXGro



new.grodatAX <- data.frame(age=seq(0, 20, by = .25))

AXfGroboot <- nlsBoot(AXfGro, niter = 999)
AXfGroint <- as_tibble(nlsBootPredict(AXfGroboot, newdata = new.grodatAX, interval = "confidence")) %>% 
  mutate(age = new.grodatAX$age)
colnames(AXfGroint) <- c("fit", "lower", "upper", "age")

AXmGroboot <- nlsBoot(AXmGro, niter = 999)
AXmGroint <- as_tibble(nlsBootPredict(AXmGroboot, newdata = new.grodatAX, interval = "confidence")) %>% 
  mutate(age = new.grodatAX$age)
colnames(AXmGroint) <- c("fit", "lower", "upper", "age")

AXGroboot <- nlsBoot(AXGro, niter = 999)
AXGroint <- as_tibble(nlsBootPredict(AXGroboot, newdata = new.grodatAX, interval = "confidence")) %>% 
  mutate(age = new.grodatAX$age)
colnames(AXGroint) <- c("fit", "lower", "upper", "age")

estsAXm <- AXmGroboot$coefboot
estsAXf <- AXfGroboot$coefboot

plot(Linf~K, data = estsAXm, main = "Males")
plot(Linf~K, data = estsAXf, main = "Females")
estsAXfdf <- as.data.frame(estsAXf)
estsAXmdf <- as.data.frame(estsAXm)
AXfplot <- ggplot(estsAXfdf, aes(K, Linf)) + geom_point()
AXbiv <- AXfplot + geom_point(data = estsAXmdf, aes(K, Linf), color = "grey") + theme_bw()
AXbiv

vbAX <- ggplot(AX) + geom_point(aes(x=age, y=length, color = finalsex, shape = finalsex), size=2, alpha = 0.5) + xlab("Age (yrs)") + ylab("Fork length (mm)")  + scale_shape_manual(breaks = c("female", "male", "NA"), values=c(19, 19, 4)) + scale_color_manual(breaks = c("female", "male", "NA"), values=c("darkgrey", "black", "black")) +
  geom_line(data=AXmGroint, aes(x = age, y = fit ), linetype = "longdash")+
  geom_line(data=AXfGroint, aes(x = age, y = fit ), linetype = "dotted")+ 
  geom_line(data=AXGroint, aes(x = age, y = fit ))+ 
  xlim(0, 12) + ylim(0, 450)+
  geom_ribbon(data=AXGroint, aes(x=age, ymin=lower, ymax=upper), alpha=0.5, inherit.aes=F, fill="grey") + ggtitle("Acanthurus xanthopterus") + theme_classic()
vbAX


grid.arrange(vbAB, vbAG, vbAL, vbAN, vbAO, vbAT, vbAX, vbCSt, vbNT, vbNU, vbNV, vbZV, nrow = 4)



################################################
############AGE-BASED CATCH CURVES##############
################################################

library(segmented)
#create a function to calculatue the mode
mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#Naso unicornis data frame generation
max_ageNU <- 22 #max age
subNU <- NU[,-c(1,2)] #trimming data frame for ages only
subNU <- subset(subNU, subNU$age >=  mode(subNU$age)) # trimming ages below the modal age
freq_countNU <- table(subNU$age) #take a frequency count per age class
ln_freq_countNU <- as.numeric(log(freq_countNU)) #log the frequency count
above_modeNU <- rbind(data.frame(freq_countNU,
                                 ln_freq_countNU,
                                 age = sort(unique(subNU$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.NU <- lm(ln_freq_countNU ~ age, data = above_modeNU)
segmented.NU <- segmented(lin.mod.NU, seg.Z = ~age, psi=14)

plot(above_modeNU$age,above_modeNU$ln_freq_countNU, pch=16)
plot(segmented.NU, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.NU)

NUseg <- as.data.frame(segmented.NU[["fitted.values"]])
NUseg$x <- above_modeNU$age
colnames(NUseg) <- c("fit", "x")

NUmort <-ggplot(above_modeNU, aes(x = age, y = ln_freq_countNU)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 25) + ylim(-2, 5)+ggtitle("Naso unicornis")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=NUseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.NU, segmented.NU) 


#Acanthurus nigricauda data frame generation
max_ageAN <- 17 #max age
subAN <- AN[,-c(1,2)] #trimming data frame for ages only
subAN <- subset(subAN, subAN$age >=  mode(subAN$age)) # trimming ages below the modal age
freq_countAN <- table(subAN$age) #take a frequency count per age class
ln_freq_countAN <- as.numeric(log(freq_countAN)) #log the frequency count
above_modeAN <- rbind(data.frame(freq_countAN,
                                 ln_freq_countAN,
                                 age = sort(unique(subAN$age)))) #bind all of this in a data frame

#For surgeonfishes, testing linear and biphasic mortality patterns
lin.mod.AN <- lm(ln_freq_countAN ~ age, data = above_modeAN)
segmented.AN <- segmented(lin.mod.AN, seg.Z = ~age, psi=9)

plot(above_modeAN$age,above_modeAN$ln_freq_countAN, pch=16)
plot(segmented.AN, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AN)

ANseg <- as.data.frame(segmented.AN[["fitted.values"]])
ANseg$x <- above_modeAN$age
colnames(ANseg) <- c("fit", "x")

ANmort <-ggplot(above_modeAN, aes(x = age, y = ln_freq_countAN)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 20) + ylim(-2, 5)+ggtitle("Acanthurus nigricauda")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=ANseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AN, segmented.AN) 


#Acanthurus triostegus data frame generation
max_ageAT <- 12 #max age
subAT <- AT[,-c(1,2)] #trimming data frame for ages only
subAT <- subset(subAT, subAT$age >=  mode(subAT$age)) # trimming ages below the modal age
freq_countAT <- table(subAT$age) #take a frequency count per age class
ln_freq_countAT <- as.numeric(log(freq_countAT)) #log the frequency count
above_modeAT <- rbind(data.frame(freq_countAT,
                                 ln_freq_countAT,
                                 age = sort(unique(subAT$age)))) #bind all of this in a data frame

#For surgeonfishes, testing linear and biphasic mortality patterns
lin.mod.AT <- lm(ln_freq_countAT ~ age, data = above_modeAT)
segmented.AT <- segmented(lin.mod.AT, seg.Z = ~age, psi=6)

plot(above_modeAT$age,above_modeAT$ln_freq_countAT, pch=16)
plot(segmented.AT, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AT)

ATseg <- as.data.frame(segmented.AT[["fitted.values"]])
ATseg$x <- above_modeAT$age
colnames(ATseg) <- c("fit", "x")

ATmort <-ggplot(above_modeAT, aes(x = age, y = ln_freq_countAT)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 15) + ylim(-2, 5)+ggtitle("Acanthurus triostegus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=ATseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AT, segmented.AT) 



#Naso tonganus data frame generation
max_ageNT <- 21 #max age
subNT <- NT[,-c(1,2)] #trimming data frame for ages only
subNT <- subset(subNT, subNT$age >=  mode(subNT$age)) # trimming ages below the modal age
freq_countNT <- table(subNT$age) #take a frequency count per age class
ln_freq_countNT <- as.numeric(log(freq_countNT)) #log the frequency count
above_modeNT <- rbind(data.frame(freq_countNT,
                                 ln_freq_countNT,
                                 age = sort(unique(subNT$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.NT <- lm(ln_freq_countNT ~ age, data = above_modeNT)
segmented.NT <- segmented(lin.mod.NT, seg.Z = ~age, psi=14)

plot(above_modeNT$age,above_modeNT$ln_freq_countNT, pch=16)
plot(segmented.NT, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.NT)

NTseg <- as.data.frame(segmented.NT[["fitted.values"]])
NTseg$x <- above_modeNT$age
colnames(NTseg) <- c("fit", "x")

NTmort <-ggplot(above_modeNT, aes(x = age, y = ln_freq_countNT)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 22) + ylim(-2, 4)+ggtitle("Naso tonganus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=NTseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.NT, segmented.NT) 




#Naso vlamingii data frame generation
max_ageNV <- 9 #max age
subNV <- NV[,-c(1,2)] #trimming data frame for ages only
subNV <- subset(subNV, subNV$age >=  mode(subNV$age)) # trimming ages below the modal age
freq_countNV <- table(subNV$age) #take a frequency count per age class
ln_freq_countNV <- as.numeric(log(freq_countNV)) #log the frequency count
above_modeNV <- rbind(data.frame(freq_countNV,
                                 ln_freq_countNV,
                                 age = sort(unique(subNV$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.NV <- lm(ln_freq_countNV ~ age, data = above_modeNV)
segmented.NV <- segmented(lin.mod.NV, seg.Z = ~age, psi=4)

plot(above_modeNV$age,above_modeNV$ln_freq_countNV, pch=16)
plot(segmented.NV, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.NV)

NVseg <- as.data.frame(segmented.NV[["fitted.values"]])
NVseg$x <- above_modeNV$age
colnames(NVseg) <- c("fit", "x")

NVmort <-ggplot(above_modeNV, aes(x = age, y = ln_freq_countNV)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 10) + ylim(-2, 4)+ggtitle("Naso vlamingii")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=NVseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.NV, segmented.NV) 



#Zebrasoma velifer data frame generation
subZV <- ZV[,-c(1,2)] #trimming data frame for ages only
subZV <- subset(subZV, subZV$age >=  mode(subZV$age)) # trimming ages below the modal age
freq_countZV <- table(subZV$age) #take a frequency count per age class
ln_freq_countZV <- as.numeric(log(freq_countZV)) #log the frequency count
above_modeZV <- rbind(data.frame(freq_countZV,
                                 ln_freq_countZV,
                                 age = sort(unique(subZV$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.ZV <- lm(ln_freq_countZV ~ age, data = above_modeZV)
segmented.ZV <- segmented(lin.mod.ZV, seg.Z = ~age, psi=12)

plot(above_modeZV$age,above_modeZV$ln_freq_countZV, pch=16)
plot(segmented.ZV, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.ZV)

ZVseg <- as.data.frame(segmented.ZV[["fitted.values"]])
ZVseg$x <- above_modeZV$age
colnames(ZVseg) <- c("fit", "x")

ZVmort <-ggplot(above_modeZV, aes(x = age, y = ln_freq_countZV)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 23) + ylim(-1, 2.5)+ggtitle("Zebrasoma velifer")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=ZVseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.ZV, segmented.ZV) 




#Acanthurus blochii data frame generation
subAB <- AB[,-c(1,2)] #trimming data frame for ages only
subAB <- subset(subAB, subAB$age >=  mode(subAB$age)) # trimming ages below the modal age
freq_countAB <- table(subAB$age) #take a frequency count per age class
ln_freq_countAB <- as.numeric(log(freq_countAB)) #log the frequency count
above_modeAB <- rbind(data.frame(freq_countAB,
                                 ln_freq_countAB,
                                 age = sort(unique(subAB$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.AB <- lm(ln_freq_countAB ~ age, data = above_modeAB)
segmented.AB <- segmented(lin.mod.AB, seg.Z = ~age, psi=6)

plot(above_modeAB$age,above_modeAB$ln_freq_countAB, pch=16)
plot(segmented.AB, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AB)

ABseg <- as.data.frame(segmented.AB[["fitted.values"]])
ABseg$x <- above_modeAB$age
colnames(ABseg) <- c("fit", "x")

ABmort <-ggplot(above_modeAB, aes(x = age, y = ln_freq_countAB)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 9) + ylim(-1, 2.5)+ggtitle("Acanthurus blochii")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=ABseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AB, segmented.AB) 



#Acanthurus guttatus data frame generation
subAG <- AG[,-c(1,2)] #trimming data frame for ages only
subAG <- subset(subAG, subAG$age >=  mode(subAG$age)) # trimming ages below the modal age
freq_countAG <- table(subAG$age) #take a frequency count per age class
ln_freq_countAG <- as.numeric(log(freq_countAG)) #log the frequency count
above_modeAG <- rbind(data.frame(freq_countAG,
                                 ln_freq_countAG,
                                 age = sort(unique(subAG$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.AG <- lm(ln_freq_countAG ~ age, data = above_modeAG)
segmented.AG <- segmented(lin.mod.AG, seg.Z = ~age, psi=12)

plot(above_modeAG$age,above_modeAG$ln_freq_countAG, pch=16)
plot(segmented.AG, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AG)

AGseg <- as.data.frame(segmented.AG[["fitted.values"]])
AGseg$x <- above_modeAG$age
colnames(AGseg) <- c("fit", "x")

AGmort <-ggplot(above_modeAG, aes(x = age, y = ln_freq_countAG)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 23) + ylim(-1, 3.5)+ggtitle("Acanthurus guttatus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=AGseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AG, segmented.AG) 



#Acanthurus lineatus data frame generation
subAL <- AL[,-c(1,2)] #trimming data frame for ages only
subAL <- subset(subAL, subAL$age >=  mode(subAL$age)) # trimming ages below the modal age
freq_countAL <- table(subAL$age) #take a frequency count per age class
ln_freq_countAL <- as.numeric(log(freq_countAL)) #log the frequency count
above_modeAL <- rbind(data.frame(freq_countAL,
                                 ln_freq_countAL,
                                 age = sort(unique(subAL$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.AL <- lm(ln_freq_countAL ~ age, data = above_modeAL)
segmented.AL <- segmented(lin.mod.AL, seg.Z = ~age, psi=12)

plot(above_modeAL$age,above_modeAL$ln_freq_countAL, pch=16)
plot(segmented.AL, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AL)

ALseg <- as.data.frame(segmented.AL[["fitted.values"]])
ALseg$x <- above_modeAL$age
colnames(ALseg) <- c("fit", "x")

ALmort <-ggplot(above_modeAL, aes(x = age, y = ln_freq_countAL)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 20) + ylim(-1, 3.5)+ggtitle("Acanthurus lineatus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=ALseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AL, segmented.AL) 



#Acanthurus olivaceus data frame generation
subAO <- AO[,-c(1,2)] #trimming data frame for ages only
subAO <- subset(subAO, subAO$age >=  mode(subAO$age)) # trimming ages below the modal age
freq_countAO <- table(subAO$age) #take a frequency count per age class
ln_freq_countAO <- as.numeric(log(freq_countAO)) #log the frequency count
above_modeAO <- rbind(data.frame(freq_countAO,
                                 ln_freq_countAO,
                                 age = sort(unique(subAO$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.AO <- lm(ln_freq_countAO ~ age, data = above_modeAO)
segmented.AO <- segmented(lin.mod.AO, seg.Z = ~age, psi=12)

plot(above_modeAO$age,above_modeAO$ln_freq_countAO, pch=16)
plot(segmented.AO, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AO)

AOseg <- as.data.frame(segmented.AO[["fitted.values"]])
AOseg$x <- above_modeAO$age
colnames(AOseg) <- c("fit", "x")

AOmort <-ggplot(above_modeAO, aes(x = age, y = ln_freq_countAO)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 15) + ylim(-1, 3.5)+ggtitle("Acanthurus olivaceus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=AOseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AO, segmented.AO) 




#Acanthurus xanthopterus data frame generation
subAX <- AX[,-c(1,2)] #trimming data frame for ages only
subAX <- subset(subAX, subAX$age >=  mode(subAX$age)) # trimming ages below the modal age
freq_countAX <- table(subAX$age) #take a frequency count per age class
ln_freq_countAX <- as.numeric(log(freq_countAX)) #log the frequency count
above_modeAX <- rbind(data.frame(freq_countAX,
                                 ln_freq_countAX,
                                 age = sort(unique(subAX$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.AX <- lm(ln_freq_countAX ~ age, data = above_modeAX)
segmented.AX <- segmented(lin.mod.AX, seg.Z = ~age, psi=7)

plot(above_modeAX$age,above_modeAX$ln_freq_countAX, pch=16)
plot(segmented.AX, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.AX)

AXseg <- as.data.frame(segmented.AX[["fitted.values"]])
AXseg$x <- above_modeAX$age
colnames(AXseg) <- c("fit", "x")

AXmort <-ggplot(above_modeAX, aes(x = age, y = ln_freq_countAX)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 15) + ylim(-1, 3)+ggtitle("Acanthurus xanthopterus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=AXseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.AX, segmented.AX) 




#Ctenochaetus striatus data frame generation
subCSt <- CSt[,-c(1,2)] #trimming data frame for ages only
subCSt <- subset(subCSt, subCSt$age >=  mode(subCSt$age)) # trimming ages below the modal age
freq_countCSt <- table(subCSt$age) #take a frequency count per age class
ln_freq_countCSt <- as.numeric(log(freq_countCSt)) #log the frequency count
above_modeCSt <- rbind(data.frame(freq_countCSt,
                                 ln_freq_countCSt,
                                 age = sort(unique(subCSt$age)))) #bind all of this in a data frame

#testing linear and biphasic mortality patterns
lin.mod.CSt <- lm(ln_freq_countCSt ~ age, data = above_modeCSt)
segmented.CSt <- segmented(lin.mod.CSt, seg.Z = ~age, psi=12)

plot(above_modeCSt$age,above_modeCSt$ln_freq_countCSt, pch=16)
plot(segmented.CSt, conf.level=0.95, shade = TRUE, add=T)
abline(lin.mod.CSt)

CStseg <- as.data.frame(segmented.CSt[["fitted.values"]])
CStseg$x <- above_modeCSt$age
colnames(CStseg) <- c("fit", "x")

CStmort <-ggplot(above_modeCSt, aes(x = age, y = ln_freq_countCSt)) + xlab("Age (yrs)") + ylab("LN Frequency")  +
  geom_point(col = "black") + xlim(0, 20) + ylim(-1, 3.5)+ggtitle("Ctenochaetus striatus")+
  stat_smooth(method = "lm", col = "black", size = .5, level = 0.95) + theme_classic() +
  geom_hline(yintercept = 0, color = "black", size = .5) +
  geom_line(data=CStseg, aes(x = x, y = fit), size = .75)

anova(lin.mod.CSt, segmented.CSt) 


grid.arrange(ABmort, AGmort, ALmort, ANmort, AOmort, ATmort, AXmort, CStmort, NTmort, NUmort, NVmort, ZVmort, nrow = 4)



################################################
############OTOLITH WEIGHT BY AGE###############
################################################

NUotwt <- ggplot(NUdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "loess", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Naso unicornis")

ANotwt <- ggplot(ANdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "loess", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus nigricauda")

ATotwt <- ggplot(ATdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus triostegus")

ABotwt <- ggplot(ABdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus blochii")

AGotwt <- ggplot(AGdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus guttatus")

ALotwt <- ggplot(ALdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "loess", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus lineatus")

AOotwt <- ggplot(AOdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus olivaceus")

AXotwt <- ggplot(AXdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Acanthurus xanthopterus")

CStotwt <- ggplot(CStdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Ctenochaetus striatus")

NTotwt <- ggplot(NTdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Naso tonganus")

NVotwt <- ggplot(NVdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Naso vlamingii")

ZVotwt <- ggplot(ZVdat, aes(x = otolith_w, y = age)) + xlab("Otolith weight (g)") + ylab("Age (years)") + 
  geom_point(col = "black", alpha = 0.25) + stat_smooth(method = "lm", col = "black", size = 0.5, level = 0.95) + theme_classic() +
  ggtitle("Zebrasoma velifer")

grid.arrange(ABotwt, AGotwt, ALotwt, ANotwt, AOotwt, ATotwt, AXotwt, CStotwt, NTotwt, NUotwt, NVotwt, ZVotwt, nrow = 4)
