library(readxl)
library(BBmisc)
my_data<- read_excel("Mercury Final.xlsx")
View(my_data)

keeps<-c("Site", "Depth", "Length", "LogLength", 
         "Weight", "Maturity", "Morph", "Adjdel13C", 
         "del15N", "del34S","WetHg (ppm)")
newdata<- my_data[keeps]
names(newdata)[11]<- "WetHg"
View(newdata)
attach(newdata)

hist(Length)
shapiro.test(Length)

hist(LogLength)
shapiro.test(LogLength)

hist(WetHg)
shapiro.test(WetHg)

newdata[,11] <- log(newdata[11])
names(newdata)[11]<-"LWetHg"

hist(LWetHg)
shapiro.test(LWetHg)

shapiro.test(log(del34S))
shapiro.test(log(del15N))
shapiro.test(Adjdel13C)

hist(Weight)
shapiro.test(Weight)
hist(log(Weight))
shapiro.test(log(Weight))
nweight<-normalize(Weight, method='range', range=c(0,1))
nweight<-nweight[nweight>0.0000000]
lnweight<-log(nweight)
hist(lnweight)
shapiro.test(lnweight)

newdata[,10] <- log(newdata[10])
names(newdata)[10]<-"Ldel34S"

newdata[,9] <- log(newdata[9])
names(newdata)[9]<-"Ldel15N"

newdata[,5] <- log(normalize(newdata[5], method='range', range=c(0,1)))
names(newdata)[5]<-"LWeight"
hist(newdata$LWeight)

newdata<- subset(newdata, LWeight!="-Inf")
newdata<- subset(newdata, Maturity!="Unknown")
View(newdata)
attach(newdata)

library(car)
library(multcompView)

class(Maturity)
maturmod<-lm(LWetHg~Maturity, data=newdata)
maturmod2<-aov(LWetHg~Maturity, data=newdata)
plot(maturmod2)
shapiro.test(residuals(maturmod2))

MTukey<- TukeyHSD(maturmod2)
plot(TukeyHSD(maturmod2))
plot(MTukey, las=1 , col="brown")
MTukey
maturmod2<- Anova(maturmod2, type = "III")
#Ho: mean Lwethg is the same across all maturities

class(Depth)
depthmod<-aov(LWetHg~Depth)
depthmodtype3<-Anova(depthmod, type="III")
summary(depthmod)
boxplot(LWetHg~Depth)
plot(depthmod)
shapiro.test(residuals(depthmod))
DTukey<-TukeyHSD(depthmod)
DTukey

class(Morph)
morphmod<-aov(LWetHg~Morph)
summary(morphmod)
boxplot(LWetHg~Morph)
plot(morphmod)
shapiro.test(residuals(morphmod))
plot(TukeyHSD(morphmod))

Morph<-as.factor(Morph)
levels(Morph)
Depth<-as.factor(Depth)
Maturity<-as.character(Maturity)


library(dplyr)
aicdata<-select(newdata,-Site)
aicdata<-select(aicdata,-Length)
attach(aicdata)
View(aicdata)

entiremod <- lm(LWetHg~., data=aicdata)
summary(entiremod)
plot(entiremod)
shapiro.test(residuals(entiremod))
library(car)
vif(entiremod)

newmod<-lm(LWetHg~Depth + LogLength + Adjdel13C + Ldel15N + Ldel34S + Maturity + Morph, 
           data=aicdata)
summary(newmod)
vif(newmod)

categmod<-lm(LWetHg~Maturity + Depth + Morph, data=aicdata)
summary(categmod)

#AIC analysis
mod1 <- lm(LWetHg~Ldel15N, data = aicdata)
mod2 <- lm(LWetHg~LogLength, data = aicdata)
mod3 <- lm(LWetHg~Depth, data = aicdata)
mod4 <- lm(LWetHg~LogLength + Depth, data = aicdata)
mod5 <- lm(LWetHg~Ldel15N + Ldel34S, data = aicdata)
mod6 <- lm(LWetHg~Ldel15N + Ldel34S + Depth, data = aicdata)
mod7 <- lm(LWetHg~Ldel15N + Ldel34S + Depth + LogLength, data=aicdata)
mod8 <- lm(LWetHg~Ldel15N + Ldel34S + Depth + LogLength + Morph, data=aicdata)
AICc(mod8)
library(MuMIn)

allaic <- c(AICc(mod1), AICc(mod2), AICc(mod3), AICc(mod4), AICc(mod5), AICc(mod6), AICc(mod7), AICc(mod8))
allaic
delaic <- allaic - min(allaic)
likelihood <- exp(-0.5 * delaic)
aicweight <- likelihood/sum(likelihood)
eratio <- max(aicweight)/aicweight

aictable <- data.frame(AICc = allaic, 
                       deltaAIC = delaic,
                       ModelLikelihood = likelihood,
                       ModelWeight = aicweight,
                       EvidenceRatio = eratio)
rownames(aictable) <- c("Model 1", "Model 2","Model 3",
                        "Model 4","Model 5","Model 6","Model 7", "Model 8")
round(aictable, digits = 3)
aictable

library(MASS)

stepwiseAIC <- stepAIC(mod8, direction = "backward")

