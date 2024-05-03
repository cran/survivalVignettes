## ----include = FALSE------------------------------------------------
options(continue="  ", width=70, prompt="> ")
options(contrasts=c("contr.treatment", "contr.poly")) #ensure default
options(show.significant.stars = FALSE) #statistical intelligence

knitr::opts_chunk$set(
  collapse = TRUE, warning=FALSE, error=FALSE, tidy=FALSE,
  fig.asp=.75, fig.align="center",
  comment = "#>", highlight=TRUE, echo=TRUE, prompt=FALSE,
  fig.width=7, fig.height=5.5
)

library(survival)
library(geepack)
library(survey)

## ----lung1----------------------------------------------------------
lfit0 <- survfit(Surv(time/365.25, status) ~ 1, data=lung)
lfit1 <- survfit(Surv(time/365.25, status) ~ ph.ecog, data=lung)
print(lfit1, rmean=2.5)
plot(lfit1, lty=1:4, xlab="Years since enrollment", ylab="Survival")
legend("topright", c("PS 0", "PS 1", "PS 2", "PS 3"), 
       lty=1:4, bty='n')

## ----lung2----------------------------------------------------------
pmean <- pseudo(lfit0, times=2.5, type="RMST")
ldata <- data.frame(lung, pmean= c(pmean), 
                    id=1:nrow(lung))  
afit1 <- lm(pmean ~ ph.ecog + sex + age, data=ldata)
round(summary(afit1)$coefficients, 3)

## ----geese----------------------------------------------------------
afit2 <- geese(pmean ~ ph.ecog + sex + age, data=ldata, 
               subset = (!is.na(ph.ecog)))
round(summary(afit2)$mean, 3)


ldesign <- svydesign(data=ldata, id= ~id, weights=NULL,
              variables= ~ . -id)
afit3 <- svyglm(pmean ~ ph.ecog + sex + age, design=ldesign)
round(summary(afit3)$coefficients, 3)

## ----strata---------------------------------------------------------
# There is 1 missing value for inst so pseudo by default will have 227 values instead of 228
# When creating pmean4, you need to put the values in the correct spot.
# Specifying na.exclude expands the results to have the proper length

lfit4 <- survfit(Surv(time/365.25, status) ~ inst, data=lung, na.action=na.exclude)
pmean4 <- pseudo(lfit4, time=2.5, type="auc")
ldata$pmean4 <- c(pmean4)
afit4a <- lm(pmean4 ~ ph.ecog + sex + age, data=ldata)

# survey package does not allow a missing strata, and one obs in the lung data
# has institution missing
temp <- subset(ldata, !is.na(inst))
ldesign4 <- svydesign(data=temp,  id= ~id, weights=NULL, strata= ~inst,
                      variables= ~ .-id -inst)
afit4b <- svyglm(pmean4 ~ ph.ecog + sex + age, design=ldesign4)
round(summary(afit4a)$coefficients, 3)
round(summary(afit4b)$coefficients[1:4,], 3)

## ----test1----------------------------------------------------------
with(aml, Surv(time, status))
fit1 <- survfit(Surv(time, status) ~1, aml)
rr1 <- resid(fit1, times=c(12,24))
pv1 <- pseudo(fit1, times= c(12, 24))
round(pv1[1:8,], 4)

## ----test2----------------------------------------------------------
pdata <- pseudo(fit1, times= c(12, 24), data.frame=TRUE)
lfit1 <- lm(pseudo ~1, pdata, subset= (time==12))
lfit2 <- lm(pseudo ~1, pdata, subset= (time==24))

summary(lfit1)$coefficients
summary(lfit2)$coefficients

summary(fit1, time=c(12,24))

## ----test3----------------------------------------------------------
cfit3 <- coxph(Surv(time, status) ~x, data=aml)
cfit3

pdata <- cbind(pdata, x= aml$x)   # original data + pseudo-values
lfit3 <- lm(pseudo ~ x + factor(time), data=pdata)
summary(lfit3)$coefficients

sfit3 <- survfit(cfit3, newdata=list(x= c("Maintained", "Nonmaintained")))
summary(sfit3, times=c(12,24))

## ----test4, echo=FALSE----------------------------------------------
temp <- matrix(0, 2, 4, dimnames=list(c("Maintainance", "Nonmaintainance"),
                                      c("Cox, 12", "Pseudo, 12", 
                                        "Cox, 24", "Pseudo, 24")))
temp[,c(1,3)] <- 1- t(summary(sfit3, times=c(12,24))$surv)
dummy <-  expand.grid(x= c("Maintained", "Nonmaintained"),
                      time=c(12, 24))
temp[,c(2,4)] <- 1- predict(lfit3, newdata= dummy)
temp <- rbind(temp, "Difference" = temp[2,]- temp[1,])
round(temp,2)

## ----svy------------------------------------------------------------
# Naive lm
summary(lfit3)$coefficients
#
#survey
# When an id is constructed rather than supplied, pseudo() purposely gives
#  it a non-standard name to avoid confusion with any user's variable names.
#  But that name is a PITA to use in formulas
flag <- toupper(names(pdata))=='(ID)'
pdata$id <- pdata[,flag] 
pdesign <- svydesign(id= ~ id, varibles= ~ pseudo + x + time, weights=NULL,
                     data=pdata)
sfit <- svyglm(pseudo ~ x + factor(time), design=pdesign)
summary(sfit)$coefficients
#
# GEE 1
gfit1 <- geese(pseudo ~ x + factor(time), data=pdata, id= id) # wrong answer
summary(gfit1)$mean
#
# GEE 2
pdata <- pdata[order(pdata$id),]  # group id values to be together
gfit2 <- geese(pseudo ~ x + factor(time), data=pdata, id= id)
summary(gfit2)$mean

## ----cdata----------------------------------------------------------
cdata <- subset(colon, etype==1, -etype)  # time to recurrence
temp  <- subset(colon, etype==2)
cdata$status <- pmax(cdata$status, temp$status)

trtsurv <- survfit(Surv(time, status) ~ rx, cdata, id=id)
plot(trtsurv, fun="event", xscale= 365.25, col = 1:3, lwd=2,
     xlab= "Years from randomization", ylab="Treatment failure")
legend(5*365, .25, levels(cdata$rx), col=1:3, lwd=2, bty='n')

ccox <- coxph(Surv(time, status) ~ rx + sex + age + extent + node4+
                  obstruct + perfor + adhere,  data=cdata)

## ----colon2---------------------------------------------------------
csurv <- survfit(Surv(time, status)  ~1, data=cdata, id=id)
cptemp <- pseudo(csurv, time= 1:6 * 365.25, data.frame=TRUE)
cptemp$pdeath <- 1- cptemp$pseudo
cptemp$year <- round(cptemp$time/365.25)

cpdat <- merge(cptemp, 
               subset(x=cdata,select=c(id, rx, extent, node4)), 
               by="id")

## -------------------------------------------------------------------
hist(cpdat$pdeath, nclass=50, main=NULL)

## ----c6-------------------------------------------------------------
temp1 <- summary(trtsurv, times= 1:6*365.25)
temp2 <- matrix(temp1$surv, nrow=6)
temp3 <- rbind("Obs vs 5FU"= temp2[,1]- temp2[,2],
            "Lev vs 5FU"= temp2[,3]- temp2[,2]) 

cat("absolute difference\n")
round(temp3*100, 1)
#
cat("logit difference\n") 
logit <- function(x) log(1/(1-x))
temp4 <- rbind("Obs vs 5FU"= logit(temp2[,1])- logit(temp2[,2]),
            "Lev vs 5FU"= logit(temp2[,3])- logit(temp2[,2]))
round(temp4*100, 1)

## ----times----------------------------------------------------------
onetime <- matrix(0, 6, 4,
                  dimnames=list(paste("Time", 1:6), 
                                c("coefficient", "se.glm", "coef", "se.survey")))
cpdesign <- svydesign(id= ~id, weights=NULL, data=cpdat,
                      variables = ~ . - id)
for (i in 1:6) {
    fit1 <- lm(pdeath ~ rx + extent + node4, data=cpdat, subset= (year==i))
    onetime[i,1:2] <- summary(fit1)$coefficients[3,1:2]
    sfit1 <- svyglm(pdeath ~ rx+ extent + node4, design=cpdesign, 
                    subset= (year==i))
    onetime[i,3:4] <- summary(sfit1)$coefficients[3, 1:2]
}

onetime <- cbind(onetime, "coef/se"= onetime[,3]/onetime[,4])
round(onetime[, c(1,2,4,5)],3)

## ----twotime--------------------------------------------------------
pfun <- function(x,d=2) printCoefmat(x, digits=d, P.values=TRUE,
                                 has.Pvalue=TRUE, signif.stars= FALSE)

sfit1 <- svyglm(pdeath ~ rx + extent + node4,
                        design= cpdesign, subset=(year==4))
pfun(summary(sfit1)$coefficients[2:5,])

#
sfit3 <- svyglm(pdeath ~ rx + extent + node4 + factor(year), 
                        design= cpdesign, 
                subset= (year==2 | year==4 | year==6))
pfun(summary(sfit3)$coefficients[2:5,],)

sfit6 <- svyglm(pdeath ~ rx + extent + node4 + factor(year),
                design= cpdesign)
pfun(summary(sfit6)$coefficients)

## ----catch----------------------------------------------------------
tryCatch( svyglm(pdeath ~rx + extent + node4 + factor(year), design=cpdesign,
              family= gaussian(link = "logit")),
         error= function(e) e)

## ----link1----------------------------------------------------------
sfit6b <- svyglm(pdeath ~rx + extent + node4 + factor(year), design=cpdesign,
              family= gaussian(link = "logit"), 
              etastart= pmax(.05, pmin(.95, cpdat$pdeath)))

## ----blogit---------------------------------------------------------
blogit <- function(edge=.05) {
    new <- make.link("logit")
    new$linkfun <- function(mu) { 
        x <- (pmax(edge, pmin(mu, 1-edge)))
        log(x/(1-x))
    }
    new
}
sfit3c <- svyglm(pdeath ~rx + extent + node4 + factor(year), design=cpdesign,
              family= gaussian(link = blogit()), subset=(year %in% c(2,4,6)))
pfun(summary(sfit3c)$coefficients)

## ----loglog---------------------------------------------------------
bcloglog <- function(edge=.05) {
    new <- make.link("cloglog")
    new$linkfun <- function(mu) {
        x <- (pmax(edge, pmin(mu, 1-edge)))
        log(-log(x))
    }
    new$name <- "bcloglog"
    new
}
sfit3d <- svyglm(pseudo ~rx + extent + node4 + factor(year), design=cpdesign,
              family= gaussian(link = bcloglog()), subset=(year %in% c(2,4,6)))
pfun(summary(sfit3d)$coefficients)
cfit <- coxph(Surv(time, status) ~ rx + extent + node4, cdata)                 
round(summary(cfit)$coefficients, 3)

## ----mgus2cr--------------------------------------------------------
mdata <- mgus2
mdata$etime <- with(mdata, ifelse(pstat==1, ptime, futime))
mdata$event <- factor(with(mdata, ifelse(pstat==1, 1, 2*death)), 0:2,
                      c("censor", "PCM", "Death"))
mdata$age10 <- mdata$age/10  # age in decades

msurv <- survfit(Surv(etime, event) ~1, data=mdata, id=id)
plot(msurv, lty=1:2, xscale=12, 
     xlab="Years from Diagnosis", ylab="Probability in state")
text(c(345, 340), c(.18, .73), c("PCM", "Death w/o PCM"))

## ----mgus2p---------------------------------------------------------
mps <- pseudo(msurv, times=12*(1:30), type="pstate", data.frame=TRUE)
mdata6 <- merge(mdata, subset(mps, state== "PCM"), by="id")
mdata6$year <- mdata6$time/12

mdesign <- svydesign(data=mdata6, id=~id, weights=NULL, variables= ~. -id)
mpfit1 <- svyglm(pseudo ~ age10 + sex + mspike + factor(year), design = mdesign,
                 family = gaussian(link =bcloglog()))
mpfit2 <- svyglm(pseudo ~ age10 + sex + mspike + factor(year), design = mdesign,
                 subset = (year %in% (3* 1:10)),
                 family = gaussian(link =bcloglog()))
mpfit3 <- svyglm(pseudo ~ age10 + sex + mspike + factor(year), design = mdesign,
                 subset = (year %in% (6* 1:5)),
                 family = gaussian(link =bcloglog()))

fgdata <- finegray(Surv(etime, event) ~ .,  mdata, etype="PCM")
fgfit <-  coxph(Surv(fgstart, fgstop, fgstatus) ~ age10 + sex + mspike, 
                data=fgdata, weights = fgwt)

temp <- cbind(summary(fgfit)$coef[1:3, c(1,3)],
              summary(mpfit1)$coef[2:4, 1:2],
              summary(mpfit2)$coef[2:4,1:2],
              summary(mpfit3)$coef[2:4,1:2])
colnames(temp) <- c("FGcoef", "FGstd", "Coef30", "Std30", 
                    "Coef10", "Std10", "Coef5", "Std5")

round(temp[,c(1,3,5,7, 2,4,6,8)], 3)

