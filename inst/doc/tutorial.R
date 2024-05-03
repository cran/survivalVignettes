## ----echo=FALSE-----------------------------------------------------
options(continue="  ", width=70, prompt="> ")
options(contrasts=c("contr.treatment", "contr.poly")) #ensure default
options(show.significant.stars = FALSE) #statistical intelligence
knitr::opts_chunk$set(
  collapse = TRUE, warning=FALSE, error=FALSE, tidy=FALSE,
  fig.asp=.75, fig.align="center",
  comment = "#>", highlight=TRUE, echo=TRUE, prompt=FALSE,
  fig.width=7, fig.height=5.5
)

library("survival")

## ----initial--------------------------------------------------------
data(aidssi, package="mstate")  # AIDS data set
data(ebmt3, package="mstate")   # transplant data set

## ----aids, fig.height=3, fig.width=3--------------------------------
oldpar <- par(mar=c(0,0,0,0))
states <- c("Event free", "AIDS", "SI")
smat <- matrix(0, 3, 3, dimnames=list(states, states))
smat[1,2] <- smat[1,3] <- 1
statefig(1:2, smat)
par(oldpar)

## ----fig-putter1, fig.cap='Correct Aalen-Johansen curves for the cumulative incidence of AIDS and SI'----
aidssi$event <- factor(aidssi$status, 0:2, c("censored", "AIDS", "SI"))

# The correct Aalen-Johansen curves
ajfit <- survfit(Surv(time, event) ~1, data = aidssi)
ajfit$transitions
plot(ajfit, xmax = 13, col = 1:2, lwd = 2,
     xlab = "Years from HIV infection", ylab = "Probability")
legend(8, .2, c("AIDS", "SI"), lty = 1, lwd = 2, col = 1:2, bty = 'n')

## ----fig-T2, fig.cap="This figure (T2) shows the estimated survival curve for AIDS and probability of SI appearance, based on the naive Kaplan-Meier estimate."----
# re-create figure T2
# KM curves that censor the other endpoint (a bad idea)
bad1 <- survfit(Surv(time, event=="AIDS") ~ 1, data = aidssi)
bad2 <- survfit(Surv(time, event=="SI") ~1, data = aidssi)

plot(bad1, conf.int = FALSE, xmax = 13, 
     xlab = "Years from HIV infection", ylab = "Probability")
lines(bad2, conf.int = FALSE, fun = 'event', xmax = 13)
text(c(8,8), c(.8, .22), c("AIDS", "SI"))

## ----pstack---------------------------------------------------------
pstack <- function(fit, top=FALSE, ...) {
    temp <- survfit0(fit)   # add the point at time 0
    if (is.matrix(temp$pstate))  # usual case
        temp$pstate <- t(apply(temp$pstate, 1, cumsum))
    else if (is.array(temp$pstate)) 
        temp$pstate <- aperm(apply(temp$pstate, 1:2, cumsum), c(2,3,1))
    # this works because we don't change any other aspect of the survfit
    #  object, but only modify the probabilities.
    if (top) plot(temp, noplot="", ...)
    else plot(temp, noplot=temp$states[length(temp$states)], ...)
}

## ----fig-T3, fig.cap="This figure (T3) shows estimates of probabilities of AIDS and SI appearance, based on the naive Kaplan-Meier (grey) and on cumulative incidence/Aalen-Johansen functions (black)"----
# re-create figure T3
pstack(ajfit[c(2,1,3)], col=1, xmax=13, lwd=2, 
     xlab="Years from HIV infection", ylab="Probability")
lines(bad1, conf.int=FALSE, col="lightgray")
lines(bad2, conf.int=FALSE, fun='event', col='lightgray')
text(c(4, 8,8), c(.5, .85, .15), c("Event free", "AIDS", "SI"), col=1)

## ----fig-T4, fig.cap="This figure (T4) shows cumulative incidence curves of AIDS and SI appearance.  The cumulative incidence functions are stacked; the distance between two curves represent the probability of the different events."----
pstack(ajfit[c(2,3,1)], xmax=13, lwd=2, col=1, ylim=c(0,1),
        xlab="Years from HIV infection", ylab="Probability")
text(c(11, 11, 11), c(.2, .55, .9), c("AIDS", "SI", "Event free"))

## ----fig-cuminc, fig.cap="Cumulative incidence fit using the cumlative hazard function"----
plot(ajfit, cumhaz=TRUE, xmax=13, col=1:2, lty=2,
     xlab="Years from HIV infection", ylab="Cumulative incidence")
lines(bad1, cumhaz=TRUE, conf.int=FALSE)
lines(bad2, cumhaz=TRUE, col=2, conf.int=FALSE)

## ----cfit-----------------------------------------------------------
cfit0 <- coxph(Surv(time, event) ~ ccr5, data = aidssi, id = patnr)
print(cfit0, digits=2)

cfit1 <- coxph(Surv(time, event=="AIDS") ~ ccr5, data = aidssi)
print(cfit1, digits=2)

cfit2 <- coxph(Surv(time, event=="SI") ~ ccr5, data = aidssi)
print(cfit2, digits=2)

## ----stack----------------------------------------------------------
temp <- subset(aidssi, select= c(patnr, time, ccr5))
temp1 <- data.frame(temp, status= 1*(aidssi$event=="AIDS"), cause="AIDS")
temp2 <- data.frame(temp, status= 1*(aidssi$event=="SI"),   cause="SI")
stack <- rbind(temp1, temp2)

cfit3 <- coxph(Surv(time, status) ~ ccr5 * strata(cause), data=stack)
print(cfit3, digits=2)

## ----stack2---------------------------------------------------------
stack$ccr5.1 <- (stack$ccr5=="WM") * (stack$cause == "AIDS")
stack$ccr5.2 <- (stack$ccr5=="WM") * (stack$cause == "SI")
cfit3b <- coxph(Surv(time, status) ~ ccr5.1 + ccr5.2 + strata(cause), data = stack)
cfit3b$coef

temp <- cbind(cfit0=cfit0$loglik, cfit3=cfit3$loglik, cfit3b=cfit3b$loglik)
rownames(temp) <- c("beta=0", "beta=final")
temp

## ----common---------------------------------------------------------
common1 <- coxph(Surv(time, status) ~ ccr5 + strata(cause), data=stack)
print(common1, digits=2)

common1b <- coxph(list( Surv(time, event) ~ 1, 
                        1:2 + 1:3 ~ ccr5/common ),
                  data=aidssi, id=patnr)

## ----common2--------------------------------------------------------
common2 <- coxph(Surv(time, status) ~ ccr5, data = stack)
all.equal(common2$coef, common1$coef)

## ----common3--------------------------------------------------------
# reprise common1 and common2, using the breslow option
test1 <- coxph(Surv(time, status) ~ ccr5 + strata(cause), stack,
               ties='breslow')
test2 <- coxph(Surv(time, status) ~ ccr5, stack, ties='breslow')
all.equal(test2$loglik + test2$nevent * log(2),  test1$loglik)
all.equal(test2$coef, test1$coef)

test3 <- coxph(Surv(time, status) ~ ccr5 + cause, stack, ties='breslow')
test3
all.equal(test3$coef[1], test1$coef)

## ----aidscurve------------------------------------------------------
# re-create figure T5 in a single panel
dummy <- data.frame(ccr5=c("WW", "WM"))
pred.aj <- survfit(cfit0, newdata=dummy)
dim(pred.aj)
pred.aj$states

## ----fig-T5, fig.cap="Figure T5 showing predicted curves for AIDS and SI stratified by subjects with CCR5 wild-type (WW) and mutant (WM)."----
oldpar <- par(mfrow=c(1,2))
plot(pred.aj[,,"AIDS"], lwd=2, col=c("black", "gray"), 
     xmax=13, ylim=c(0,.5),
     xlab="Years from HIV infection", ylab="Probability of AIDS")
text(c(9.5, 10), c(.3, .1), c("WW", "WM"))
plot(pred.aj[,,"SI"], lwd=2, col=c("black", "gray"), 
      xmax=13, ylim=c(0,.5),
    xlab="Years from HIV infection", ylab="Probability of SI")
text(c(8.5, 9), c(.33, .25), c("WW", "WM"))
par(oldpar)

## ----finegray-------------------------------------------------------
fdata1 <- finegray(Surv(time, event) ~ ., data = aidssi, etype = 'AIDS')
fgfit1 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ ccr5, data = fdata1,
                weight = fgwt)
fgfit1

fdata2 <- finegray(Surv(time, event) ~., aidssi, etype="SI")
fgfit2 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ ccr5, fdata2,
                weight = fgwt)
fgfit2

## ----fig-T8, fig.cap="Reproduction of Figure T8 showing cumulative incidence functions for AIDS (left) and SI appearance (right), for CCR5 wild-type (WW) and mutant (WM), based on the Fine and Gray model."----
# re-create figure T8: Fine-Gray curves
fgsurv1<-survfit(fgfit1, newdata=dummy)
fgsurv2<-survfit(fgfit2, newdata=dummy)

oldpar <- par(mfrow=c(1,2), mar=c(4.1, 3.1, 3.1, 1)) #leave room for title
plot(fgsurv1, col=1:2, lty=c(1,1,2,2), lwd=2, xmax=13,
     ylim=c(0, .5),fun='event',
     xlab="Years from HIV infection", ylab="Probability")
title("AIDS")
plot(fgsurv2, col=1:2, lty=c(1,1,2,2), lwd=2, xmax=13,
     ylim=c(0, .5), fun='event',
     xlab="Years from HIV infection", ylab="Probability")
title("SI appearance")    
par(oldpar)

## ----fig-T9, fig.cap="Reproduction of Figure T9 showing non-parametric cumulative incidence functions for AIDS (left) and SI appearance (right), for CCR5 wild-type (WW) and mutant (WM)"----
# re-create figure T9: curves by CCR type
aj2 <- survfit(Surv(time, event) ~ ccr5, data = aidssi)
oldpar <- par(mfrow=c(1,2))
plot(aj2[,"AIDS"], xmax=13, col=1:2, lwd=2, ylim=c(0, .5),
     xlab="Years from HIV infection", ylab="Probability of AIDS")
text(c(10, 10), c(.35, .07), c("WW", "WM"))

plot(aj2[,"SI"], xmax=13, col=1:2, lwd=2, ylim=c(0, .5), 
     xlab="Years from HIV infection", ylab="Probability of SI")
text(c(8, 8), c(.34, .18), c("WW", "WM"))
par(oldpar)

## ----fig-T13, fig.height=3, fig.width=4, fig.cap="Reproduce Figure T13"----
oldpar <- par(mar=c(0,0,0,0))
states <- c("Transplant", "Platelet recovery", 
            "Relapse or death")
tmat <- matrix(0, 3,3, dimnames=list(states, states))
tmat[1,2] <- tmat[1,3] <- tmat[2,3] <- 1 # arrows
statefig(cbind((1:3)/4, c(1,3,1)/4), tmat)
text(c(.3, .5, .7), c(.5, .3, .5), c(1169, 458, 383))
par(oldpar)

## ----tableT2--------------------------------------------------------
table(ebmt3$dissub)
table(ebmt3$drmatch)
table(ebmt3$tcd)
table(ebmt3$age)

## ----data1----------------------------------------------------------
temp <- subset(ebmt3, select = -c(prtime, prstat, rfstime, rfsstat))
edata <- tmerge(temp, ebmt3, id, 
                rstat = event(rfstime, rfsstat),
                pstat = event(prtime, prstat),
                priorpr = tdc(prtime))
print(edata[15:20,-(3:5)])

# Check that no one had recovery and death on the same day
with(edata, table(rstat, pstat))

# Create the factor outcome
edata$event <- with(edata, factor(pstat + 2*rstat, 0:2,
                           labels = c("censor", "PR", "RelDeath")))
levels(edata$drmatch) <- c("Match", "Mismatch")
survcheck(Surv(tstart, tstop, event) ~1, edata, id=id)

## ----fig-data1b-----------------------------------------------------
surv1 <- survfit(Surv(tstart, tstop, event) ~ 1, edata, id=id)
surv1$transitions   # matches the Frequencies on page C5
plot(surv1, col=1:2, xscale=365.25, lwd=2, 
     xlab="Years since transplant", ylab="Fraction in state")
legend(1000, .2, c("Platelet recovery", "Death or Relapse"), 
       lty=1, col=1:2, lwd=2, bty='n')

## ----efit1----------------------------------------------------------
efit1 <- coxph(Surv(tstart, tstop, event) ~ dissub + age + drmatch + tcd,
               id=id, data=edata, ties='breslow')
print(efit1, digits=2)

## ----fig-T14, fig.cap="Reproduce Figure T14"------------------------
# a data set containing the "reference" categories 
rdata <- data.frame(dissub="AML", age="<=20", drmatch="Match", tcd="No TCD")
esurv1 <- survfit(efit1, newdata=rdata)
plot(esurv1, cumhaz=TRUE, lty=1:3, xscale=365.25, xmax=7*365.35,
     xlab="Years since transplant", ylab="Cumulative hazard")
legend(365, .8, c("Transplant to platelet recovery (1:2)",
                "Transplant to death (1:3)",
                "Platelet recovery to death (2:3)"), lty=1:3, bty='n')

## ----efit2----------------------------------------------------------
efit2 <- coxph(list(Surv(tstart, tstop, event) ~ dissub + age + drmatch + tcd,
                    0:state("RelDeath") ~ 1 / shared),
                    id=id, data=edata, ties='breslow')
print(coef(efit2, type='matrix'), digits=2)

## ----efit3----------------------------------------------------------
prtime <- ifelse(edata$priorpr==1, edata$tstart, 0)/365.25
efit3 <-  coxph(list(Surv(tstart, tstop, event) ~ dissub + age + drmatch + tcd,
                    0:state("RelDeath") ~ 1/ shared,
                    "PR":"RelDeath" ~ prtime), 
                    id=id, data=edata, ties='breslow')
print(coef(efit3, type='matrix'), digits=2)

## -------------------------------------------------------------------
edummy <- expand.grid(age="<=20", dissub="AML", drmatch="Mismatch",
                      tcd=c("No TCD", "TCD"), priorpr=1)
ecurve2 <- survfit(efit2, newdata= edummy)
plot(ecurve2, col=c(1,1,2,2,3,3), lty=1:2, lwd=2, xscale=365.25,
     noplot=NULL, 
     xlab="Years since transplant", ylab="Predicted probabilities")
legend(700, .9, c("Currently alive in remission, no PR", "Currently in PR",
               "Relapse or death"), col=1:3, lwd=2, bty='n')
text(700, .95, "Solid= No TCD, dashed = TCD", adj=0)

## ----fourstate, fig.height=3, fig.width=5---------------------------
oldpar <- par(mar=c(0,0,0,0))
state4 <- c("Transplant", "Platelet recovery", "Relapse or death (1)",
            "Relapse or death (2)")
cmat <- matrix(0, 4, 4, dimnames = list(state4, state4))
cmat[1,2] <- cmat[1,3] <- cmat[2,4] <- 1
statefig(c(1,2,1), cmat)
par(oldpar)

## ----fourstate2-----------------------------------------------------
etemp <- as.numeric(edata$event)
etemp <- ifelse(etemp==3 & edata$priorpr==1, 4, etemp)
edata$event4 <- factor(etemp, 1:4, c("censor", "PR", "RelDeath1", 
                                     "RelDeath2")) 
survcheck(Surv(tstart, tstop, event4) ~ 1, edata, id=id)

efit4 <- coxph(list(Surv(tstart, tstop, event4) ~ dissub + age + drmatch + tcd,
                    1:3 + 2:4 ~ 1/ shared),
                    id=id, data=edata, ties='breslow')
efit4$cmap
# some of the coefficient names change, but not the values
all.equal(coef(efit4), coef(efit2), check.attributes= FALSE)

## ----fig-T15, fig.cap="Figure T15"----------------------------------
edummy <- expand.grid(dissub="AML", age= "<=20", drmatch="Match",
                      tcd=c("No TCD", "TCD"), priorpr=1)
ecurve4 <- survfit(efit4, newdata=edummy)

oldpar <- par(mfrow=c(1,2), mar=c(4.1, 3.1, 3.1, .1))

pstack(ecurve4[,1,c(2,4,3,1)],
       xscale=365.25, ylim=c(0,1),
       xlab="Years since transplant", ylab="Predicted probabilities")
text(rep(4*365, 4), c(.35, .54, .66, .9), cex=.7, 
     c("Alive in remission, PR", "Relapse or death after PR",
       "Relapse or death without PR", "Alive in remission, no PR"))
title("No TCD")

pstack(ecurve4[,2,c(2,4,3,1)],
       xscale=365.25, ylim=c(0,1),
       xlab="Years since transplant", ylab="Predicted probabilities")
text(rep(4*365, 4), c(.35, .65, .8, .95), cex=.7, 
     c("Alive in remission, PR", "Relapse or death after PR",
       "Relapse or death without PR", "Alive in remission, no PR"))
title("TCD")
par(oldpar)

