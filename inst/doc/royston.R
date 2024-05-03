## ----setup, include = FALSE-----------------------------------------
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

## ----table1, echo= TRUE---------------------------------------------
rotterdam2 <- subset(rotterdam, nodes>0)
table(rotterdam2$size)
table(rotterdam2$meno)
table(rotterdam2$hormon)
round(c(mean(rotterdam2$age), sd(rotterdam2$age)), 1)
round(c(mean(rotterdam2$nodes),sd(rotterdam2$nodes)), 1)
round(c(mean(rotterdam2$pgr),  sd(rotterdam2$pgr)), 1)
round(c(mean(rotterdam2$er), sd(rotterdam2$er)),1)

## ----table2, echo=TRUE----------------------------------------------
y7 <- round(7*365.25)  # 7 years or 84 months
r7 <- rotterdam2
r7$recur <- ifelse(r7$rtime > y7, 0, r7$recur)
r7$rtime <- pmin(r7$rtime, y7)
r7$death <- ifelse(r7$dtime > y7, 0, r7$death)
r7$dtime <- pmin(r7$dtime, y7)

r7$rfstime <- with(r7, ifelse(recur==1, rtime, dtime))  # time to recur or death
r7$rfs <- with(r7, pmax(death, recur)) 
table(r7$rfs)
agefun <- function(x) cbind((x/100)^3, (x/100)^3 * log(x/100))
cfit <- coxph(Surv(rfstime, rfs) ~ agefun(age) + meno + size + 
                    I(1/sqrt(nodes)) + I(er/1000) + hormon, 
              data= r7, ties="breslow")
cbind("cfit"= round(coef(cfit),3), 
      "paper"= c(1.07, 9.13, 0.46, 0.23, 0.31, -1.74, -0.34, -0.35))

## ----termplot-------------------------------------------------------
termplot(cfit, term=1, ylab="Estimated age effect", col.term=1, col.se=1,
         se=TRUE)

## ----gbsg1, echo=TRUE-----------------------------------------------
table(cut(gbsg$size, c(0, 20, 50, 150), c("<=20", "20-50", ">50")))
table(gbsg$meno)
table(gbsg$hormon)
round(c(mean(gbsg$age), sd(gbsg$age)), 1)
round(c(mean(gbsg$nodes),sd(gbsg$nodes)), 1)
round(c(mean(gbsg$pgr), sd(gbsg$pgr)), 1)
round(c(mean(gbsg$er),  sd(gbsg$er)), 1)

## ----gbsg2----------------------------------------------------------
gbsg2 <- gbsg
gbsg2$size <- cut(gbsg$size, c(0, 20, 50, 150), c("<=20", "20-50", ">50"))
gbsg2$rfs <- gbsg2$status

gbsg2$PI <- predict(cfit, newdata=gbsg2)

## ----hist1, echo=FALSE----------------------------------------------
PI1 <- predict(cfit) - mean(predict(cfit))
oldpar <- par(mfrow=c(2,1), mar=c(6,6,1,1))
hist(PI1, breaks=seq(-1.6, 1.6, by=.1), xlab="PI in derivation data", main=NULL)
abline(v = quantile(PI1, c(.16, .5 , .84)), col=2)

PI2 <- gbsg2$PI - mean(predict(cfit))
hist(PI2, breaks=seq(-1.6, 1.6, by=.1), xlab="PI in validation data", main=NULL)
abline(v = quantile(PI2, c(.16, .5 , .84)), col=2)
par(oldpar)

## ----fig2-----------------------------------------------------------
grp1 <- cut(PI1, quantile(PI1, c(0, .16, .5 , .84, 1)), include.lowest=TRUE)
km1 <- survfit(Surv(rfstime, rfs) ~ grp1, r7)
plot(km1, xscale=365.25, xlab="Years since enrollment", 
     ylab="Relapse free survival")
grp2 <- cut(PI2, quantile(PI2, c(0, .16, .5 , .84, 1)), include.lowest=TRUE) 
km2 <- survfit(Surv(rfstime, rfs) ~ grp2, gbsg2)
lines(km2, col=2)

