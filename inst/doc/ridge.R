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

## ----rdata----------------------------------------------------------
set.seed(1954)  # force reproducability
library(survival)

n <- nrow(lung)
snp <- matrix(rbinom(20*n, 1, p=.1), nrow=n)
snpdata <- cbind(lung, data.frame(snp))
dim(snpdata)

## ----pass1----------------------------------------------------------
cfit1 <- coxph(Surv(time, status) ~ age + sex + ridge(X1, X2, X3, X4, X5, X6,
                     X7, X8, X9, X10, X11, X12, X13, X14, X15, X16, X17, X18, 
		     X19, X20, theta=.1), data=snpdata)

## ----pass2----------------------------------------------------------
xlist <-  paste0("X", 1:20)
myform <- paste("Surv(time, status) ~ age + sex + ridge(",
                paste(xlist, collapse= ", "), ", theta=.1)")
cfit2 <- coxph(formula(myform), data=snpdata)

all.equal(cfit1$loglik, cfit2$loglik)

## ----pass3----------------------------------------------------------
cfit3 <- coxph(Surv(time, status) ~ age + sex + ridge(snp, theta=.1),
      	 data = snpdata)

## ----pass4, error=TRUE----------------------------------------------
newsnp <- matrix(rbinom(20*4, 1, p=.12), nrow=4)
prdata <- data.frame(age= c(50, 65, 48, 70), sex= c(1, 1, 2,1),
                     newsnp)
predict(cfit3, newdata=prdata)

## ----pass5----------------------------------------------------------
prdata <- data.frame(age= c(50, 65, 48, 70), sex= c(1, 1, 2,1),
                     newsnp)
snpsave <- snp
snp <- newsnp
age <- prdata$age
sex <- prdata$sex

new <- predict(cfit3)
length(new)

snp <- snpsave  # restore the status quo
rm(age, sex)

## ----pass6----------------------------------------------------------
prmat <- cbind(age= c(50, 65, 48, 70), sex=c(1,1,2,1), newsnp)
drop(prmat %*% coef(cfit3)) # simplify results to vector
#alternate (center)
drop(prmat %*% coef(cfit3)) - sum(coef(cfit3)* cfit3$mean)

## ----pass7----------------------------------------------------------
# new data type
ridgemat <- function(x) {
    class(x) <- c("ridgemat", class(x))
    x
}
as.data.frame.ridgemat <- function(x, ...) 
    as.data.frame.model.matrix(as.matrix(x), ...)

snpdata2 <- cbind(lung, snp= ridgemat(snp))
names(snpdata2)

cfit4 <- coxph(Surv(time, status) ~ age + sex + ridge(snp, theta=.1),
      	     data= snpdata2)
     

prdata <- data.frame(age= c(50, 65, 48, 70), sex= c(1, 1, 2,1),
                       snp= ridgemat(newsnp))
predict(cfit4, newdata=prdata)

