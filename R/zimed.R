#' The effect decomposition of mediation analysis for zero-inflated models in microbiome data
#'
#' @param M_mat a data frame or matrix of the microbiome count data. Rows of the matrix represent samples and columns are the taxa
#' @param Treat a binary vector of group indicators
#' @param Outcome a continuous outcome
#' @param method methods for testing the H0 hypothesis, including joint, HDMT or DACT
#' @param ci.method methods for calculating confidence interval, including bootstrap or delta
#' @return NIE the estimated natural indirect effect
#' @return NIE.p the calculated p value for NIE
#' @return NIE.ci the calculated confidence interval for NIE
#' @return NIEA the estimated natural indirect effect through changes of abundance
#' @return NIEA.p the calculated p value for NIEA
#' @return NIEA.ci the calculated confidence interval for NIEA
#' @return NIEP the estimated natural indirect effect through changes of absence/presence
#' @return NIEP.p the calculated p value for NIEP
#' @return NIEP.ci the calculated confidence interval for NIEP

#' @examples
#' Treat <- exposure
#' Outcome <- Y_mat
#' zimed.result <- zimed(M_mat,Treat,Outcome,method="joint",ci.method="delta")
#' @export
#' @import pscl



zimed <- function(M_mat,Treat,Outcome,method,ci.method){
  beta.p <- c()
  for(j in 1:ncol(M_mat)){
    Y_data <- as.data.frame(cbind(Treat,M_mat[,j],Outcome))
    colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
    lm.res <-  try(lm(Outcom ~ ., data = Y_data),silent=TRUE)
    beta.p[j] <-summary(lm.res)$coefficients[3,4]
    }

  wald.p <- c()
  for(j in 1:ncol(M_mat)){
  Y_data <- as.data.frame(cbind(Treat,M_mat[,j],Outcome))
  colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
  nb.res <-  try(glm.nb(Mediator ~ Treat, data = Y_data),silent=TRUE)
  scr.t <- try(overdisp_scoretest(OTU_table=Y_data$Mediator,group=Y_data$Treat), silent = TRUE)
  if(sum(Y_data$Mediator==0)==0){
    if(inherits(b,"error")){wald.p[j]<-NA}
    else{
      wald.p[j] <- summary(nb.res)$coefficients[2,4]
    }
  }
  else{
    scr.t <- try(overdisp_scoretest(OTU_table=Y_data$Mediator,group=Y_data$Treat), silent = TRUE)
    if(scr.t <= 0.05){
      zinb.res <-  try(zeroinfl(Mediator ~ Treat | Treat, data = Y_data,dist = c("negbin"),link = c("logit")), silent = TRUE)
      if(inherits(d, "try-error")){
        wald.p[j] <- NA
      }
      else{
        alpha_hat<- summary(zinb.res)$coefficients$count[2,1]
        gamma_hat <- summary(zinb.res)$coefficients$zero[2,1]

        alpha_cov<- summary(zinb.res)$coefficients$count[2,2]
        gamma_cov <- summary(zinb.res)$coefficients$zero[2,2]
        num <- nrow(summary(zinb.res)$coefficients$count)-3

        cov<- matrix(c(vcov(zinb.res)[2,2],vcov(zinb.res)[2,4+num],vcov(zinb.res)[4+num,2],vcov(zinb.res)[4+num,4+num]),2,2)
        wald.t <- try(t(matrix(c(alpha_hat,gamma_hat)))%*%ginv(cov)%*%matrix(c(alpha_hat,gamma_hat)), silent = TRUE)
        if(inherits(wald.t, "try-error")){
          wald.p[j] <- NA
        }
        else{
          wald.p[j] <- 1-pchisq(as.numeric(wald.t),df=2)
        }

      }
    }
    else{
      zip.res <-  try(zeroinfl(Mediator ~ Treat | Treat, offset=tij_mat,data = Y_data,dist = c("poisson"),link = c("logit")), silent = TRUE)

      print("ZIP")
      alpha_hat<- summary(zip.res)$coefficients$count[2,1]
      gamma_hat <- summary(zip.res)$coefficients$zero[2,1]
      alpha_cov<- summary(zip.res)$coefficients$count[2,2]
      gamma_cov <- summary(zip.res)$coefficients$zero[2,2]
      num <- nrow(summary(d)$coefficients$count)-3


      cov<- matrix(c(vcov(zip.res)[2,2],vcov(zip.res)[2,4+num],vcov(zip.res)[4+num,2],vcov(zip.res)[4+num,4+num]),2,2)
      wald.t <- try(t(matrix(c(alpha_hat,gamma_hat)))%*%ginv(cov)%*%matrix(c(alpha_hat,gamma_hat)), silent = TRUE)
      wald.p[j] <- 1-pchisq(as.numeric(wald.t),df=2)
    }

  }
  }
  if(method=="joint"){
  bind.p <- rbind(p.adjust(beta.p,"fdr"),p.adjust(wald.p,"fdr"))
  final.p <- apply(bind.p, 2, max)
  }
  if(method=="HDMT"){
    bind.p <- cbind(wald.p,beta.p)
    nullprop <- try(HDMT::null_estimation(bind.p,lambda=0.5),silent = TRUE)
    final.p  <- HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10, nullprop$alpha1,nullprop$alpha2,input_pvalues,exact=0)
  }
  if(method=="DACT"){
    bind.p <- cbind(wald.p,beta.p)
    final.p <- try(DACT(bind.p[,1],bind.p[,2],"JC"),silent = TRUE)
  }

  id.test <- which(final.p<0.05)

  for(j in 1:length(id.test)){
    Y_data <- as.data.frame(cbind(Treat,M_mat[,id.test[j]],Outcome))
    colnames(Y_data) <- c("Treat", "Mediator", "Outcom")
    zinb.res <-  try(zeroinfl(Mediator ~ Treat| Treat, data = Y_data,dist = c("negbin"),link = c("logit")), silent = TRUE)
    lm.res <-  try(lm(Outcom ~ ., data = Y_data),silent=TRUE)
    bet.val <-summary(lm.res)$coefficients[3,1]
    alpha.val<- summary(zinb.res)$coefficients$count[1:2,1]
    gamma.val <- summary(zinb.res)$coefficients$zero[,1]

    beta.sd <- summary(lm.res)$coefficients[3,2]

    beta.p <- summary(lm.res)$coefficients[3,4]
    alpha.p<- summary(zinb.res)$coefficients$count[2,4]
    gamma.p <- summary(zinb.res)$coefficients$zero[2,4]

    nie.p <- final.p[id.test[[j]]]
    niea.p <- max(beta.p,alpha.p)
    niep.p <- max(beta.p,gamma.p)



    niea <- bet.val*(1/(1+exp(gamma.val[1])))*
            (exp(alpha.val[1]+alpha.val[2])-
             exp(alpha.val[1]))

    niep <- bet.val*(exp(alpha.val[1]+alpha.val[2]))*
            (1/(1+exp(gamma.val[1]+gamma.val[2]))-
             1/(1+exp(gamma.val[1])))

    nie <- niea+niep
    if(ci.method=="bootstrap"){
    nie.perm=niea.perm=niep.perm=c()
    for(jk in 1:1000){
      set.seed(jk)
      Y_data.perm <-Y_data[sample(1:nrow(Y_data),replace = TRUE),]

      zinb.perm <-  try(zeroinfl(Mediator ~ Treat| Treat, data = Y_data.perm,dist = c("negbin"),link = c("logit")), silent = TRUE)
      lm.perm <- summary(glm(Outcom ~ Mediator+Treat, family = "gaussian", data = Y_data.perm))$coefficients

      alpha.val<- summary(zinb.perm)$coefficients$count[1:2,1]
      gamma.val <- summary(zinb.perm)$coefficients$zero[,1]
      bet.val <- lm.perm[2, 1]
      niea.perm[jk] <- bet.val*(1/(1+exp(gamma.val[1])))*
                       (exp(alpha.val[1]+alpha.val[2])-
                        exp(alpha.val[1]))
      niep.perm[jk] <- bet.val*(exp(alpha.val[1]+alpha.val[2]))*
                       (1/(1+exp(gamma.val[1]+gamma.val[2]))-
                        1/(1+exp(gamma.val[1])))

      nie.perm[jk] <- niea.perm[jk]+niep.perm[jk]

    }

    nie.ci <- quantile(nie.perm[!is.na(nie.perm)], probs = c(0.025, 0.975))
    niea.ci  <- quantile(niea.perm[!is.na(nie1.perm)], probs = c(0.025, 0.975))
    niep.ci <- quantile(niep.perm[!is.na(nie2.perm)], probs = c(0.025, 0.975))
    }

    if(ci.method=="delta"){

      var.mat <- as.matrix(bdiag(beta.sd^2,vcov(zinb.res)))

      grad_beta2 <-  (1/(1+exp(gamma.val[1]+gamma.val[2])))*
                     exp(alpha.val[1]+alpha.val[2])-
                     (1/(1+exp(gamma.val[1])))*
                     exp(alpha.val[1])

      grad_alpha0 <- bet.val*((1/(1+exp(gamma.val[1]+gamma.val[2])))*
                                (exp(alpha.val[1]+alpha.val[2]))-
                                (1/(1+exp(gamma.val[1])))*
                                (exp(alpha.val[1])))


      grad_alpha1 <- bet.val*(1/(1+exp(gamma.val[1]+gamma.val[2])))*
                     (exp(alpha.val[1]+alpha.val[2]))



      grad_gamma0<- bet.val*(exp(gamma.val[1]+gamma.val[2])*
                            (-exp(gamma.val[1]+gamma.val[2])/(1+exp(gamma.val[1]+gamma.val[2]))^2)-
                             exp(alpha.val[1])*
                             (-exp(gamma.val[1])/(1+exp(gamma.val[1]))^2))

      grad_gamma1 <-  bet.val*exp(alpha.val[1]+alpha.val[2])*
                      (-exp(gamma.val[1]+gamma.val[2])/(1+exp(gamma.val[1]+gamma.val[2]))^2)


      nie.vec <- matrix(c(grad_beta2,grad_alpha0,grad_alpha1,grad_gamma0,grad_gamma1))

      nie.var <- t(nie.vec)%*%var.mat%*%nie.vec
      nie.ci <- c(nie-1.96*sqrt(nie.var),nie+1.96*sqrt(nie.var))




      grad_beta2 <-  (1/(1+exp(gamma.val[1]+gamma.val[2])))*
                     (exp(alpha.val[1]+alpha.val[2])-
                     exp(alpha.val[1]))

      grad_alpha0 <- bet.val*(1/(1+exp(gamma.val[1]+gamma.val[2])))*
                     (exp(alpha.val[1]+alpha.val[2])-
                     exp(alpha.val[1]))


      grad_alpha1 <- bet.val*(1/(1+exp(gamma.val[1]+gamma.val[2])))*
                     (exp(alpha.val[1]+alpha.val[2]))



      grad_gamma0<- bet.val* (exp(alpha.val[1]+alpha.val[2])-
                                exp(alpha.val[1]))*
                             (-exp(gamma.val[1]+gamma.val[2])/(1+exp(gamma.val[1]+gamma.val[2]))^2)


      grad_gamma1 <- grad_gamma0

      niea.vec <- matrix(c(grad_beta2,grad_alpha0,grad_alpha1,grad_gamma0,grad_gamma1))

      niea.var <- t(niea.vec)%*%var.mat%*%niea.vec
      niea.ci <- c(niea-1.96*sqrt(niea.var),niea+1.96*sqrt(niea.var))




      grad_beta2 <- exp(alpha.val[1]+alpha.val[2])*
                     (1/(1+exp(gamma.val[1]+gamma.val[2]))-1/(1+exp(gamma.val[1])))


      grad_alpha0 <- bet.val*exp(alpha.val[1])*(1/(1+exp(gamma.val[1]+gamma.val[2]))-1/(1+exp(gamma.val[1])))


      grad_alpha1 <- grad_alpha0



      grad_gamma0<- bet.val*exp(alpha.val[1]+alpha.val[2])*
                    (-exp(gamma.val[1]+gamma.val[2])/((1+exp(gamma.val[1]+gamma.val[2]))^2)+
                     exp(gamma.val[1])/((1+exp(gamma.val[1]))^2))

      grad_gamma1 <- bet.val*exp(alpha.val[1]+alpha.val[2])*
                     (-exp(gamma.val[1]+gamma.val[2])/(1+exp(gamma.val[1]+gamma.val[2]))^2)


      niep.vec <- matrix(c(grad_beta2,grad_alpha0,grad_alpha1,grad_gamma0,grad_gamma1))
      niep.var <- t(niep.vec)%*%var.mat%*%niep.vec
      niep.ci <- c(niep-1.96*sqrt(niep.var),niep+1.96*sqrt(niep.var))
    }


  }
  results <- list(NIE=nie,NIE.p =nie.p, NIE.ci=nie.ci,NIEA=niea,NIEA.p =niea.p,NIEA.ci=niea.ci,NIEP=niep,NIEP.p =niep.p,NIEP.ci=niep.ci,)
  return(results)
}
