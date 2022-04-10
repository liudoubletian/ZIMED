overdisp_scoretest <- function(OTU_table, group){
  library(pscl)
  library(MASS)
  
  
  nsamples = length(OTU_table)
  ntest_pool = 1
 
  df = data.frame(OTU_table, group)
  colnames(df) =c("otu", "group")
    if (sum(OTU_table==0) >= 1) {
      tryCatch({
        m1 <- zeroinfl(otu ~ group  | group ,data=df)
        para_logPoi= summary(m1)$coefficients[[1]][,1]
        para_zeroinf= summary(m1)$coefficients[[2]][,1]
      }, error = function(e){})
    }else{
      tryCatch({
        m1 <- glm(otu ~ group, family = poisson)
        para_logPoi = summary(m1)$coefficients[,1]
        para_zeroinf = c(0,0)
      }, error = function(e){})
    }
 
  
  # MLE of mu
  B = rbind(rep(1, nsamples), group)
  mu = t(matrix(para_logPoi)) %*% B
  mu = exp(mu)
  # MLE of pi
  A = rbind(rep(1, nsamples), group)
  expo = exp(t(matrix(para_zeroinf)) %*% A)
  pi = expo / (expo + 1)
  # expo can be Inf in some cases, cause pi be NA
  pi[is.na(pi)] = 1
  # adjust for OTUs without 0 count
  #pi[1:num_nonzero_OTU,] = 0
  
  # Compute score function for each OTU
  index_0 = OTU_table
  index_0[index_0 == 0] <- 0.1
  index_0[index_0 >= 1] <- 0
  index_0 = index_0*10
  
  U = ((OTU_table-mu)^2-OTU_table) - c(index_0)*mu^2*pi/(pi+(1-pi)*exp(-mu) + 1e-300)
  U = 0.5 * sum(U)
  # compute inverse of Fisher Information matrix
  inv_I = mu^2*(2*(1-pi)-mu^2*pi*(1-pi/(pi+(1-pi)*exp(-mu) + 1e-300)))
  inv_I = 0.25*sum(inv_I)
  
  # test statistic: sqrt of Score statistic
  test_stat = U*sqrt(inv_I)
  # calculate p-values
  overdisp_p = 1-pnorm(test_stat)
  
  
  return(overdisp_p)
}
