loglik.cond = function(X, y, beta, family) {
    K = dim(beta)[2]
    link = cbind(1, X) %*% beta
    yrep = repmat.cond(y, 1, K)
    if (family == "gaussian")
        return(apply((yrep - link)^2, 2, sum))
    if (family == "poisson")
        return(apply(exp(link) - yrep * link, 2, sum))
    if (family == "binomial")
        return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
}

repmat.cond = function(X, m, n) {
    ## R equivalent of repmat (matlab)
    X = as.matrix(X)
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}

getdf.cond = function(coef.beta) {
    apply(abs(coef.beta) > 1e-10, 2, sum)
}


margcoef.cond <- function(x, y, exposure, condind = NULL, family, null.model = FALSE, iterind) {
  n = dim(x)[1]
  p = dim(x)[2]
  ones = rep(1, n)
  candind = setdiff(1:p, condind)
  if (iterind == 0) {
    if (family == "cox"){margcoef = abs(cor(x, y[, 1]))} else {margcoef = abs(cor(x, y))}
  } else {
    if (null.model == TRUE) {
      if (is.null(condind) == TRUE) {
        x = x[sample(1:n), ]
      }
      if (is.null(condind) == FALSE) {
        x[, candind] = x[sample(1:n), candind]
      }
    }
    margcoef = sapply(candind, mg.cond, x, y, ones, family, condind,exposure)
  }
  return(margcoef)
}

mg.cond <- function(index, x = x, y = y, ones = ones, family = family, condind = condind, exposure) {
  margfit = switch(family, gaussian = summary.glm(glm.fit(cbind(ones, x[, index], x[, condind],exposure), y, family = gaussian()))$coefficients[2,4],
                   binomial = summary.glm(glm.fit(cbind(ones, x[, index], x[, condind],exposure), y, family = binomial()))$coefficients[2,4], poisson = summary.glm(glm.fit(cbind(ones,
                                                                                                                                            x[, index], x[, condind],exposure), y, family = poisson()))$coefficients[2,4], cox = summary(coxph(y ~ cbind(x[, index], x[,
                                                                                                                                                                                                                                            condind],exposure)))$coefficients[1,5])
}


obtain.ix0.cond <- function(x, y, exposure, s1, s2, family, nsis, iter, varISIS, perm, q, greedy, greedy.size, iterind) {
  if (iter == FALSE) {
    margcoef = margcoef.cond(x, y, exposure=exposure, family = family, null.model = FALSE, iterind = iterind)
    rankcoef = sort(margcoef, decreasing = FALSE, index.return = TRUE)
    ix0 = rankcoef$ix[1:nsis]
  } else {
    if (varISIS == "vanilla") {
      margcoef = margcoef.cond(x, y, exposure=exposure, family = family, null.model = FALSE, iterind = iterind)
      rankcoef = sort(margcoef, decreasing = FALSE, index.return = TRUE)
      if (perm == FALSE)
        ix0 = rankcoef$ix[1:floor((2/3) * nsis)] else {
          count = 0
          repeat {
            count = count + 1
            randcoef = margcoef.cond(x, y, exposure=exposure, family = family, null.model = TRUE, iterind = iterind)
            if (length(which(margcoef >= quantile(randcoef, q))) > 0)
              break
            if(count > 10)
              break
          }
          if (greedy == FALSE) {
            if (length(which(margcoef >= quantile(randcoef, q))) >= 2) {
              length1 = length(which(margcoef >= quantile(randcoef, q)))
              above.thresh = rankcoef$ix[1:length1]
              ix0 = rankcoef$ix[1:floor((2/3) * nsis)]
              ix0 = sort(intersect(ix0, above.thresh))
            } else ix0 = rankcoef$ix[1:2]
          } else {
            if (greedy.size == 1)
              ix0 = rankcoef$ix[1:2] else ix0 = rankcoef$ix[1:greedy.size]
          }
        }
    } else {
      if(family == 'cox'){
        margcoef1 = margcoef.cond(x[s1, ], y[s1,], exposure=exposure[s1,], family = family, null.model = FALSE, iterind = iterind)
        margcoef2 = margcoef.cond(x[s2, ], y[s2,], exposure=exposure[s2,], family = family, null.model = FALSE, iterind = iterind)
      } else{
        margcoef1 = margcoef.cond(x[s1, ], y[s1], exposure=exposure[s1], family = family, null.model = FALSE, iterind = iterind)
        margcoef2 = margcoef.cond(x[s2, ], y[s2], exposure=exposure[s2], family = family, null.model = FALSE, iterind = iterind)

      }

      rankcoef1 = sort(margcoef1, decreasing = FALSE, index.return = TRUE)
      rankcoef2 = sort(margcoef2, decreasing = FALSE, index.return = TRUE)
      if (perm == FALSE) {
        if (varISIS == "aggr") {
          ix01 = rankcoef1$ix[1:floor((2/3) * nsis)]
          ix02 = rankcoef2$ix[1:floor((2/3) * nsis)]
          ix0 = sort(intersect(ix01, ix02))
          if (length(ix0) <= 1)
            ix0 = int.size.k.cond(rankcoef1$ix, rankcoef2$ix, 2)
        }
        if (varISIS == "cons") {
          iensure = intensure.cond(floor((2/3) * nsis), l1 = rankcoef1$ix, l2 = rankcoef2$ix, k = floor((2/3) * nsis))
          ix01 = rankcoef1$ix[1:iensure]
          ix02 = rankcoef2$ix[1:iensure]
          ix0 = sort(intersect(ix01, ix02))
        }
      } else {
        count = 0
        repeat {
          count  = count + 1
          randcoef1 = margcoef.cond(x[s1, ], y[s1], exposure=exposure[s1], family = family, null.model = TRUE, iterind = iterind)
          randcoef2 = margcoef.cond(x[s2, ], y[s2], exposure=exposure[s2], family = family, null.model = TRUE, iterind = iterind)
          if (length(which(margcoef1 >= quantile(randcoef1, q))) > 0 && length(which(margcoef2 >= quantile(randcoef2,
                                                                                                           q))) > 0)
            break
          if(count > 10) break
        }
        if (greedy == FALSE) {
          length1 = length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 = length(which(margcoef2 >= quantile(randcoef2, q)))
          above.thresh.1 = rankcoef1$ix[1:length1]
          above.thresh.2 = rankcoef2$ix[1:length2]
          ix01 = rankcoef1$ix[1:floor((2/3) * nsis)]
          ix02 = rankcoef2$ix[1:floor((2/3) * nsis)]
          ix01 = sort(intersect(ix01, above.thresh.1))
          ix02 = sort(intersect(ix02, above.thresh.2))
          ix0 = sort(intersect(ix01, ix02))
          if (length(ix0) <= 1)
            ix0 = int.size.k.cond(rankcoef1$ix, rankcoef2$ix, 2)
        } else {
          if (greedy.size == 1)
            ix0 = int.size.k.cond(rankcoef1$ix, rankcoef2$ix, 2) else ix0 = int.size.k.cond(rankcoef1$ix, rankcoef2$ix, greedy.size)
        }
      }
    }
  }
  return(ix0)
}

obtain.newix.cond <- function(x, y, exposure, ix1, candind, s1, s2, family, pleft, varISIS, perm, q, greedy, greedy.size, iterind) {
  if (varISIS == "vanilla") {
    margcoef = margcoef.cond(x, y, exposure=exposure, ix1, family = family, null.model = FALSE, iterind = iterind)
    rankcoef = sort(margcoef, decreasing = FALSE, index.return = TRUE)
    if (perm == FALSE) {
      if (pleft > 0)
        newix = candind[rankcoef$ix[1:pleft]] else newix = NULL
    } else {
      randcoef = margcoef.cond(x, y, exposure=exposure, ix1, family = family, null.model = TRUE, iterind = iterind)
      if (length(which(margcoef >= quantile(randcoef, q))) > 0) {
        if (greedy == FALSE) {
          length1 = length(which(margcoef >= quantile(randcoef, q)))
          above.thresh = candind[rankcoef$ix[1:length1]]
          newix = candind[rankcoef$ix[1:pleft]]
          newix = sort(intersect(newix, above.thresh))
        } else newix = candind[rankcoef$ix[1:greedy.size]]
      } else newix = NULL
    }
  } else {
    margcoef1 = margcoef.cond(x[s1, ], y[s1], exposure=exposure[s1], ix1, family = family, null.model = FALSE, iterind = iterind)
    margcoef2 = margcoef.cond(x[s2, ], y[s2], exposure=exposure[s2], ix1, family = family, null.model = FALSE, iterind = iterind)
    rankcoef1 = sort(margcoef1, decreasing = FALSE, index.return = TRUE)
    rankcoef2 = sort(margcoef2, decreasing = FALSE, index.return = TRUE)
    if (perm == FALSE) {
      if (pleft > 0) {
        if (varISIS == "aggr") {
          newix1 = candind[rankcoef1$ix[1:pleft]]
          newix2 = candind[rankcoef2$ix[1:pleft]]
          newix = sort(intersect(newix1, newix2))
        }
        if (varISIS == "cons") {
          iensure = intensure.cond(pleft, l1 = rankcoef1$ix, l2 = rankcoef2$ix, k = pleft)
          newix1 = candind[rankcoef1$ix[1:iensure]]
          newix2 = candind[rankcoef2$ix[1:iensure]]
          newix = sort(intersect(newix1, newix2))
        }
      } else newix = NULL
    } else {
      randcoef1 = margcoef.cond(x[s1, ], y[s1], exposure=exposure[s1], ix1, family = family, null.model = TRUE, iterind = iterind)
      randcoef2 = margcoef.cond(x[s2, ], y[s2], exposure=exposure[s2], ix1, family = family, null.model = TRUE, iterind = iterind)
      if (length(which(margcoef1 >= quantile(randcoef1, q))) > 0 && length(which(margcoef2 >= quantile(randcoef2,                                                                          q))) > 0) {
        if (greedy == FALSE) {
          length1 = length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 = length(which(margcoef2 >= quantile(randcoef2, q)))
          above.thresh.1 = candind[rankcoef1$ix[1:length1]]
          above.thresh.2 = candind[rankcoef2$ix[1:length2]]
          newix1 = candind[rankcoef1$ix[1:pleft]]
          newix2 = candind[rankcoef2$ix[1:pleft]]
          newix1 = sort(intersect(newix1, above.thresh.1))
          newix2 = sort(intersect(newix2, above.thresh.2))
          newix = sort(intersect(newix1, newix2))
        } else {
          length1 = length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 = length(which(margcoef2 >= quantile(randcoef2, q)))
          newix1 = candind[rankcoef1$ix[1:length1]]
          newix2 = candind[rankcoef2$ix[1:length2]]
          iensure = intensure.cond(greedy.size, l1 = newix1, l2 = newix2, k = greedy.size)
          if(is.null(iensure)) newix = NULL
          else newix = sort(intersect(newix1[1:iensure], newix2[1:iensure]))
        }
      } else newix = NULL
    }
  }
  return(newix)
}


intensure.cond <- function(i, l1, l2, k) {
  for(j in i:length(l1)){
    if (length(intersect(l1[1:j], l2[1:j])) >= k)
      return(j)
  }
 # if (length(intersect(l1[1:i], l2[1:i])) >= k)
#    return(i) else return(intensure(i + 1, l1, l2, k))
}

int.size.k.cond <- function(l1, l2, k) {
  iensure = intensure.cond(k, l1 = l1, l2 = l2, k = k)
  ix01 = l1[1:iensure]
  ix02 = l2[1:iensure]
  ix0 = sort(intersect(ix01, ix02))
  return(ix0)
}

calculate.nsis.cond <- function(family, varISIS, n, p) {
  if (varISIS == "aggr")
    nsis = floor(n/log(n)) else {
      if (family == "gaussian") {
        nsis = floor(n/log(n))
      }
      if (family == "binomial") {
        nsis = floor(n/(4 * log(n)))
      }
      if (family == "poisson") {
        nsis = floor(n/(2 * log(n)))
      }
      if (family == "cox") {
        nsis = floor(n/(4 * log(n)))
      }
    }
  if (p < n)
    nsis = p
  return(nsis)
}
SIS.cond <- function(x, y, exposure,family = c("gaussian", "binomial", "poisson", "cox"), penalty = c("SCAD", "MCP", "lasso"),
                concavity.parameter = switch(penalty, SCAD = 3.7, 3), tune = c("bic", "ebic", "aic", "cv"), nfolds = 10,
                type.measure = c("deviance", "class", "auc", "mse", "mae"), gamma.ebic = 1, nsis = NULL, iter = TRUE, iter.max = ifelse(greedy ==
                                                                                                                                          FALSE, 10, floor(nrow(x)/log(nrow(x)))), varISIS = c("vanilla", "aggr", "cons"), perm = FALSE, q = 1,
                greedy = FALSE, greedy.size = 1, seed = NULL, standardize = TRUE) {

  #this.call = match.call()
  # family = match.arg(family)
  # penalty = match.arg(penalty)
  # tune = match.arg(tune)
  # type.measure = match.arg(type.measure)
  # varISIS = match.arg(varISIS)

  if (is.null(x) || is.null(y))
    stop("The data is missing!")
  if (class(concavity.parameter) != "numeric")
    stop("concavity.parameter must be numeric!")
  if (class(nfolds) != "numeric")
    stop("nfolds must be numeric!")
  if (!is.null(seed) &  class(seed) != "numeric")
    stop("seed must be numeric!")

  if (family == "cox" && penalty %in% c("SCAD", "MCP"))
    stop("Cox model currently not implemented with selected penalty")

  if (type.measure %in% c("class", "auc") && family %in% c("gaussian", "poisson", "cox"))
    stop("'class' and 'auc' type measures are only available for logistic regression")

  if (type.measure %in% c("class", "auc", "mse", "mae") && penalty %in% c("SCAD", "MCP"))
    stop("Only 'deviance' is available as type.measure for non-convex penalties")

  fit = switch(family, gaussian = sisglm.cond(x, y, "gaussian", penalty,exposure, concavity.parameter, tune, nfolds, type.measure,
                                         gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize), binomial = sisglm.cond(x,
                                                                                                                                                        y, "binomial", penalty, exposure,concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max,
                                                                                                                                                        varISIS, perm, q, greedy, greedy.size, seed, standardize), poisson = sisglm.cond(x, y,"poisson", penalty,exposure,
                                                                                                                                                                                                                                    concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q,
                                                                                                                                                                                                                                    greedy, greedy.size, seed, standardize), cox = sisglm.cond(x, y, "cox", penalty,exposure, concavity.parameter, tune,
                                                                                                                                                                                                                                                                                          nfolds, type.measure, gamma.ebic, nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed,
                                                                                                                                                                                                                                                                                          standardize))
  #fit$call = this.call
  #class(fit) = c(class(fit), "SIS")
  return(fit)
}

sisglm.cond <- function(x, y, family, penalty, exposure, concavity.parameter, tune, nfolds, type.measure, gamma.ebic, nsis,
                   iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1 = NULL, s2 = NULL, split.tries = 0) {

  storage.mode(x) = "numeric"
  n = dim(x)[1]
  p = dim(x)[2]
  models = vector("list")
  if (is.null(nsis) == TRUE)
    nsis = calculate.nsis.cond(family = family, varISIS = varISIS, n = n, p = p)
  if (is.null(s1) == TRUE) {
    if(!is.null(seed)){
      set.seed(seed)
    }
    split.sample = sample(1:n)
    s1 = split.sample[1:ceiling(n/2)]
    s2 = setdiff(split.sample, s1)
  }
  old.x = x#cbind(x,exposure)
  if (standardize == TRUE) {
    x = scale(x)
  }
  iterind = 0

  if (iter == TRUE) {
    ix0 = sort(obtain.ix0.cond(x = x, y = y, s1 = s1, s2 = s2,exposure=exposure, family = family, nsis = nsis, iter = iter, varISIS = varISIS,
                          perm = perm, q = q, greedy = greedy, greedy.size = greedy.size, iterind = iterind))
    sis.ix0 = ix0
    repeat {
      iterind = iterind + 1
      cat("Iter", iterind, ", screening: ", ix0, "\n")
      if(length(ix0) == 1 & penalty == 'lasso'){
        ix0 = c(ix0, p+1-ix0)
      }
      pen.ind = ix0
      selection.fit = tune.fit.cond(old.x[,ix0,drop = FALSE], y, exposure,family , penalty , concavity.parameter, tune, nfolds , type.measure , gamma.ebic)
      coef.beta = selection.fit$beta
      a0 = selection.fit$a0

      lambda  = selection.fit$lambda
      lambda.ind = selection.fit$lambda.ind
      ix1 = sort(ix0[selection.fit$ix])
      if (length(ix1) == 0) {
        split.tries = split.tries + 1
        split.sample = sample(1:n)
        s1 = split.sample[1:ceiling(n/2)]
        s2 = setdiff(split.sample, s1)
        cat("Sample splitting attempt: ", split.tries, "\n")
        if (split.tries >= 20) {
          cat("No variables remaining after ", split.tries, " sample splitting attempts! \n")
          cat("You can try a more conservative variable screening approach! \n")
        } else return(sisglm.cond(old.x, y, exposure=exposure,family, penalty, concavity.parameter, tune, nfolds, type.measure, gamma.ebic,
                             nsis, iter, iter.max, varISIS, perm, q, greedy, greedy.size, seed, standardize, s1, s2, split.tries))
      }


      cat("Iter", iterind, ", selection: ", ix1, "\n")
      if (length(ix1) >= nsis || iterind >= iter.max) {
        ix0 = ix1
        if (length(ix1) >= nsis)
          cat("Maximum number of variables selected \n")
        if (iterind >= iter.max)
          cat("Maximum number of iterations reached \n")
        break
      }

      models[[iterind]] = ix1
      flag.models = 0
      if (iterind > 1) {
        for (j in 1:(iterind - 1)) {
          if (identical(models[[j]], ix1) == TRUE)
            flag.models = 1
        }
      }
      if (flag.models == 1) {
        cat("Model already selected \n")
        break
      }

      candind = setdiff(1:p, ix1)
      pleft = nsis - length(ix1)
      newix = sort(obtain.newix.cond(x = x, y = y,exposure=exposure, candind = candind, ix1 = ix1, s1 = s1, s2 = s2, family = family,
                                pleft = pleft, varISIS = varISIS, perm = perm, q = q, greedy = greedy, greedy.size = greedy.size,
                                iterind = iterind))
      cat("Iter", iterind, ", conditional-screening: ", newix, "\n")
      ix1 = sort(c(ix1, newix))
      if (setequal(ix1, ix0)) {
        flag.models = 1
      }
      ix0 = ix1
      if(length(ix1) == 0) break
    }  # end repeat

  } else {
    # end if(iter==TRUE)
    ix0 = sort(obtain.ix0.cond(x = x, y = y, exposure=exposure,s1 = s1, s2 = s2, family = family, nsis = nsis, iter = iter, varISIS = varISIS,
                          perm = perm, q = q, greedy = greedy, greedy.size = greedy.size, iterind = iterind))
    sis.ix0 = ix0
    if(length(ix0) == 1 & penalty == 'lasso'){
      ix0 = c(ix0, p+1-ix0)
    }
    pen.ind = ix0
    selection.fit = tune.fit.cond(old.x[,ix0,drop = FALSE], y, exposure,family , penalty , concavity.parameter, tune, nfolds , type.measure , gamma.ebic)
    coef.beta = selection.fit$beta
    a0 = selection.fit$a0
    lambda = selection.fit$lambda
    lambda.ind = selection.fit$lambda.ind
    ix1 = sort(ix0[selection.fit$ix])
  }



  if (family == "cox") {
    if (length(ix1) > 0){
      names(coef.beta) = paste("X", ix1, sep = "")
    }
  }  else {
    coef.beta = c(a0, coef.beta)
    if(length(ix1)>0){
      names(coef.beta) = c("(Intercept)", paste("X", ix1, sep = ""))
    }
  }

  return(list(sis.ix0 = sis.ix0, ix = ix1, coef.est = coef.beta, fit = selection.fit$fit, lambda = lambda, lambda.ind = lambda.ind, ix0 = pen.ind))
}
tune.fit.cond <- function(x, y, exposure,family = c("gaussian", "binomial", "poisson", "cox"), penalty = c("SCAD", "MCP", "lasso"), concavity.parameter = switch(penalty, SCAD = 3.7, 3), tune = c("cv", "aic", "bic", "ebic"), nfolds = 10,
                          type.measure = c("deviance", "class", "auc", "mse", "mae"), gamma.ebic = 1) {

  if(is.null(exposure)){x <- x} else{x <- cbind(x,exposure)}
  if (is.null(x) || is.null(y))
    stop("The data is missing!")

  #this.call = match.call()
  # family = match.arg(family)
  # penalty = match.arg(penalty)
  if (class(concavity.parameter) != "numeric")
    stop("concavity.parameter must be numeric!")
  #tune = match.arg(tune)
  if (class(nfolds) != "numeric")
    stop("nfolds must be numeric!")
  #type.measure = match.arg(type.measure)


  if (tune == "cv") {
    if (penalty == "lasso" ) {
      cv.fit = cv.glmnet(x, y, family = family, type.measure = type.measure, nfolds = nfolds)
      coef.beta = coef(cv.fit, s = "lambda.1se")
      reg.fit = cv.fit$glmnet.fit
      lambda = cv.fit$lambda.1se
      lambda.ind = which(cv.fit$lambda == cv.fit$lambda.1se)
    } else if (family != 'cox') {
      cv.fit = cv.ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter, nfolds = nfolds)
      cv.1se.ind = min(which(cv.fit$cve<cv.fit$cve[ cv.fit$min]+cv.fit$cvse[ cv.fit$min]))
      coef.beta = cv.fit$fit$beta[, cv.1se.ind]  # extract coefficients at a single value of lambda, including the intercept
      reg.fit = cv.fit$fit

      lambda = cv.fit$lambda[cv.1se.ind]
      lambda.ind = cv.1se.ind
    } else {
      cv.fit = cv.ncvsurv(x, y, family = family, penalty = penalty, gamma = concavity.parameter, nfolds = nfolds)
      cv.1se.ind = min(which(cv.fit$cve<cv.fit$cve[ cv.fit$min]+cv.fit$cvse[ cv.fit$min]))
      coef.beta = cv.fit$fit$beta[, cv.1se.ind]  # extract coefficients at a single value of lambda
      reg.fit = cv.fit$fit

      lambda = cv.fit$lambda[cv.1se.ind]
      lambda.ind = cv.1se.ind
    }
  } else {
    n = nrow(x)
    if (penalty == "lasso" ) {
      reg.fit = glmnet(x, y, family = family)
      coef.beta = rbind(reg.fit$a0,as.matrix(reg.fit$beta))  # extract coefficients at all values of lambda,  including the intercept
      dev = deviance(reg.fit)
      reg.df = reg.fit$df
    } else {
      if(family != 'cox'){

        reg.fit = ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter)
        coef.beta = reg.fit$beta  # extract coefficients at all values of lambda, including the intercept
        dev = loglik.cond(x, y, coef.beta, family = family)
        reg.df = getdf.cond(coef.beta[-1, , drop = FALSE])
      } else {
        reg.fit = ncvsurv(x, y, family = family, penalty = penalty, gamma = concavity.parameter)
        coef.beta = reg.fit$beta  # extract coefficients at all values of lambda, including the intercept
        dev = 2*reg.fit$loss
        reg.df = getdf.cond(coef.beta)
      }
    }

    if (tune == "aic") {
      obj = dev + 2 * reg.df
    }
    if (tune == "bic") {
      obj = dev + log(n) * reg.df
    }
    if (tune == "ebic") {
      obj = dev + log(n) * reg.df + 2 * gamma.ebic * log(choose(dim(x)[2], reg.df))
    }
    lambda.ind = which.min(obj)
    coef.beta = coef.beta[, lambda.ind]
    lambda = reg.fit$lambda[lambda.ind]
  }

  if(family != 'cox'){
    a0 = coef.beta[1]
    coef.beta = coef.beta[-1]
  } else{
    a0 = NULL
    coef.beta = as.vector(coef.beta)
  }
  if(is.null(exposure)){coef.beta <- coef.beta} else{coef.beta <- coef.beta[-ncol(x)]}
  ix = which(coef.beta != 0)
  beta = coef.beta[ix]
  return(list(ix = ix, a0 = a0, beta = beta, fit = reg.fit, lambda = lambda, lambda.ind = lambda.ind))
}
