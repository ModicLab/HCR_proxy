#Lasso and Selective inference

cvfit_wrapper <- function(mtx, rhs, maxit_tries = 3) {
  print("### Initialising cvfit_wrapper!")
  #returns best cvfit under circumstances
  cvfit <- NULL
  thresh <- 1E-28
  torep <- T
  i <- 1
  while (torep == T & i <= maxit_tries) {
    print(paste0("### Initialising cvfit_wrapper at iteration = ",i," !"))
    torep <- F
    maxit <- 1E7 * 10^(i-1)
    i <- i + 1
    FLI.warning <- NA_character_
    warning.found <- F
    
    withCallingHandlers({
      #print(paste0("cvfit at maxit=",maxit,", thresh=",as.character(thresh)))
      set.seed(1337)
      print("### Initialising cv.glmnet!")
      cvfit <- cv.glmnet(mtx, rhs, standardize = FALSE, thresh = thresh, maxit = maxit,
                         # lambda = seq(0.001,1,by=0.001),
                         nfolds = 10)
    },
    warning = function(w){
      if (grepl("lambda value not reached after maxit", conditionMessage(w))){
        print(conditionMessage(w))
        torep <- T
      } else {
        print(conditionMessage(w))
      }
    }
    )
  }
  cvfit$thresh <- thresh
  cvfit$maxit <- maxit
  return(cvfit)
}

fixedLassoInf_wrapper <- function(mtx, rhs, alphax, tailarea_rtol, n_tries, scales) {
  print("### Initialising fixedLassoInf_wrapper!")
  cvfit <- cvfit_wrapper(mtx,rhs, maxit_tries = 5)
  #return(cvfit)
  if(is.null(cvfit)){
    print("###### cv.glmnet method not converged!")
  } else {
    print("### cv.glmnet method successfully converged!")
  }
  
  # #####EXPERIMENTAL#######
  # if(cvfit[["cvm"]][cvfit[["lambda"]]==min(cvfit[["lambda"]])]-1*cvfit[["cvsd"]][cvfit[["lambda"]]==min(cvfit[["lambda"]])] <
  #    (cvfit[["cvm"]][cvfit[["lambda"]]==cvfit[["lambda.min"]]])) {
  #   lambda.min <- min(cvfit[["lambda"]])
  # } else {
  #   lambda.min <- cvfit$lambda.min
  # }
  #lambda.min <- cvfit$lambda.1se
  
  lambda.min <- cvfit$lambda.min
  print(paste0("Lambda (",round(min(cvfit[["lambda"]]),3), "<", round(lambda.min,5), "<",round(max(cvfit[["lambda"]]),3), ") in range: ", ifelse(cvfit$lambda.min %in% cvfit[["lambda"]], "TRUE","FALSE")))
  
  estimates.min <- coef(cvfit[["glmnet.fit"]], s = lambda.min, exact = TRUE, x = mtx,
                        y = rhs, thresh = cvfit$thresh, maxit = cvfit$maxit, standardize = FALSE)

  estimates.scales <- ifelse(rownames(estimates.min) %in% names(scales), scales[rownames(estimates.min)], 1.0)
  rescaled.glmnet.estimates.min <- estimates.min / estimates.scales
  
  temp_tol.kkt <- 0.01
  res.min <- NULL
  bits <- 100
  FLI.error <- T
  if(!all(estimates.min[-1] == 0)){
    for (i in 1:n_tries) {
      print(paste0("### [lambda.min] Running fixedLassoInf at bits = ", as.character(bits)))
      FLI.warning <- NA_character_
      withCallingHandlers({
        withRestarts({
          #print(lambda.min * length(rhs))
          set.seed(1337)
          #print(lambda.min*length(rhs))
          #print(lambda.min)
          #print(cvfit[["lambda"]])
          #Here, lambda=lambda.min*length(rhs) is recommended but incompatible with calculations
          res.min <- fixedLassoInf(x=mtx, y=rhs, beta=estimates.min[-1], lambda=lambda.min, verbose = F,
                                   tol.beta = 0.01,
                                   alpha = alphax, tol.kkt=temp_tol.kkt, bits=bits, family = "gaussian",
                                   sigma = estimateSigma(mtx, rhs, intercept=TRUE, standardize=FALSE)$sigmahat)
          FLI.error <- F
          res.min$converged <- F
        }, muffleStop=function() {
          message("fixedLassoInf not converged! - muffleStop")
        })
      },
      warning = function(w){FLI.warning <- conditionMessage(w)},
      error = function(e){
        # FLI.error <- conditionMessage(e)
        # print(FLI.error)
        invokeRestart("muffleStop")
      })
      
      print(FLI.error)
      res.min$tailarea_rerrs <- abs(res.min$tailarea - 0.5*alphax) / (0.5*alphax)
      res.min$n.tailarea_tries <- i
      res.min$bits.tailarea_tries <- bits
      
      print(res.min$tailarea)
      print(res.min$tailarea_rerrs)
      
      if (!all(res.min$tailarea_rerrs < tailarea_rtol)){
        print("### [lambda.min] !all(res.min$tailarea_rerrs < tailarea_rtol)")
        bits <- bits + 200 # increase precision for the next try
        res.min <- NULL
      } else if (grepl("Solution beta does not satisfy the KKT conditions", FLI.warning)){
        print("### [lambda.min] Solution beta does not satisfy the KKT conditions")
        bits <- bits + 200 # increase precision for the next try
        res.min <- NULL
        #computes new cvfit with lower thresh parameter?
      } else {
        if(!is.na(FLI.warning)){
          print(FLI.warning)
        }
        if(FLI.error){
          print("### [lambda.min] fixedLassoInf NOT converged!")
          res.min$converged <- F
        } else {
          print("### [lambda.min] fixedLassoInf converged!")
          res.min$converged <- T
        }
        break
      }
    }
  } else {
    print("### [lambda.min] All coefs == 0!")
    res.min$converged <- T
  }
  return (list(cvf = cvfit, at_lambda_min=list(lambda.min = lambda.min, estimates.min=estimates.min,rescaled.glmnet.estimates.min=rescaled.glmnet.estimates.min, res.min=res.min)))
}


