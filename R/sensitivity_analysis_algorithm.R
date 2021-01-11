#' Sensitivity Assessment to Unmeasured Confounding with Multiple Treatments
#'
#'This function implements the nested multiple imputation using Bayesian Additive Regression Trees (BART)
#'
#' @param covariates Dataframe including all the covariates
#' @param y Numeric vector for the binary outcome
#' @param A Numeric vector for the treatment indicator
#' @param alpha Priors for sensitivity parameters
#' @param n_p Number of nested imputations to conduct
#' @param nposterior Number of posterior samples, default is 1000
#' @param sensitivity_correction Whether to use sensitivity correction algorithm, default is TRUE
#'
#' @return A list of dataframes for each ATE between different treatments.
#' If number of treatments = 3, it contains
#' \item{ATE12:}{A dataframe with number of rows =
#' n_p * nrow(alpha) and number of columns = length(y)}
#' \item{ATE23:}{A dataframe with number of rows =
#' n_p * nrow(alpha) and number of columns = length(y)}
#'\item{ATE13:}{A dataframe with number of rows =
#' n_p * nrow(alpha) and number of columns = length(y)}
#' @export
#'
#' @examples
#'\donttest{sample_size = 10
#'x1 = rbinom(sample_size, 1, prob=0.4)
#'x2 = rbinom(sample_size, 1, prob=0.5)
#'lp.A = 0.2 * x1 + 0.4 * x2 + rnorm(sample_size, 0, 0.1)
#'lp.B = -0.3 * x1 + 0.8 * x2 + rnorm(sample_size, 0, 0.1)
#'lp.C = 0.1 * x1 + 0.5 * x2 + rnorm(sample_size, 0, 0.1)
#' # calculate the true probability of assignment
#'p.A1 <- exp(lp.A)/(exp(lp.A)+exp(lp.B)+exp(lp.C))
#'p.A2 <- exp(lp.B)/(exp(lp.A)+exp(lp.B)+exp(lp.C))
#'p.A3 <- exp(lp.C)/(exp(lp.A)+exp(lp.B)+exp(lp.C))
#'p.A <- matrix(c(p.A1,p.A2,p.A3),ncol = 3)
#'A = NULL
#'for (i in 1:sample_size) { # assign treatment
#'  A[i] <- sample(c(1, 2, 3),
#'                 size = 1,
#'                 replace = TRUE,
#'                 prob = p.A[i, ])
#'}
#'table(A)
#'# set the binary outcome
#'Y2 = 0.3 * x1 + 0.2 * x1 * x2 + 1.3 * x2
#'Y1 = -0.6 * x1 + 0.5 * x2 + 0.3 * x1 * x2
#'Y0 = -0.8 * x1 - 1.2 * x2 + 1.5 * x2 * x1
#'Y2 = rbinom(sample_size, 1, exp(Y2)/(1+exp(Y2)))
#'Y1 = rbinom(sample_size, 1, exp(Y1)/(1+exp(Y1)))
#'Y0 = rbinom(sample_size, 1, exp(Y0)/(1+exp(Y0)))

#'dat = cbind(Y0, Y1, Y2, A)
#'Yobs <- apply(dat, 1, function(x)
#'  x[1:3][x[4]]) #observed when trt is received
#'n = 1
#'alpha = cbind(
#'  runif(n, mean(Y0[A ==1])-mean(Y0[A ==2]) - 0.001, mean(Y0[A ==1])-mean(Y0[A ==2]) + 0.001),
#'  runif(n, mean(Y1[A ==2])-mean(Y1[A ==1]) - 0.001, mean(Y1[A ==2])-mean(Y1[A ==1]) + 0.001),
#'  runif(n, mean(Y1[A ==2])-mean(Y1[A ==3]) - 0.001, mean(Y1[A ==2])-mean(Y1[A ==3]) + 0.001),
#'  runif(n, mean(Y0[A ==1])-mean(Y0[A ==3]) - 0.001, mean(Y0[A ==1])-mean(Y0[A ==3]) + 0.001),
#'  runif(n, mean(Y2[A ==3])-mean(Y2[A ==1]) - 0.001, mean(Y2[A ==3])-mean(Y2[A ==1]) + 0.001),
#'  runif(n, mean(Y2[A ==3])-mean(Y2[A ==2]) - 0.001, mean(Y2[A ==3])-mean(Y2[A ==2]) + 0.001))
#'y <- Yobs
#'n_p <- 1
#'sample_gap <- 10
#'sensitivity_analysis_result <- sensitivity_analysis(cbind(x1, x2), Yobs,
#'A, alpha, n_p = 1, sensitivity_correction = TRUE)
#'mean(sensitivity_analysis_result$ATE_23)
#'mean(sensitivity_analysis_result$ATE_13)
#'mean(sensitivity_analysis_result$ATE_12)}


sensitivity_analysis = function(covariates, y, A, alpha, n_p, nposterior = 1000, sensitivity_correction = TRUE){
  # change the type of y and A as the input parameter of bart function
  covariates = as.matrix(covariates)
  y = as.numeric(y)
  A = as.integer(A)
  A_unique_length <- length(unique(A))
  alpha = as.matrix(alpha)
  n_alpha = nrow(alpha)
  # if we don`t need to adjust the sensitivity, then a simple bart model will be fitted and return a list of ATE
  if (!sensitivity_correction) {
    # fit binary bart
    bart_mod = BART::pbart(x.train = cbind(covariates, A),  y.train = y,  ndpost = nposterior, printevery = 10000)
    for (j in 1:(A_unique_length)){
      assign(paste0("predict_",j), pnorm(BART::pwbart(cbind(covariates, A = j), bart_mod$treedraws)))
    }

    for (i in 1:(A_unique_length-1)){
      for (j in (i + 1):A_unique_length){
        assign(paste0("ATE_",i,j), rowMeans(eval(parse(text =(paste0("predict_",i)))) - eval(parse(text =(paste0("predict_",j))))))
      }
    }
    result <- NULL
    for (i in 1:(A_unique_length-1)){
      for (j in (i + 1):A_unique_length){

        assign(paste0("ATE_",i,j), list(eval(parse(text =(paste0("ATE_",i,j))))))
        assign(paste0("ATE_",i,j), stats::setNames(eval(parse(text =(paste0("ATE_",i,j)))), paste0("ATE_",i,j)))
        result <- c(result, (eval(parse(text =(paste0("ATE_",i,j))))))
      }
    }
    return(result)
  }



  #### from here, the function of sensitivity correction will begin.

  # fit the treatment assigment model, to use gap-sampling, we over sample n * sample_gap samples, and select a sample per sample_gap =10 turns
  assign_mod = BART::mbart2(x.train = covariates, as.integer(as.factor(A)), x.test = covariates, ndpost = n_p * 10, nskip = 1000)

  # assign the estimated assignment probability to each sample, the size is (n, #treatment, sample_size)
  p = array(assign_mod$prob.test[seq(1, nrow(assign_mod$prob.test), 10),], dim = c(n_p, A_unique_length, length(A)))

  # start to calculate causal effect by n * n times
  # for (j in 1:(A_unique_length)){
  #   assign(paste0("causal_effect_",j), matrix(NA, nrow = n_alpha, ncol = n_p * nposterior))
  # }
  step = 1
  train_x = cbind(covariates, A)
  for (m in 1:(A_unique_length-1)){
    for (n in (m + 1):A_unique_length){
      assign(paste0("ATE_",m,n), c())
    }
  }
  for (i in 1:n_alpha) {
    for (m in 1:(A_unique_length-1)){
      for (n in (m + 1):A_unique_length){
        assign(paste0("ate_",m,n), c())
      }
    }
    for (j in 1:n_p) {
      print(paste("step :", step, "/", n_alpha*n_p))
      sort(unique(train_x[, "A"]))
      # correct the binary outcome based on A, alpha, p

        train_y = ifelse(train_x[, "A"] == sort(unique(train_x[, "A"]))[1], y - (unlist(alpha[i, 1]) * p[j, 2, ] + unlist(alpha[i, 4]) * p[j, 3, ]),
                         ifelse(train_x[, "A"] == sort(unique(train_x[, "A"]))[2], y - (unlist(alpha[i, 2]) * p[j, 1, ] + unlist(alpha[i, 3]) * p[j, 3, ]),
                                y - (unlist(alpha[i, 5]) * p[j, 1, ] + unlist(alpha[i, 6]) * p[j, 2, ]))) # A = 3

      # fit the bart model to estimate causal effect
      bart_mod = BART::wbart(x.train = cbind(covariates, A),  y.train = train_y,  ndpost = nposterior, printevery = 10000)

      for (z in 1:(A_unique_length)){
        assign(paste0("predict_",z), pnorm(BART::pwbart(cbind(covariates, A = z), bart_mod$treedraws)))
      }

      for (m in 1:(A_unique_length-1)){
        for (n in (m + 1):A_unique_length){
          assign(paste0("ate_",m,n), c(eval(parse(text =(paste0("ate_",m,n)))), rowMeans(eval(parse(text =(paste0("predict_",m)))) - eval(parse(text =(paste0("predict_",n)))))))
        }
      }

      step = step + 1
    }
    for (m in 1:(A_unique_length-1)){
      for (n in (m + 1):A_unique_length){
        assign(paste0("ATE_",m,n), rbind(eval(parse(text =(paste0("ATE_",m,n)))),
                                         eval(parse(text =(paste0("ate_",m,n))))))
      }
    }
  }
  # return the result
  result <- NULL
  for (i in 1:(A_unique_length-1)){
    for (j in (i + 1):A_unique_length){

      assign(paste0("ATE_",i,j), list(eval(parse(text =(paste0("ATE_",i,j))))))
      assign(paste0("ATE_",i,j), stats::setNames(eval(parse(text =(paste0("ATE_",i,j)))), paste0("ATE_",i,j)))
      result <- c(result, (eval(parse(text =(paste0("ATE_",i,j))))))
    }
  }
  return(result)
}
