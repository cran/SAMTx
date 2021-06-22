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
#'for (m in 1:sample_size) { # assign treatment
#'  A[m] <- sample(c(1, 2, 3),
#'                 size = 1,
#'                 replace = TRUE,
#'                 prob = p.A[m, ])
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
#'mean(sensitivity_analysis_result$ATE_12)
#'mean(sensitivity_analysis_result$ATE_02)
#'mean(sensitivity_analysis_result$ATE_01)}


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
    predict_1 = pnorm(BART::pwbart(cbind(covariates, A = sort(unique(A))[1]), bart_mod$treedraws))
    predict_2 = pnorm(BART::pwbart(cbind(covariates, A = sort(unique(A))[2]), bart_mod$treedraws))
    predict_3 = pnorm(BART::pwbart(cbind(covariates, A = sort(unique(A))[3]), bart_mod$treedraws))
    causal_effect_1 = rowMeans(predict_1 - predict_2)
    causal_effect_2 = rowMeans(predict_2 - predict_3)
    causal_effect_3 = rowMeans(predict_1 - predict_3)
    return(list(ATE_01 = causal_effect_1, ATE_12 = causal_effect_2, ATE_02 = causal_effect_3))
  }
   #### from here, the function of sensitivity correction will begin.

  # Algorithm 1.1: Fit the multinomial probit BART model to the treatment A
  A_model = BART::mbart2(x.train = covariates, as.integer(as.factor(A)), x.test = covariates, ndpost = n_p * 10, nskip = 1000)


  # Algorithm 1.1: Estimate the generalized propensity scores for each individual
  p = array(A_model$prob.test[seq(1, nrow(A_model$prob.test), 10),], dim = c(n_p, 3, length(A)))

  # Algorithm 1.2: start to calculate causal effect by M1 * M2 times
  causal_effect_1 = matrix(NA, nrow = n_alpha * n_p, ncol = nposterior)
  causal_effect_2 = matrix(NA, nrow = n_alpha * n_p, ncol = nposterior)
  causal_effect_3 = matrix(NA, nrow = n_alpha * n_p, ncol = nposterior)
  step = 1
  train_x = cbind(covariates, A)
  for (j in 1:n_p) {
    # Algorithm 1.2: Draw M1 generalzied propensity scores from the posterior predictive distribution of the A model for each individual
    p_draw_1 <- p[j, 1, ]
    p_draw_2 <- p[j, 2, ]
    p_draw_3 <- p[j, 3, ]
    for (m in 1:n_alpha) {
      # Algorithm 1.2: Draw M2 values from the prior distribution of each of the sensitivity paramaters alpha for eacg treatment
      print(paste("step :", step, "/", n_alpha*n_p))
      sort(unique(train_x[, "A"]))
      # Algorithm 1.3: Compute the adjusted outcomes y_corrected for each treatment for each M1M2 draws
      y_corrected = ifelse(
        train_x[, "A"] == sort(unique(train_x[, "A"]))[1],
        y - (unlist(alpha[m, 1]) * p_draw_2  + unlist(alpha[m, 4]) * p_draw_3),
        ifelse(
          train_x[, "A"] == sort(unique(train_x[, "A"]))[2],
          y - (unlist(alpha[m, 2]) * p_draw_1 + unlist(alpha[m, 3]) * p_draw_3),
          y - (unlist(alpha[m, 5]) * p_draw_1 + unlist(alpha[m, 6]) * p_draw_2)
        )
      )
      # Algorithm 1.4: Fit a BART model to each set of M1*M2 sets of observed data with the adjusted outcomes y_corrected
      bart_mod = BART::wbart(x.train = cbind(covariates, A),  y.train = y_corrected,  ndpost = nposterior, printevery = 10000)

      predict_1 = BART::pwbart(cbind(covariates, A = sort(unique(A))[1]), bart_mod$treedraws)
      predict_2 = BART::pwbart(cbind(covariates, A = sort(unique(A))[2]), bart_mod$treedraws)
      predict_3 = BART::pwbart(cbind(covariates, A = sort(unique(A))[3]), bart_mod$treedraws)
      causal_effect_1[((m - 1) * n_p) + j, ] = rowMeans(predict_1 - predict_2)
      causal_effect_2[((m - 1) * n_p) + j, ] = rowMeans(predict_2 - predict_3)
      causal_effect_3[((m - 1) * n_p) + j, ] = rowMeans(predict_1 - predict_3)
      step = step + 1
    }

  }
  # Final combined adjusted causal effect
  list(ATE_01 = causal_effect_1, ATE_12 = causal_effect_2, ATE_02 = causal_effect_3)
}
