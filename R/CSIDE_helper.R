choose_sigma_gene <- function(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode,
                              verbose = F, n.iter = 100, MIN_CHANGE = 0.001,
                              MAX_ITER_SIGMA = 10, PRECISION.THRESHOLD = .01,
                              spase = F) {
  sigma_s_best <- sigma_init
  if(spase){
    phi_vals <- (5:95)/100
    sigma_vals <- (1 - phi_vals) / phi_vals
    sigma_s_best <- "1"
  } else
    sigma_vals <- names(Q_mat_all)
  alpha1 <- NULL; alpha2 <- NULL;
  MIN_ITERATIONS <- 15
  n.iter.tot <- 0
  for(iter in 1:MAX_ITER_SIGMA) {
    last_sigma <- sigma_s_best
    sigma_curr <- as.character(sigma_s_best)
    if(!spase)
      set_likelihood_vars(Q_mat_all[[sigma_curr]], X_vals, sigma = sigma_curr)

    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI, test_mode, verbose = verbose, n.iter = n.iter,
                                  MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                                  alpha1_init = alpha1, alpha2_init = alpha2,
                                  MIN_ITERATIONS = MIN_ITERATIONS, r = as.numeric(sigma_curr), spase = spase)
    n.iter.tot <- n.iter.tot + res$n.iter
    alpha1 <- res$alpha1; alpha2 <- res$alpha2
    lambda <- as.vector(res$prediction)
    res_val <- numeric(length(sigma_vals))
    names(res_val) <- as.character(sigma_vals)
    for(sigma_s in sigma_vals) {
      if(!spase)
        set_likelihood_vars(Q_mat_all[[as.character(sigma_s)]], X_vals, sigma = sigma_s)
      res_val[as.character(sigma_s)] <- calc_log_l_vec(lambda, as.vector(t(Y)), spase = spase, nUMI = nUMI, r = sigma_s)
    }
    sigma_s_best <- names(which.min(res_val))
    if(sigma_s_best == last_sigma) {
      break
    }
  }
  res$n.iter <- n.iter.tot
  return(list(sigma_s_best = sigma_s_best, res = res))
}

mysweept <- function(tX2,tlk, K) {
  g_2 <- tX2[rep(1:dim(tX2)[1],K),] * tlk[rep(1:K, each = dim(tX2)[1]), ]
  return(g_2)
}

sweep1t <- function(tX1, lambda) {
  g_1 <- Rfast::eachrow(tX1, lambda,oper = "*")
  return(g_1)
}

sweep2t <- function(tX1, tdl, k) {
  if(dim(tX1)[1] > 0) {
    X1_Q <- Rfast::eachrow(tX1, tdl[k,],oper = "*")
  } else {
    X1_Q <- tX1
  }
  return(X1_Q)
}

sweep3t_all <- function(tX2, tdl, K) {
  X2_Q <- tX2[rep(1:dim(tX2)[1],K),] * tdl[rep(1:K, each = dim(tX2)[1]), ]
  return(X2_Q)
}

construct_hess_fast <- function(X1,X2,lambda,lambda_k, K, d1_d2,
                                dlambda_k = NULL, dlambda = NULL,
                                tlambda_k = NULL, tlambda = NULL, spase = F) {
  tX1 <- t(X1); tX2 <- t(X2)
  if(spase) {
    g_1 <- sweep1t(tX1, dlambda)
    tlk <- t(tlambda_k)
    g_2 <- mysweept(tX2,t(dlambda_k), K)
    X1_Q <- Rfast::eachrow(tX1, tlambda*d1_d2$d1_vec, '*')
  } else {
    g_1 <- sweep1t(tX1, lambda)
    tlk <- t(lambda_k)
    g_2 <- mysweept(tX2,tlk, K)
    X1_Q <- Rfast::eachrow(tX1, lambda*d1_d2$d1_vec, '*')
  }
  tdl <- Rfast::eachrow(tlk, d1_d2$d1_vec, '*')
  grad <- rbind(g_1, g_2)
  grad_Q <- Rfast::eachrow(grad, d1_d2$d2_vec, '*')
  H1 <- -grad_Q %*% t(grad)
  grad_1 <- matrix(rowSums(X1_Q), 1, dim(X1)[2])
  H2_11 <- X1_Q %*% X1
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  H2_12 <- matrix(0,nrow = L1, ncol = L2*K)
  for(k in 1:K) {
    X1_Q <- sweep2t(tX1, tdl, k)
    H2_12[,(1 + L2*(k-1)):(L2*k)] <- X1_Q %*% X2
  }
  H2 <- matrix(0, nrow = L1 + L2*K, ncol = L1 + L2*K)
  if(L1 > 0) {
    H2[(1:L1),1:L1] <- H2_11
    H2[1:L1,(L1+1):(L1+L2*K)] <- H2_12
    H2[(L1+1):(L1+L2*K),1:L1] <- t(H2_12)
  }
  X2_Q <- sweep3t_all(tX2, tdl, K)
  grad_2 <- matrix(rowSums(X2_Q),dim(X2)[2],K)
  H2m <- X2_Q %*% X2
  for(k in 1:K) {
    H2[(L1 + 1 + (k-1)*L2):(L1 + k*L2), (L1 + 1 + (k-1)*L2):(L1 + k*L2)] <-
      H2m[(1 + (k-1)*L2):(k * L2),]# X2_Q %*% X2
  }
  H <- (H1 - H2)
  return(list(H = H, grad_1 = grad_1, grad_2 = grad_2))
}

solveIRWLS.effects_trust <- function(Y, X1, X2, my_beta, test_mode, verbose = FALSE,
                                     n.iter = 200, MIN_CHANGE = .01, PRECISION.THRESHOLD = .05,
                                     alpha1_init = NULL, alpha2_init = NULL, MIN_ITERATIONS = 15,
                                     spase = F, nUMI = NULL, r = NULL){
  lam_threshold = 1e-8;
  beta_succ <- 1.1; beta_fail <- 0.5;
  gamma <- 0.8 # prev gamma = 0.1 / 0.25
  epsilon_2 <- 1e-6 * length(Y);
  delta <- 0.1 #trust region radius
  if(spase)
    init_val <- 0
  else {
    init_val <- min(-5, log(10/median(rowSums(my_beta))))
    Y[Y > K_val] <- K_val
  }
  L1 <- dim(X1)[2]; L2 <- dim(X2)[2]
  n_cell_types = dim(my_beta)[2]
  if(is.null(alpha1_init))
    alpha1 <- numeric(dim(X1)[2])
  else
    alpha1 <- alpha1_init
  if(is.null(alpha2_init)) {
    alpha2 <- matrix(0,nrow = dim(X2)[2], ncol = n_cell_types) # initialize it to be the previous cell type means
    alpha2[1,] <- init_val
    if(test_mode == 'categorical') {
      alpha2[,] <- init_val # multi mode
    }
  } else {
    alpha2 <- alpha2_init
  }
  pred_decrease_vals <- numeric(n.iter)
  K <- dim(my_beta)[2]
  S_X_a <- sweep(X2 %*% (alpha2), 1, X1 %*% (alpha1),'+')
  S_X_a <- pmax(-50,pmin(50,S_X_a))
  if(spase) {
    lambda_k <- expit(S_X_a)*my_beta #J by K
  } else
    lambda_k <- exp(S_X_a)*my_beta #J by K
  lambda <- rowSums(lambda_k)
  lambda[lambda < lam_threshold] <- lam_threshold
  if(spase)
    lambda[lambda > 1 - lam_threshold] <- 1 - lam_threshold
  d1_d2 <- get_d1_d2(Y, lambda, spase = spase, nUMI = nUMI, r = r)
  if(spase) {
    dlambda_k <- dexpit(S_X_a)*my_beta #J by K
    dlambda <- rowSums(dlambda_k)
    tlambda_k <- texpit(S_X_a)*my_beta #J by K
    tlambda <- rowSums(tlambda_k)
    prev_ll <- calc_log_l_vec.spASE(lambda,Y, nUMI, r)
  } else {
    prev_ll <-  -sum(d1_d2$d0_vec)
    dlambda_k <- NULL; dlambda <- NULL
    tlambda_k <- NULL; tlambda <- NULL
  }
  error_vec <- (1:dim(my_beta)[2]) == 0
  for(itera in 1:n.iter) {
    H_list <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2,
                                  dlambda_k = dlambda_k, dlambda = dlambda,
                                  tlambda_k = tlambda_k, tlambda = tlambda, spase = spase)
    H <- H_list$H
    if(spase) {
      grad_1 <- d1_d2$d1_vec %*% sweep(X1, 1,dlambda,'*')
      d1_lam <- sweep(dlambda_k, 1, d1_d2$d1_vec, '*')
      grad_2 <- t(X2) %*% d1_lam
    } else {
      grad_1 <- H_list$grad_1; grad_2 <- H_list$grad_2
    }
    d_vec_o <- c(grad_1, grad_2)
    #if (itera == 55) {
    #  ca = 1
    #}
    D_mat_o <- psd(H)
    norm_factor <- norm(D_mat_o,"2")
    D_mat <- D_mat_o / norm_factor
    d_vec <- (d_vec_o / norm_factor) #- 2 * lambda_reg * c(alpha1, alpha2)
    #epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
    A <- cbind(diag(dim(D_mat)[2]), -diag(dim(D_mat)[2]))
    bzero <- rep(-delta,2*dim(D_mat)[2]); lambda_reg <- 1e-7
    D_mat <- D_mat + diag(dim(D_mat)[1])*lambda_reg # avoid numerical errors
    solution <-  quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution #* 0.01 # CHANGE THIS

    predicted_decrease = -(0.5*t(solution) %*% D_mat_o %*% solution - sum(d_vec_o*solution))

    alpha1_new <- alpha1 + solution[1:L1]
    alpha2_new <- alpha2 + matrix(solution[(L1 + 1):length(solution)], nrow = L2, ncol = K)
    S_X_a <- sweep(X2 %*% (alpha2_new), 1, X1 %*% (alpha1_new),'+')
    S_X_a <- pmax(-50,pmin(50,S_X_a))
    if(spase) {
      lambda_k_new <- expit(S_X_a)*my_beta #J by K
    } else
      lambda_k_new <- exp(S_X_a)*my_beta #J by K
    error_vec <- is.na(colMeans(lambda_k_new)) | error_vec
    if(spase)
      lambda_k_new[is.na(lambda_k_new)] <- 0.5
    else
      lambda_k_new[is.na(lambda_k_new)] <- 1
    lambda_new <- rowSums(lambda_k_new)
    lambda_new[lambda_new < lam_threshold] <- lam_threshold
    if(spase) {
      lambda_new[lambda_new > 1 - lam_threshold] <- 1 - lam_threshold
      log_l <- calc_log_l_vec.spASE(lambda_new,Y, nUMI, r)
    } else {
      d_all <- get_d1_d2(Y, lambda_new, spase = spase, nUMI = nUMI, r = r)
      log_l <-  -sum(d_all$d0_vec)
    }
    true_decrease = prev_ll - log_l

    pred_decrease_vals[itera] <- predicted_decrease
    if(true_decrease >= (gamma) * predicted_decrease) {
      delta = beta_succ*delta
      alpha1 <- alpha1_new
      alpha2 <- alpha2_new
      lambda_k <- lambda_k_new
      lambda <- lambda_new
      lambda[lambda < lam_threshold] <- lam_threshold
      prev_ll <- log_l
      if(spase) {
        d1_d2 <- get_d1_d2(Y, lambda, spase = spase, nUMI = nUMI, r = r)
        dlambda_k <- dexpit(S_X_a)*my_beta #J by K
        dlambda <- rowSums(dlambda_k)
        tlambda_k <- texpit(S_X_a)*my_beta #J by K
        tlambda <- rowSums(tlambda_k)
      } else
        d1_d2 <- d_all
    } else {
      delta = min(1,beta_fail*delta)
    }
    if(delta < MIN_CHANGE || (itera >= MIN_ITERATIONS &&
                              max(pred_decrease_vals[(itera - MIN_ITERATIONS+1):itera]) < min(epsilon_2))) {
      break
    }
  }
  H <- construct_hess_fast(X1,X2,lambda,lambda_k, K, d1_d2,
                           dlambda_k = dlambda_k, dlambda = dlambda,
                           tlambda_k = tlambda_k, tlambda = tlambda, spase = spase)$H
  eps = 1e-6
  if(min(eigen(H)$values) < eps) {
    I <- solve(psd(H, epsilon = eps))
  } else {
    I <- solve(H)
  }
  precision <- abs(solve(D_mat) %*% d_vec)
  precision[is.na(diag(I))] <- PRECISION.THRESHOLD + 100
  converged_vec <- check_converged_vec(X1,X2,my_beta, itera, n.iter,
                                       error_vec, precision, PRECISION.THRESHOLD)
  names(error_vec) <- colnames(my_beta)
  return(list(alpha1 = alpha1, alpha2 = alpha2, converged = any(converged_vec), I = I, d = d_vec_o,
              n.iter = itera, log_l = prev_ll, precision = precision, prediction = lambda,
              converged_vec = converged_vec, error_vec = error_vec))
}

estimate_effects_trust <- function(Y, X1, X2, my_beta, nUMI, test_mode, verbose = F,
                                   n.iter = 200, MIN_CHANGE = 1e-3, PRECISION.THRESHOLD = 0.05,
                                   alpha1_init = NULL, alpha2_init = NULL, MIN_ITERATIONS = 15, spase = F, r = NULL) {
  if(!spase)
    my_beta <- sweep(my_beta,1, nUMI, '*')
  solveIRWLS.effects_trust(Y,X1,X2, my_beta, test_mode, verbose = verbose,
                           n.iter = n.iter, MIN_CHANGE = MIN_CHANGE,
                           PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                           alpha1_init = alpha1_init, alpha2_init = alpha2_init,
                           MIN_ITERATIONS = MIN_ITERATIONS, spase = spase, r = r, nUMI = nUMI)
}

estimate_gene_wrapper <- function(Y,X1,X2,my_beta, nUMI, sigma_init, test_mode, verbose = F, n.iter = 200, MIN_CHANGE = 1e-3,
                                  sigma_gene = T, PRECISION.THRESHOLD = 0.05, spase = F) {
  if(sigma_gene)
    return(choose_sigma_gene(sigma_init, Y, X1, X2, my_beta, nUMI,test_mode, verbose = verbose,
                             n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD, spase = spase))
  else {
    res <- estimate_effects_trust(Y,X1,X2,my_beta, nUMI,test_mode, verbose = verbose,
                                  n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, PRECISION.THRESHOLD = PRECISION.THRESHOLD, spase = spase, r = as.numeric(sigma_init))
    return(list(sigma_s_best = -1, res = res))
  }
}
