library("BayesianTools")
library("corrplot")
library("readr")
library("stats")
library(IDPmisc)
library(coda)
library(sm)
library(zoo)
library("vioplot")
library(cairoDevice)

#functions--------------------------------------------------
rm(list = ls())
##generate initial condition--------------------------------------------------
generate_init_condi <- function(r0,
                                Di = 2.9,
                                Dp = 2.3,
                                De = 2.9,
                                Dq = c(21,13,6,5,3),
                                alpha = 0.55,
                                Dh = 30,
                                N = 8600000,
                                flowN = c(339395, 340112, 341185, 341185, 341185)
) {
  
  stopifnot(r0>=0 & r0<=1 & Di>=0 & Dp>=0 & De>=0 & all(Dp>=0) & alpha>=0 & alpha<=1 & Dh>=0 & N>=0 & all(flowN>=0))
  
  ## N            : population size
  ## H0           : initial number of hospitalized cases based on the reports
  ## R0           : initial number of removed individuals
  ## De           : latent period
  ## r0           : initial ascertainment rate
  ## realData     : real data from openZH
  R0 <- 0
  H0 <- 27
  
  #where -----/ is the path location of the data file
  #where ...../ is the path location where you want to see the plots
  
  realData_all <-read.csv("-----/cases.csv")[ ,c( 'CH')]
  daily_new_case_all<- rep(realData_all[1],264)
  for (i in 2:264) 
  {
    daily_new_case_all[i] <-realData_all[i]- realData_all[i-1]
  }
  #smoothing
  d<-lowess(daily_new_case_all, f=10/264)
  daily_new_case_all<-d$y
  #cases for the 25th of February
  daily_new_case <- daily_new_case_all[-c(1,2,3,4,5,6,7)]
  #start on the 4th of March to have sufficient number of infections to run the code
  
  mar4_idx = 8
  
  ##
  E0 <- sum(daily_new_case_all[(mar4_idx+round(Dp)):(mar4_idx+round(Dp)+round(De)-1)]) / r0 
  P0 <- sum(daily_new_case_all[mar4_idx:(mar4_idx+round(Dp)-1)]) / r0                    
  I0 <- sum(daily_new_case_all[(mar4_idx-round(Di)):(mar4_idx-1)])   
  A0 <- I0 * (1 - r0) / r0
  S0 <- N - E0 - P0 - I0 - A0 - H0 - R0
  init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, H = H0, R = R0), 0)
  
  ## helper function
  # transform variables to a form that SEIRpred can use
  # so that SEIRpred can be re-used as much as possible
  transform_var_main_5stage=function(pars) {
    
    b_vec <- pars[1:5]
    ##
    r1 <- pars[6]
    r2 <- 1 / (1 + (1 - r1) / (r1 * exp(pars[7])))
    r3 <- 1 / (1 + (1 - r2) / (r2 * exp(pars[8])))
    r4 <- 1 / (1 + (1 - r3) / (r3 * exp(pars[9])))
    r5 <- 1 / (1 + (1 - r4) / (r4 * exp(pars[10])))
    r_vec <- c(r1,r2,r3,r4,r5)
    
    return(list(b_vec, r_vec))
  }
  
  return(list(Di=Di,
              Dp=Dp,
              De=De,
              Dq=Dq,
              alpha=alpha,
              Dh=Dh,
              N=N,
              flowN=flowN,
              daily_new_case = daily_new_case, 
              daily_new_case_all = daily_new_case_all, 
              init_states = init_states,
              days_to_fit=1:255,
              stage_intervals=list(
                c(start=1, end=23),
                c(start=24, end=100),
                c(start=101, end=195),
                c(start=196, end=241),
                c(start=242, end=235)
              ),
              var_trans_fun=transform_var_main_5stage,
              par_lower = c(b1= 0, b2 = 0, b3 = 0, b4 = 0, b5 = 0,  r1 = 0, delta3 = -10, delta4 = -10, delta5 = -10, delta6 = -10),
              par_upper = c(b1 = 2, b2 = 2, b3 = 2, b4 = 2,b5 = 2,  r1 = 1, delta3 = 10, delta4 = 10, delta5 = 10, delta6 = 10)))
  # TODO: please confirm the following:
  # boundaries for delta3-5 will not be used, they are here merely to meet the formality imposed by runMCMC
}

# get_init_sets_list is an alias of generate_init_condi in order not to break exsiting code
get_init_sets_list = generate_init_condi

delta_mean <- 0
delta_sd <- 1
beta_shape1 <- 5
beta_shape2 <- 10

# delta_mean <- 0
# delta_sd <- 1
# beta_shape1 <- 7.3
# beta_shape2 <- 24.6

##SEIR Fitting--------------------------------------------------

default_pars_density <- function(pars) {
  d_vec <- rep(NA,10)
  ##b1,b2, b3, b4
  for(i in c(1:5)) {
    d_vec[i] <- dunif(pars[i], 0, 2, log = T)
  }
  ## r1
  r1 <- pars[6]
  d_vec[6] <- dbeta(r1, beta_shape1, beta_shape2, log = T)
  ## r2
  delta_3 <- pars[7]
  d_vec[7] <- dnorm(delta_3, delta_mean, delta_sd, log = T)
  ## r3
  delta_4 <- pars[8]
  d_vec[8] <- dnorm(delta_4, delta_mean, delta_sd, log = T)
  ## r4
  delta_5 <- pars[9]
  d_vec[9] <- dnorm(delta_5, delta_mean, delta_sd, log = T)
  ## r5
  delta_6 <- pars[10]
  d_vec[10] <- dnorm(delta_6, delta_mean, delta_sd, log = T)

  return(sum(d_vec))
}

default_pars_sampler <- function(n = 1) {
  s_vec <- matrix(NA, n,10)
  ## b1,b2, b3, b4,b5,b6
  for(i in c(1:5)) {
    s_vec[, i] <- runif(n, 0, 2) 
  }
  ## r1 
  s_vec[, 6] <- rbeta(n, beta_shape1, beta_shape2)
  ## r2
  s_vec[, 7] <- rnorm(n, delta_mean, delta_sd)
  ## r3
  s_vec[, 8] <- rnorm(n, delta_mean, delta_sd)
  ## r4
  s_vec[, 9] <- rnorm(n, delta_mean, delta_sd)
  ## r5
  s_vec[, 10] <- rnorm(n, delta_mean, delta_sd)
  
  
  return(s_vec)
}

SEIRfitting=function(init_sets_list, 
                     randomize_startValue=F, 
                     startValue=NA, 
                     output_ret=T, 
                     run_id=0, 
                     skip_MCMC=F, 
                     panel_B_R_ylim=6,
                     plot_combined_fig=T,
                     pars_density=default_pars_density,
                     pars_sampler=default_pars_sampler,
                     pars_name=c("b1", "b2", "b3", "b4","b5","r1", "delta3", "delta4", "delta5","delta6"),
                     calc_clearance=F,
                     n_burn_in=0,
                     n_iterations=10000) {
  if (randomize_startValue & !is.na(startValue)) {
    print("startValue will be ignored since randomize_startValue is set to TRUE!")
  } else if (!randomize_startValue & is.na(startValue)) {
    print("Please specify a startValue since you have set randomize_startValue to FALSE! Exiting!")
    q(save="no")
  }
  
  onset_obs <- init_sets_list$daily_new_case
  init_states <- init_sets_list$init_states
  n_pars = length(pars_name)
  n_stage = length(init_sets_list$stage_intervals)
  
  
  
  ################################################################################################################
  ## likelihood function
  NegBinModel<-function(y,mu,nu)
  {
    x<-lgamma(y+nu)-lgamma(nu)-lgamma(y+1)+nu*log(nu)+y*log(mu)-(nu+y)*log(mu+nu)
    return(exp(x))
  }
  
  loglh_func <- function(pars){
    ypred <- SEIRpred(pars, init_settings = init_sets_list)
    ypred <- ypred[, "Onset_expect"]
    
    # meant to suppress warnings when ypred is negative
    suppressWarnings(p <- NegBinModel(onset_obs, ypred,1))
    
    if(any(p == 0) || any(is.nan(p))){
      logL <- -Inf
    }else{
      logL <- sum(log10(p))
    }
    return(logL)
  }
  ##1
  u<-'u'
  write_file(u, paste0("...../1.txt"))
  # loglh_func(pars = c(1.4, 0.4, 0.1, 0.1, 0.5, -1, 0, 0))
  
  ## Create BayesianSetup and settings, lower/upper for parameters:
  
  pars_prior <- createPrior(density = pars_density, sampler = pars_sampler, 
                            lower = init_sets_list$par_lower, upper = init_sets_list$par_upper)
  
  if (!skip_MCMC) {
    bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
    
    if (randomize_startValue) {
      startValue=pars_sampler()
      while (is.infinite(loglh_func(startValue))) {
        startValue=pars_sampler()
      }
    }
    write_file(u, paste0("...../2.txt"))
    ## DRAM: Adaptive MCMC, prior optimization, delayed rejection
    mh_settings = list(startValue = startValue,
                       adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10)
    write_file(u, paste0("...../21.txt"))
    mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)
    #plot(mh_out)
    write_file(u, paste0("...../22.txt"))
    mcmc_pars_estimate <- getSample(mh_out, start = n_burn_in+2, thin = 1)  ## set start = 2002 as the burn in period
    write_file(u, paste0("...../23.txt"))
    mcmc_pars_estimate <- round(mcmc_pars_estimate, 3)
    write_file(u, paste0("...../24.txt"))
    
    colnames(mcmc_pars_estimate) <- pars_name
    
    if (output_ret) {
      write.table(mcmc_pars_estimate, paste0("...../pars_est_run_",run_id,".txt"), quote = F, row.names = F, sep = "\t")
    }
  } else {
    mcmc_pars_estimate = read.table(paste0("...../pars_est_run_",run_id,".txt"), header = T)
    pars_name = names(mcmc_pars_estimate)
  }
  write_file(u, paste0("...../3.txt"))
  summary_string=paste0(paste(pars_name, collapse = ","), "\n")
  write_file(u, paste0("...../4.txt"))
  par_str=list()
  for (i_par in 1:n_pars) {
    par_str[[i_par]]=paste0(round(mean(mcmc_pars_estimate[,i_par]),2), " (",
                            round(quantile(mcmc_pars_estimate[,i_par],0.025),2)," - " ,
                            round(quantile(mcmc_pars_estimate[,i_par],0.975),2), ")")
  }
  
  write_file(u, paste0("...../5.txt"))
  summary_string = paste0(summary_string, paste(par_str,collapse = ", "),"\n\n")
  estRt_mat <- apply(mcmc_pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_sets_list))
  write_file(u, paste0("...../6.txt"))
  
  summary_string = paste0(summary_string, paste0("stage",1:n_stage,collapse=","), "\n")
  
  r_str=list()
  
  if (n_stage>1) {
    for (i_stage in 1:n_stage) {
      r_str[[i_stage]]=paste0(round(mean(estRt_mat[i_stage,]),2), " (",
                              round(quantile(estRt_mat[i_stage,],0.025),2)," - " ,
                              round(quantile(estRt_mat[i_stage,],0.975),2), ")")
    }
  } else {
    r_str[[1]]=paste0(round(mean(estRt_mat),2), " (",
                      round(quantile(estRt_mat,0.025),2)," - " ,
                      round(quantile(estRt_mat,0.975),2), ")")
  }
  write_file(u, paste0("...../7.txt"))
  summary_string = paste0(summary_string, paste(r_str,collapse = ", "),"\n\n")
  
  if (calc_clearance) {
    write_file(u, paste0("...../71.txt"))
    clearance_date = Findzero(mcmc_pars_estimate, init_sets_list)
    write_file(u, paste0("...../72.txt"))
    summary_string = paste0(summary_string, paste(names(clearance_date), collapse = ", "))
    write_file(u, paste0("...../73.txt"))
    summary_string = paste0(summary_string, "\n", paste(clearance_date, collapse = ", "), "\n")
  }
  
  write_file(summary_string, paste0("...../summary_run_",run_id,".txt"))
  write_file(u, paste0("...../8.txt"))
  # cairo_pdf(paste0("...../par_cor_run_",run_id,".pdf"),width=10,height=10)
  # correlationPlot_modified(mcmc_pars_estimate, scaleCorText = F)
  # dev.off()
  write_file(u, paste0("...../9.txt"))
  png(paste0("...../par_hist_run_",run_id,".png"))
  par(mfrow = c(2, 4))
  for(i in 1:n_pars) {
    hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
    rm(i)
  }
  dev.off()
  write_file(u, paste0("...../10.txt"))
  png(paste0("...../par_traj_run_",run_id,".png"), width=1000, height=500)
  par(mfrow = c(2, 4))
  for(i in 1:n_pars) {
    plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
    rm(i)
  }
  dev.off()
  write_file(u, paste0("...../11.txt"))
  if (plot_combined_fig) {
    SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim)
  }
  
  par(mfrow = c(1, 1))
  # corrplot(cor(mcmc_pars_estimate))
  # pairs(mcmc_pars_estimate)
  
}
##SEIR Deterministic-----------------
SEIRpred <- function(pars, init_settings) {
  tmp_ret=init_settings$var_trans_fun(pars)
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  De <- init_settings$De
  Dq_vec <- init_settings$Dq
  alpha <- init_settings$alpha
  Dh <- init_settings$Dh
  N <- init_settings$N
  flowN_vec <- init_settings$flowN
  init_states <- init_settings$init_states
  days_to_fit <- init_settings$days_to_fit
  ## ODE function based on deterministic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    Dq <- stage_pars[3]
    n <- stage_pars[4]
    ## old states number: c(S, E, P, I, A, H, R)
    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    I <- states_old[4]
    A <- states_old[5]
    H <- states_old[6]
    R <- states_old[7]
    ## new values
    S_new <- S - b * S * (alpha * P + I + alpha * A) / N + n - n * S / N
    E_new <- E + b * S * (alpha * P + I + alpha * A) / N - E / De - n * E / N
    P_new <- P +  E / De  - P / Dp - n * P / N
    I_new <- I + r * P / Dp - I / Di - I / Dq
    A_new <- A + (1 - r) * P / Dp - A / Di - n * A / N
    H_new <- H + I / Dq - H / Dh
    R_new <- R + H / Dh + (A + I) / Di - n * R / N
    Onset_expect <- r * P / Dp
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
  
  myold_states <- init_states
  
  for (i_stage in 1:n_stage) {
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  return(states_mat)
}
##Estimating R0 for the five period----------------------
estimate_R <- function(pars, init_settings) {
  tmp_ret=init_settings$var_trans_fun(pars)
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  ##
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  Dq_vec <- init_settings$Dq
  alpha <- init_settings$alpha
  N <- init_settings$N
  flowN_vec <- init_settings$flowN
  #
  R0_est <- rep(NA, n_stage)
  for(i in 1:n_stage) {
    b <- b_vec[i]
    r <- r_vec[i]
    Dq <- Dq_vec[i]
    n <- flowN_vec[i]
    R0_est[i] <- alpha * b / (1 / Dp + n / N) + (1 - r) * alpha * b / (1 / Di + n / N) + r * b / (1 / Di + 1 / Dq)
    rm(i, b, r, Dq, n)
  }
  return(R0_est)
}
##Find the date when 0 case occurs
Findzero <- function(pars_estimate, init_settings) {
  idx_to_date <- function(d) {
    d <- as.numeric(d)
    mydate <- paste(rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"),c(31, 29, 31, 30, 31, 30, 31, 31, 30)),
                    c(1:31, 1:29, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30))
    return(mydate[d])
  }
  ##
  init_settings$days_to_fit <- 1:200
  ##
  est_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, c("I", "P", "E", "A")])
  est_mat_I <- est_mat[1:200, ]
  est_mat_EPIA <- est_mat[1:200, ] + est_mat[201:400, ] + est_mat[401:600, ] + est_mat[601:800, ]
  ## I = 0
  est <- apply(est_mat_I, 2, function(x) which(x == 0)[1])
  lowdate <- idx_to_date(round(quantile(est, 0.025), 0))
  meandate <- idx_to_date(round(mean(est), 0))
  upperdate <- idx_to_date(round(quantile(est, 0.975), 0))
  intdate1 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  rm(est)
  # E + P + I+ A = 0
  est <- apply(est_mat_EPIA, 2, function(x) tail(which(x != 0), 1) + 1)
  lowdate <- idx_to_date(round(quantile(est, 0.025), 0))
  meandate <- idx_to_date(round(mean(est), 0))
  upperdate <- idx_to_date(round(quantile(est, 0.975), 0))
  intdate2 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  rm(est)
  #
  intdate <- c(intdate1, intdate2)
  names(intdate) <- c("I=0", "E+P+I+A=0")
  return(intdate)
}
## Stochastic SEIR model--------------
SEIRsimu <- function(pars, init_settings, num_periods = 5) {
  tmp_ret=init_settings$var_trans_fun(pars)
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
 
  ##
  startp1<-1
  endp1<-23
  
  startp2<- 24
  endp2<-100
  
  startp3<-101
  endp3<-195
  
  startp4<-196
  endp4<-241
  
  startp5<-242
  endp5<-255
  
  
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  De <- init_settings$De
  Dq_vec <- init_settings$Dq
  alpha <- init_settings$alpha
  Dh <- init_settings$Dh
  N <- init_settings$N
  flowN_vec <- init_settings$flowN
  init_states <- init_settings$init_states
  days_to_fit <- init_settings$days_to_fit
  ## ODE function based on stochastic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    b <- stage_pars[1]
    r <- stage_pars[2]
    Dq <- stage_pars[3]
    #n <- stage_pars[4]
    n <- rpois(1, lambda = stage_pars[4])      ## stochastic, Poisson Distribution
    ## old states number: c(S, E, P, I, A, H, R)
    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    I <- states_old[4]
    A <- states_old[5]
    H <- states_old[6]
    R <- states_old[7]
    ## S
    ## meaning S->E, S->, S->S
    pS_vec <- c(b * (alpha * P + I + alpha * A) / N, n / N, 1 - b * (alpha * P + I + alpha * A) / N - n / N)
    sample_S <- rmultinom(1, size = S, prob = pS_vec)
    ## E
    ## meaning E->P, E->, E->E
    pE_vec <- c(1 / De, n / N, 1 - 1 / De - n / N)
    sample_E <- rmultinom(1, size = E, prob = pE_vec)
    ## P
    ## meaning P->I, P->A, P->, P->P
    pP_vec <- c(r / Dp, (1 - r) / Dp, n/N, 1 - 1 / Dp - n/N)
    sample_P <- rmultinom(1, size = P, prob = pP_vec)
    ## I
    ## meaning I->H, I->R, I->I
    pI_vec <- c(1 / Dq, 1 / Di, 1 - 1 / Dq - 1 / Di)
    sample_I <- rmultinom(1, size = I, prob = pI_vec)
    ## A
    ## meaning A->R, A->, A->A
    pA_vec <- c(1 / Di, n / N, 1 - 1 / Di - n / N)
    sample_A <- rmultinom(1, size = A, prob = pA_vec)
    ## H
    ## meaning H->R, H->H
    pH_vec <- c(1 / Dh, 1 - 1 / Dh)
    sample_H <- rmultinom(1, size = H, prob = pH_vec)
    ## R
    ## meaning R->, R->R
    pR_vec <- c(n / N, 1 - n / N)
    sample_R <- rmultinom(1, size = R, prob = pR_vec)
    ## new values
    S_new <- sample_S[3] + n
    E_new <- sample_E[3] + sample_S[1]
    P_new <- sample_P[4] + sample_E[1]
    I_new <- sample_I[3] + sample_P[1]
    A_new <- sample_A[3] + sample_P[2]
    H_new <- sample_H[2] + sample_I[1]
    R_new <- sample_R[2] + sample_I[2] + sample_A[1] + sample_H[1]
    Onset_expect <- sample_P[1]
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
  ## evovle the system according to the discretized ODEs
  stage_start <- c(startp1, startp2, startp3, startp4, startp5)               # corresponding to dates ...
  stage_end <- c(endp1, endp2, endp3,endp4,length(days_to_fit))    # corresponding to dates Jan9, Jan22, Feb1, Feb16, the last day (could vary)
  ##
  myold_states <- init_states
  for (i_stage in 1:5) {
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_start[i_stage]:stage_end[i_stage]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  if(num_periods == 5) { 
    states_mat <- states_mat
  }
  
  else if (num_periods %in% c(2,4,3)) {
    i_stage <- num_periods
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in stage_start[i_stage]:length(days_to_fit)) {
      myold_states <- states_mat[d - 1, -1]
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
    }
  }
  else if(num_periods ==1)
  {
    i_stage <- num_periods
    stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[i_stage], Dq = Dq_vec[i_stage], n = flowN_vec[i_stage])
    for (d in 2:length(days_to_fit)) {
      myold_states <- states_mat[d - 1, -1]
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
    }
  }
  
  return(states_mat)
}
##plot SEIR
SEIRplot <- function(pars_estimate, file_name, init_settings, panel_B_R_ylim=5) {
  

  
  ##
  startp1<-1
  endp1<-23
  
  startp2<- 24
  endp2<-100
  
  startp3<-101
  endp3<-195
  
  startp4<-196
  endp4<-241
  
  startp5<-242
  endp5<-255
  
  
  ##
  init_settings$days_to_fit <- 1:257
  library(vioplot)
  ##
  ##
  ##
  # onset_obs_all <- init_settings$daily_new_case_all
  # ptime <- 1:length(onset_obs_all)
  onset_obs <- init_settings$daily_new_case
  ptime <- 1:length(onset_obs)
  ##
  ##
  ##
  # mydate <- c(paste("Feb", 25:29), paste("Mar", 1:31), paste("April", 1:30),paste("May", 1:31), paste("Jun", 1:30),paste("Jul", 1:31),paste("Aug", 1:31),paste("Sep", 1:30),paste("Oct", 1:25))
  mydate <- c(paste("Mar", 4:31), paste("April", 1:30),paste("May", 1:31), paste("Jun", 1:30),paste("Jul", 1:31),paste("Aug", 1:31),paste("Sep", 1:30),paste("Oct", 1:25))
  
  #
  pdf(paste0("...../Figure_", file_name, ".pdf"), width = 9, height = 10)
  par(mar = c(4, 5, 3, 1))
  layout(matrix(c(1:6), byrow = T, nrow = 3))
  
  ##########################################   Panel A  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start A
  # plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.05), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  plot(ptime, estN_mean, xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 260,60), labels = mydate[seq(1, 260,60)])
  #
  abline(v = c(startp2,startp3,startp4,startp5 ,255), lty = 3, lwd = 2, col = "darkgrey")
  text(c(startp2,startp3,startp4,startp5, 255), par()$usr[4], labels = mydate[c(startp2,startp3,startp4,startp5, 255)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:255], rev(ptime[1:255])), c(estN_up[1:255], rev(estN_low[1:255])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:254)], rev(ptime[-c(1:254)])), c(estN_up[-c(1:254)], rev(estN_low[-c(1:254)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:254], estN_mean[1:254], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:254)], estN_mean[-c(1:254)], col = "#0072B5FF", pch = 17, cex = 0.8)
  ##
  ##
  ##
  # points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  points(ptime, onset_obs, col = "black", pch = 4, cex = 0.8)
  ##
  ##
  ##
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)
  
  ##########################################   Panel B  ##########################################################
  estRt_mat <- apply(pars_estimate, 1, function(x) estimate_R(pars = x, init_settings = init_settings))
  estRt_mat <- t(estRt_mat)
  ##
  rt_mean <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) mean(x)), 2))
  rt_low <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) quantile(x, 0.025)), 2))
  rt_up <- sprintf("%.2f", round(apply(estRt_mat, 2, function(x) quantile(x, 0.975)), 2))
  #
  #vioplot(estRt_mat[, 1], estRt_mat[, 2], estRt_mat[, 3], estRt_mat[, 4], estRt_mat[, 5], names = NA, ylim = c(-0.5, panel_B_R_ylim), col = c("#BC3C29FF","#0072B5FF", "#E18727FF", "#7876B1FF", "#FFDC91FF"), ylab = "", xlab = "")
  vioplot(estRt_mat[, 1], estRt_mat[, 2], estRt_mat[, 3], estRt_mat[, 4], estRt_mat[, 5], names = NA, ylim = c(-1.5, panel_B_R_ylim),col = c("#BC3C29FF","#0072B5FF", "#E18727FF", "#7876B1FF", "#FFDC91FF"), ylab = "", xlab = "")
  mtext("Outbreak period (2020)", side = 1, line  = 3, cex = 1.01)
  mtext(expression("R"["0"]), side = 2, line = 3, cex = 1.01)
  axis(1, at = c(1, 2, 3, 4,5), tick = F, labels = c(paste0(mydate[c(1)],"-",mydate[c(endp1)]), 
                                                     paste0(mydate[c(startp2)],"-",mydate[c(endp2)]),
                                                     paste0(mydate[c(startp3)],"-",mydate[c(endp3)]),
                                                     paste0(mydate[c(startp4)],"-",mydate[c(endp4)]),
                                                     paste0(mydate[c(startp5)],"-",mydate[c(endp5)])))
  abline(h = 1, lwd = 2, lty = 3, col = "red")
  #
  text(1, min(estRt_mat[, 1]) - 0.2, labels = rt_mean[1])
  text(1, min(estRt_mat[, 1]) - 0.45, labels = paste("(", rt_low[1], "-", rt_up[1], ")", sep = ""))
  text(2, min(estRt_mat[, 2]) - 0.2, labels = rt_mean[2])
  text(2, min(estRt_mat[, 2]) - 0.45, labels = paste("(", rt_low[2], "-", rt_up[2], ")", sep = ""))
  text(3, max(estRt_mat[, 3]) + 0.4, labels = rt_mean[3])
  text(3, max(estRt_mat[, 3]) + 0.15, labels = paste("(", rt_low[3], "-", rt_up[3], ")", sep = ""))
  text(4, max(estRt_mat[, 4]) + 0.4, labels = rt_mean[4])
  text(4, max(estRt_mat[, 4]) + 0.15, labels = paste("(", rt_low[4], "-", rt_up[4], ")", sep = ""))
  text(5, max(estRt_mat[, 5]) + 0.4, labels = rt_mean[5])
  text(5, max(estRt_mat[, 5]) + 0.15, labels = paste("(", rt_low[5], "-", rt_up[5], ")", sep = ""))
  
  
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "B", xpd = T, cex = 2)
  
  ##########################################   Panel C  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 4)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start C
  plot(ptime, estN_mean, xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  # plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.05), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 260,60), labels = mydate[seq(1, 260,60)])
  #
  abline(v = c(startp2,startp3,startp4,startp5, 255), lty = 3, lwd = 2, col = "darkgrey")
  text(c(startp2,startp3,startp4,startp5, 255), par()$usr[4], labels = mydate[c(startp2,startp3,startp4,startp5, 255)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:startp5], rev(ptime[1:startp5])), c(estN_up[1:startp5], rev(estN_low[1:startp5])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:endp4)], rev(ptime[-c(1:endp4)])), c(estN_up[-c(1:endp4)], rev(estN_low[-c(1:endp4)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:endp4], estN_mean[1:endp4], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:endp4)], estN_mean[-c(1:endp4)], col = "#0072B5FF", pch = 17, cex = 0.8)
  ##
  ##
  ##
  # points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  points(ptime, onset_obs, col = "black", pch = 4, cex = 0.8)
  ##
  ##
  ##
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "C", xpd = T, cex = 2)
  
  ##########################################   Panel D  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 3)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start D
  #  plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.01), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  plot(ptime, estN_mean, xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 260,60), labels = mydate[seq(1, 260,60)])
  #
  abline(v = c(startp2,startp3,startp4,startp5, 255), lty = 3, lwd = 2, col = "darkgrey")
  text(c(startp2,startp3,startp4,startp5,255), par()$usr[4], labels = mydate[c(startp2,startp3,startp4,startp5, 255)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:startp4], rev(ptime[1:startp4])), c(estN_up[1:startp4], rev(estN_low[1:startp4])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:endp3)], rev(ptime[-c(1:endp3)])), c(estN_up[-c(1:endp3)], rev(estN_low[-c(1:endp3)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:endp4], estN_mean[1:endp4], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:endp4)], estN_mean[-c(1:endp4)], col = "#0072B5FF", pch = 17, cex = 0.8)
  ##
  ##
  ##
  # points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  points(ptime, onset_obs, col = "black", pch = 4, cex = 0.8)
  ##
  ##
  ##
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "D", xpd = T, cex = 2)
  
  ##########################################   Panel E  ##########################################################
  estN_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 2)[, "Onset_expect"])
  estN_mean <- round(apply(estN_mat, 1, mean), 0)
  estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
  estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
  # start E
  plot(ptime, estN_mean, xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  #plot(ptime, estN_mean, ylim = c(0, max(estN_up, onset_obs_all) * 1.01), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = seq(1, 260,60), labels = mydate[seq(1, 260,60)])
  #
  abline(v = c(startp2,startp3,startp4,startp5, 255), lty = 3, lwd = 2, col = "darkgrey")
  text(c(startp2,startp3,startp4,startp5,255), par()$usr[4], labels = mydate[c(startp2,startp3,startp4,startp5, 255)], col = "darkgrey", pos = 3, xpd = T)
  #
  polygon(c(ptime[1:startp3], rev(ptime[1:startp3])), c(estN_up[1:startp3], rev(estN_low[1:startp3])), col = "#F39B7FB2", border = NA)
  polygon(c(ptime[-c(1:endp2)], rev(ptime[-c(1:endp2)])), c(estN_up[-c(1:endp2)], rev(estN_low[-c(1:endp2)])), col = "#4DBBD5B2", border = NA)
  #
  points(ptime[1:endp2], estN_mean[1:endp2], col = "#BC3C29FF", pch = 16, cex = 0.8)
  points(ptime[-c(1:endp2)], estN_mean[-c(1:endp2)], col = "#0072B5FF", pch = 14, cex = 0.8)
  ##
  ##
  ##
  # points(ptime, onset_obs_all, col = "black", pch = 4, cex = 0.8)
  points(ptime, onset_obs, col = "black", pch = 4, cex = 0.8)
  ##
  ##
  ##
  #
  legend("topleft", legend = c("Observed", "Fitted",  "Predicted"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "E", xpd = T, cex = 2)
  
  ##########################################   Panel F  ##########################################################
  estAIP_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 5)[, c("I", "A", "P")])
  estI_mat <- estAIP_mat[ptime, ]
  estA_mat <- estAIP_mat[ptime + length(ptime), ]
  estP_mat <- estAIP_mat[ptime + length(ptime) * 2, ]
  estI_mean <- apply(estI_mat, 1, mean)
  estA_mean <- apply(estA_mat, 1, mean)
  estP_mean <- apply(estP_mat, 1, mean)
  estAIP_dat <- rbind(estI_mean, estA_mean, estP_mean)
  barpos <- barplot(estAIP_dat, col = c("#BC3C29FF", "#FFDC91FF", "#0072B5FF"), xlab = "", ylab = "", border = "NA")
  mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("No. of active infectious cases", side = 2, line = 3, cex = 1.01)
  axis(1, at = barpos[seq(1, 260,60)], labels = mydate[seq(1, 260,60)])
  legend("topleft", legend = c("Presymptomatic (P)", "Unascertained (A)", "Ascertained (I)"), fill = c("#0072B5FF", "#FFDC91FF", "#BC3C29FF"), bty = "n")
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "F", xpd = T, cex = 2)
  ##figure_F finished
  dev.off()
}



#run the simulation------------------------------------------------
init_sets_list=get_init_sets_list(r0 = 0.5)

#with starting values for the MCMC to get the same displayed graph
#0.716	0.242	0.46	0.539	0.414	0.902	0.567
#	0.172	0.437	1.137	-1.835	-1.519	-0.5	-0.978

startv=c(0.716,0.242,0.46,0.539,0.414,0.902,0.567,0.172,0.437,1.137,-1.835,-1.519,-0.5,-0.97)
startValue<-matrix(NA, 1,14)
for (i in 1:12){startValue[,i]<-startv[i]}

SEIRfitting(init_sets_list, randomize_startValue = F,startValue = startValue, run_id = "main_analysis_with_starting_values", output_ret = T, skip_MCMC=F)

#randomn starting values for the MCMC
ll=SEIRfitting(init_sets_list, randomize_startValue = T,
               run_id = "main_analysis", output_ret = T, skip_MCMC=F)


SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep1", output_ret = T, skip_MCMC=F)
SEIRfitting(init_sets_list, randomize_startValue = T, run_id = "main_analysis_rep2", output_ret = T, skip_MCMC=F)