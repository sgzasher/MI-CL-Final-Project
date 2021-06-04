# Preamble ----------------------------------------------------------------
setwd("C:/Users/s-mas/OneDrive/Documents/Stanford GSB/Code and Data/ML_CI Project")
library(data.table)
library(fixest)
library(dplyr)
library(Hmisc)
library(grf)
library(rpart)
library(ggplot2)
library(neuralnet)
library(glmnet)
library(gridExtra)
library(ivdesc)
library(rlang)
library(groupdata2)
library(psych)
library(stringr)
set.seed(080699)

#######
# In places there are commented out code followed by a data read, 
# this is because some code runs on the Sherlock HPC cluster, 
# and I pull the data back onto my local machine for graphics processing.
#######

# Data --------------------------------------------------------------------
data <- fread("./data/analysis_data_owners.csv")
pscore.merge <- fread("./data/pscore_output.csv")

pscore.merge <- select(pscore.merge,
                       -starts_with("dec."),
                       -starts_with("imbens."),
                       -V1)

data <- left_join(data,
                  pscore.merge,
                  by = c("pid", "casenumber"))


# Function: Balanced Panels -----------------------------------------------
# Need to define Imbens-Rubin recursive function from pscores.R
balance_panel <- function(data, cohort.selection, outcome){
  
  # Output is data of correct cohort, balanced across panel
  # With the 1) Treatment 2) IDs 3) PScore data 4) Outcome
  ## This involves a missing-at-random assumption that is
  ## hard (read: impossible) to justify!
  
  # Read in with relevant cohort & variables
  study.data <- data[data$year_beg_fs==cohort.selection,]
  study.data <- select(study.data,
                       year_orig,
                       year_beg_fs,
                       pid,
                       casenumber,
                       evt,
                       fc_al3,
                       all_of(outcome))
  
  # 9 year panel
  study.data <- subset(study.data, evt <= 4 & evt >= -4)
  
  # With Pscores
  study.data <- left_join(study.data,
                          pscore.merge,
                          by = c("pid", "casenumber"))
  
  # Without bad propensity score observations. 
  study.data <- subset(study.data,
                       (e.hat.rf <= 0.9 & e.hat.rf >= 0.1) & 
                         (e.hat.DGT <= 0.9 & e.hat.DGT >= 0.1) & 
                         (e.hat.postlasso <= 0.9 & e.hat.postlasso >= 0.1) & 
                         (e.hat.tree <= 0.9 & e.hat.tree >= 0.1))
  
  # Balanced
  study.data <- study.data %>%
                subset(complete.cases(study.data)) %>%
                add_count(pid) %>%
                subset(n==9) %>%
                select(-n) %>%
                ungroup()
  
  # Re-apply imbens for balanced panels
  imbens.rf <- imbens_recurse(study.data, study.data$e.hat.rf, "fc_al3")
  imbens.postlasso <- imbens_recurse(study.data, study.data$e.hat.postlasso, "fc_al3")
  imbens.DGT <- imbens_recurse(study.data, study.data$e.hat.DGT, "fc_al3")
  imbens.tree <- imbens_recurse(study.data, study.data$e.hat.tree, "fc_al3")
  
  study.data$imbens.e.hat.rf <- as.numeric(cut(study.data$e.hat.rf,
                                breaks = imbens.rf[[2]],
                                include.lowest = TRUE))
  study.data$imbens.e.hat.postlasso <- as.numeric(cut(study.data$e.hat.postlasso,
                                               breaks = imbens.postlasso[[2]],
                                               include.lowest = TRUE))
  study.data$imbens.e.hat.DGT <- as.numeric(cut(study.data$e.hat.DGT,
                                               breaks = imbens.DGT[[2]],
                                                include.lowest = TRUE))
  study.data$imbens.e.hat.tree <- as.numeric(cut(study.data$e.hat.tree,
                                               breaks = imbens.tree[[2]],
                                               include.lowest = TRUE))
  

  return(study.data)
}

# Stratification: Estimation Function -------------------------------------------

strat_est <- function(data.cohort, event.year, outcome.var, pscore){
  
  # Define some variables
  bin.selection <- paste0("imbens.", pscore)
  year.selection <- event.year
  
  # Get bin info
  study.data <- data.cohort
  
  study.data$bins <- select(study.data, 
                            all_of(bin.selection))

  study.data$outcome <- select(study.data,
                               all_of(outcome.var))
  

  # Panel-ify the outcome variable
  baseline <- mean(study.data[evt==-1,]$outcome)
  
  study.data$outcome <- study.data %>%
    mutate(outcome = outcome - baseline) %>%
    select(outcome)
  
  # Relevant Data
  study.data <- study.data[study.data$evt==event.year,]
  study.data <- select(study.data,
                       outcome,
                       bins,
                       fc_al3,
                       casenumber)
  
  loop.max <- max(study.data$bins)

  # Loop over propensity score bins
  # Get ATE within bin, N, variance
  
  ate.mat <- rep(NA, loop.max)
  n.mat <- rep(NA, loop.max)
  j.mat <- rep(NA, loop.max)
  var.mat <- rep(NA, loop.max)
  
  for(i in 1:loop.max){
    

    # Stratum details
    control.data <- study.data[study.data$fc_al3==0 & study.data$bins == i,]
    treated.data <- study.data[study.data$fc_al3==1 & study.data$bins == i,]
    
    n.c <- nrow(control.data)
    n.t <- nrow(treated.data)
    
    j.c <- length(unique(control.data$casenumber))
    j.t <- length(unique(treated.data$casenumber))
  
    m.c <- mean(control.data$outcome)
    m.t <- mean(treated.data$outcome)

    # ATE within stratum
    ate.mat[i] <- m.t - m.c
    
    # N within stratum
    n.mat[i] <- n.t + n.c
    
    # J within stratum
    j.mat[i] <- j.c + j.t
    
    # Variance within stratum (See Imbens & Rubin 2015)
    # Cluster robust version (J instead of N)
    
    # Treated
    treated.data <-  treated.data %>%
                     mutate(scores = outcome) %>%
                     group_by(casenumber) %>%
                     summarize(tauj=mean(scores)) %>%
                     mutate (tau=mean(tauj)) %>%
                     mutate (diff=(tau-tauj)^2) %>%
                     ungroup()
    
    
    unique.cases.t <- as.data.frame(unique(treated.data$casenumber))
    colnames(unique.cases.t) <- c("casenumber")
    
    unique.cases.t <- semi_join(unique.cases.t,
                                select(treated.data,
                                       casenumber,
                                       diff),
                              by = "casenumber")
    
    t.se.part <- sqrt(1/(j.t*(j.t-1)) * sum(unique.cases.t$diff, na.rm = TRUE))
    
    # Control
    control.data <-  control.data %>%
                      mutate(scores = outcome) %>%
                      group_by(casenumber) %>%
                      summarize(tauj=mean(scores)) %>%
                      mutate (tau=mean(tauj)) %>%
                      mutate (diff=(tau-tauj)^2) %>%
                      ungroup()
    
    
    unique.cases.c <- as.data.frame(unique(control.data$casenumber))
    colnames(unique.cases.c) <- c("casenumber")
    
    unique.cases.c <- semi_join(unique.cases.c,
                                select(treated.data,
                                       casenumber,
                                       diff),
                                by = "casenumber")
    
    c.se.part <- sqrt(1/(j.c*(j.c-1)) * sum(unique.cases.c$diff, na.rm = TRUE))
    
    # Combine
    var.mat[i] <- t.se.part + c.se.part
  }
  
  ## Combine estimates
  ## Aware that sometimes there is an Nt or Nc=1 in a bin and we get
  ## Returned an empty variance; I don't have time for a better fix
  ## Than this right now. Dropping obs. non-random, but not dropping
  ## Many at all. 
  idx <- as.numeric(is.na(var.mat))
  n.total <- sum(n.mat[idx==0])
  j.total <- sum(j.mat[idx==0])
  ate.est <- sum(ate.mat[idx==0] * (n.mat[idx==0]/n.total))
  var.est <- sum(var.mat[idx==0] * (j.mat[idx==0]/j.total)^2)
  
  # Gaussian CI
  ci.upper <- ate.est + 1.96 * sqrt(var.est)
  ci.lower <- ate.est - 1.96 * sqrt(var.est)
  
  t.est <- ate.est/sqrt(var.est)
  
  ate.details <- c(N=n.total, 
                   J=j.total, 
                   ATE=ate.est, 
                   VAR=var.est, 
                   CIU=ci.upper, 
                   CIL=ci.lower, 
                   T=t.est)
  return(ate.details)
}


# Function: Pooled ATT Stratification -------------------------------------
pooled_strat_att <- function(data.prefix, event.year, outcome){
  
    data.2009 <- get(paste0(data.prefix, ".2009"))
    data.2010 <- get(paste0(data.prefix, ".2010"))
    data.2011 <- get(paste0(data.prefix, ".2011"))
    data.2012 <- get(paste0(data.prefix, ".2012"))
    
    # RF
    rf.2009 <- strat_est(data.2009, event.year, outcome, "e.hat.rf")
    rf.2010 <- strat_est(data.2010, event.year, outcome, "e.hat.rf")
    rf.2011 <- strat_est(data.2011, event.year, outcome, "e.hat.rf")
    rf.2012 <- strat_est(data.2012, event.year, outcome, "e.hat.rf")
    
    rf.n <- rf.2009[[1]] + rf.2010[[1]] + rf.2011[[1]] + rf.2012[[1]]
    rf.j <- rf.2009[[2]] + rf.2010[[2]] + rf.2011[[2]] + rf.2012[[2]]
    rf.att <- ((rf.2009[[1]]/rf.n)*rf.2009[[3]]) + 
              ((rf.2010[[1]]/rf.n)*rf.2010[[3]]) +
              ((rf.2011[[1]]/rf.n)*rf.2011[[3]]) +
              ((rf.2012[[1]]/rf.n)*rf.2012[[3]])
    rf.var <- ((rf.2009[[2]]/rf.j)^2*rf.2009[[4]]) + 
              ((rf.2010[[2]]/rf.j)^2*rf.2010[[4]]) +
              ((rf.2011[[2]]/rf.j)^2*rf.2011[[4]]) +
              ((rf.2012[[2]]/rf.j)^2*rf.2012[[4]])
    rf.se <- sqrt(rf.var)
    rf.t <- rf.att/rf.se
    rf.p <- 1.96*(1 - pnorm(abs(rf.t)))
    
    # RT
    rt.2009 <- strat_est(data.2009, event.year, outcome, "e.hat.tree")
    rt.2010 <- strat_est(data.2010, event.year, outcome, "e.hat.tree")
    rt.2011 <- strat_est(data.2011, event.year, outcome, "e.hat.tree")
    rt.2012 <- strat_est(data.2012, event.year, outcome, "e.hat.tree")
    
    rt.n <- rt.2009[[1]] + rt.2010[[1]] + rt.2011[[1]] + rt.2012[[1]]
    rt.j <- rt.2009[[2]] + rt.2010[[2]] + rt.2011[[2]] + rt.2012[[2]]
    rt.att <- ((rt.2009[[1]]/rt.n)*rt.2009[[3]]) + 
              ((rt.2010[[1]]/rt.n)*rt.2010[[3]]) +
              ((rt.2011[[1]]/rt.n)*rt.2011[[3]]) +
              ((rt.2012[[1]]/rt.n)*rt.2012[[3]])
    rt.var <- ((rt.2009[[2]]/rt.j)^2*rt.2009[[4]]) + 
              ((rt.2010[[2]]/rt.j)^2*rt.2010[[4]]) +
              ((rt.2011[[2]]/rt.j)^2*rt.2011[[4]]) +
              ((rt.2012[[2]]/rt.j)^2*rt.2012[[4]])
    rt.se <- sqrt(rt.var)
    rt.t <- rt.att/rt.se
    rt.p <- 1.96*(1 - pnorm(abs(rt.t)))
    
    # DGT
    DGT.2009 <- strat_est(data.2009, event.year, outcome, "e.hat.DGT")
    DGT.2010 <- strat_est(data.2010, event.year, outcome, "e.hat.DGT")
    DGT.2011 <- strat_est(data.2011, event.year, outcome, "e.hat.DGT")
    DGT.2012 <- strat_est(data.2012, event.year, outcome, "e.hat.DGT")
    
    DGT.n <- DGT.2009[[1]] + DGT.2010[[1]] + DGT.2011[[1]] + DGT.2012[[1]]
    DGT.j <- DGT.2009[[2]] + DGT.2010[[2]] + DGT.2011[[2]] + DGT.2012[[2]]
    DGT.att <-((DGT.2009[[1]]/DGT.n)*DGT.2009[[3]]) + 
              ((DGT.2010[[1]]/DGT.n)*DGT.2010[[3]]) +
              ((DGT.2011[[1]]/DGT.n)*DGT.2011[[3]]) +
              ((DGT.2012[[1]]/DGT.n)*DGT.2012[[3]])
    DGT.var <-((DGT.2009[[2]]/DGT.j)^2*DGT.2009[[4]]) + 
              ((DGT.2010[[2]]/DGT.j)^2*DGT.2010[[4]]) +
              ((DGT.2011[[2]]/DGT.j)^2*DGT.2011[[4]]) +
              ((DGT.2012[[2]]/DGT.j)^2*DGT.2012[[4]])
    DGT.se <- sqrt(DGT.var)
    DGT.t <- DGT.att/DGT.se
    DGT.p <- 1.96*(1 - pnorm(abs(DGT.t)))
    
    # Post-LASSO
    pl.2009 <- strat_est(data.2009, event.year, outcome, "e.hat.postlasso")
    pl.2010 <- strat_est(data.2010, event.year, outcome, "e.hat.postlasso")
    pl.2011 <- strat_est(data.2011, event.year, outcome, "e.hat.postlasso")
    pl.2012 <- strat_est(data.2012, event.year, outcome, "e.hat.postlasso")
    
    pl.n <- pl.2009[[1]] + pl.2010[[1]] + pl.2011[[1]] + pl.2012[[1]]
    pl.j <- pl.2009[[2]] + pl.2010[[2]] + pl.2011[[2]] + pl.2012[[2]]
    pl.att <- ((pl.2009[[1]]/pl.n)*pl.2009[[3]]) + 
              ((pl.2010[[1]]/pl.n)*pl.2010[[3]]) +
              ((pl.2011[[1]]/pl.n)*pl.2011[[3]]) +
              ((pl.2012[[1]]/pl.n)*pl.2012[[3]])
    pl.var <- ((pl.2009[[2]]/pl.j)^2*pl.2009[[4]]) + 
              ((pl.2010[[2]]/pl.j)^2*pl.2010[[4]]) +
              ((pl.2011[[2]]/pl.j)^2*pl.2011[[4]]) +
              ((pl.2012[[2]]/pl.j)^2*pl.2012[[4]])
    pl.se <- sqrt(pl.var)
    pl.t <- pl.att/pl.se
    pl.p <- 1.96*(1 - pnorm(abs(pl.t)))
    
    rf.results <- c(rf.n, rf.att, rf.var, rf.se, rf.t, rf.p)
    rt.results <- c(rt.n, rt.att, rt.var, rt.se, rt.t, rt.p)
    DGT.results <- c(DGT.n, DGT.att, DGT.var, DGT.se, DGT.t, DGT.p)
    pl.results <- c(pl.n, pl.att, pl.var, pl.se, pl.t, pl.p)
    
    returnme <- list(forest=rf.results, 
                     tree=rt.results, 
                     dgt=DGT.results, 
                     lasso=pl.results)
    return(returnme)
}
  


# Panel-Ify ---------------------------------------------------------------

panel_me <- function(data.read, outcome){
  
  baseline <- subset(data.read, pull(data.read, evt)==-1)
  baseline <- mean(pull(baseline, outcome), na.rm=TRUE)
  
  data.read[, all_of(outcome)] <- pull(data.read, outcome) - baseline
 
  
  return(data.read)

}

# Stratification: Pooled ATTs ---------------------------------------------

# For: Moved, Unpaid Collections, high school index, divorce
# Calculate ATE in period 4 for all cohorts...
# Then pool results 

# Moved
#moved.2009 <- balance_panel(data, 2009, "moved")
#moved.2010 <- balance_panel(data, 2010, "moved")
#moved.2011 <- balance_panel(data, 2011, "moved")
#moved.2012 <- balance_panel(data, 2012, "moved")
#
#moved.balanced.cohorts <- rbind(moved.2009,
#                                moved.2010,
#                                moved.2011,
#                                moved.2012)
#
#write.csv(moved.balanced.cohorts, "./data/moved_balanced_cohorts.csv")

moved.balanced.cohorts <- fread("./data/moved_balanced_cohorts.csv")
moved.2009 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2009), "moved")
moved.2010 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2010), "moved")
moved.2011 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2011), "moved")
moved.2012 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2012), "moved")

# High School Index
#school.2009 <- balance_panel(data, 2009, "schoolindex_high")
#school.2010 <- balance_panel(data, 2010, "schoolindex_high")
#school.2011 <- balance_panel(data, 2011, "schoolindex_high")
#school.2012 <- balance_panel(data, 2012, "schoolindex_high")
#
#school.balanced.cohorts <- rbind(school.2009,
#                                 school.2010,
#                                 school.2011,
#                                 school.2012)
#
#write.csv(school.balanced.cohorts, "./data/school_balanced_cohorts.csv")

school.balanced.cohorts <- fread("./data/school_balanced_cohorts.csv")
school.2009 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2009), "schoolindex_high")
school.2010 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2010), "schoolindex_high")
school.2011 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2011), "schoolindex_high")
school.2012 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2012), "schoolindex_high")


# Divorce
#divorce.2009 <- balance_panel(data, 2009, "cumulnum_divorce")
#divorce.2010 <- balance_panel(data, 2010, "cumulnum_divorce")
#divorce.2011 <- balance_panel(data, 2011, "cumulnum_divorce")
#divorce.2012 <- balance_panel(data, 2012, "cumulnum_divorce")
#
#divorce.balanced.cohorts <- rbind(divorce.2009,
#                                  divorce.2010,
#                                  divorce.2011,
#                                  divorce.2012)
#
#write.csv(divorce.balanced.cohorts, "./data/divorce.balanced.cohorts")

divorce.balanced.cohorts <- fread("./data/divorce.balanced.cohorts")
divorce.2009 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2009), "cumulnum_divorce")
divorce.2010 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2010), "cumulnum_divorce")
divorce.2011 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2011), "cumulnum_divorce")
divorce.2012 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2012), "cumulnum_divorce")

# Unpaid Collections
#unpaid.2009 <- balance_panel(data, 2009, "num00_collection_unpaid")
#unpaid.2010 <- balance_panel(data, 2010, "num00_collection_unpaid")
#unpaid.2011 <- balance_panel(data, 2011, "num00_collection_unpaid")
#unpaid.2012 <- balance_panel(data, 2012, "num00_collection_unpaid")
#
#unpaid.balanced.cohorts <- rbind(unpaid.2009,
#                                 unpaid.2010,
#                                 unpaid.2011,
#                                 unpaid.2012)
#
#write.csv(unpaid.balanced.cohorts, "./data/unpaid.balanced.cohorts")

unpaid.balanced.cohorts <- fread("./data/unpaid.balanced.cohorts")
unpaid.2009 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2009), "num00_collection_unpaid")
unpaid.2010 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2010), "num00_collection_unpaid")
unpaid.2011 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2011), "num00_collection_unpaid")
unpaid.2012 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2012), "num00_collection_unpaid")

# ATT Estimates
att.strat.moved <- pooled_strat_att("moved", 4, "moved")
att.strat.school <- pooled_strat_att("school", 4, "schoolindex_high")
att.strat.divorce <- pooled_strat_att("divorce", 4, "cumulnum_divorce")
att.strat.unpaid <- pooled_strat_att("unpaid", 4, "num00_collection_unpaid")

# AIPW: Nuisance Function Estimation --------------------------------------

# Relying on the already-estimated propensity scores for the propensity 
# score nuisance function estimation. Left then is the mu function. So
# here we need a function that trains a mu() for a given observation on
# some given data. 

covariates.prediction <- c("moved", "isowner", "log_finishedsqft_dq", "zipinc_log",
                           "schoolindex_elem", "schoolindex_middle", "schoolindex_high",
                           "cumulnum_divorce", "cumulnum_crime", "cumulnum_dui", 
                           "cumulnum_bankruptcy", "vantage", "num00_foreclosure", 
                           "num00_collection_unpaid", "num00_trade_auto", "num12_trade30", 
                           "num00_mortgage_loanmod", "num00_mortgage_open", 
                           "rtot00_mortgage_open_bal", "rtot00_mortgage_monthpay")

# Pass in e.g. moved.2009 and get the outcome model for the event year you 
# prefer. This is for the forest, tree, linear AIPW implementations. 
mu_est <- function(data.choice, outcome.selection, event.year, method){
  
  covariates <- covariates.prediction[covariates.prediction!=all_of(outcome.selection)]
  
  data.study <- left_join(subset(data.choice, evt==event.year), 
                          select(data, 
                                casenumber,
                                pid, 
                                year_orig,
                                all_of(covariates),
                                zip_year,
                                date_beg_fs),
                         by = c("casenumber", "pid", "year_orig"))

  
  # Need to spot if any variables are missing for all of the dataset
  check <- sapply(select(data.study, all_of(covariates)), function(x) sum(is.na(x)))
  index <- as.numeric(check==nrow(data.study))
  
  covariates <- covariates[index==0]
  
  #train.fmla <- formula(paste(outcome.selection, "~ fc_al3 + ",
  #                            paste(covariates, collapse="+"),
  #                            "+", paste0("(", covariates, "*", "fc_al3", ")", collapse = "+")))

  train.fmla <- formula(paste(outcome.selection, "~ fc_al3 + ",
                              paste(covariates, collapse="+")))
  
  # LASSO linear
  if(method == "lasso"){
    XX <- model.matrix.lm(train.fmla, data.study)[,-1]
    Y <- pull(data.study[complete.cases(select(data.study,
                                               all_of(covariates))),], 
              all_of(outcome.selection))

    data.study <- data.study %>%
                  add_count(casenumber) %>%
                  mutate(weight=1/n)
                  
    lasso <- cv.glmnet(x=XX, 
                       y=Y, 
                       weights=data.study[complete.cases(select(data.study,
                                                                all_of(covariates))),]$weight,  
                       alpha=1.) 
    
    covariates.postlasso <- rownames(coef(lasso, 
                                          s = 'lambda.1se'))[coef(lasso, s = 'lambda.1se')[,1]!= 0][-1] 
    
    if(length(covariates.postlasso) > 0){
      fmla.postlasso <- formula(paste(outcome.selection, "~", 
                                      paste(covariates.postlasso, collapse="+"),
                                      " | zip_year+date_beg_fs"))
    }
    
    else if(length(covariates.postlasso) == 0){
      fmla.postlasso <- formula(paste(outcome.selection, "~ 0", 
                                      " | zip_year+date_beg_fs"))
    }
    
      postlasso <- feols(fmla.postlasso, 
                       data.study, 
                       weight=data.study$weight, 
                       cluster=~casenumber)
    
    
    predict.0 <- data.study %>%
                 mutate(fc_al3 = 0)
    
    predict.1 <- data.study %>%
                 mutate(fc_al3 = 1)
    
    predictions.1 <- predict(postlasso, predict.1)
    predictions.0 <- predict(postlasso, predict.0)
    
    output <- select(data.study,
                     casenumber,
                     pid,
                     starts_with("e.hat."),
                     all_of(outcome.selection),
                     fc_al3)
    
    output$mu1 <- predictions.1
    output$mu0 <- predictions.0
    
    }
  
  # Random Tree
  else if(method == "tree"){
    
    data.study$casenumber <- as.factor(data.study$casenumber)
    data.study <- ungroup(data.study)
    data.study <- fold(data.study,
              k = 5,
              id_col="casenumber")
    
    n <- nrow(data.study)
    predictions.0 <- rep(NA, n)
    predictions.1 <- rep(NA, n)
    
    # Cross-fit mu, folding by case 
    for (idx in 1:5){
      unpruned.tree <- rpart(train.fmla, 
                             data = data.study[data.study$.folds!=idx,],
                             cp=0,
                             method="anova")
      
      cp.idx <- which.min(unpruned.tree$cptable[,"xerror"])
      cp.best <- unpruned.tree$cptable[cp.idx, "CP"]
      propensity.model <- prune(unpruned.tree, cp=cp.best)
      
      predict.0 <- data.study %>%
        mutate(fc_al3 = 0)
      
      predict.1 <- data.study %>%
        mutate(fc_al3 = 1)
      
      predictions.0[data.study$.folds==idx] <- predict(propensity.model, 
                                                      newdata=predict.0[data.study$.folds==idx,], 
                                                      type="vector")
      
      predictions.1[data.study$.folds==idx] <- predict(propensity.model, 
                                                       newdata=predict.1[data.study$.folds==idx,], 
                                                       type="vector")
    }
    
    output <- select(data.study,
                           casenumber,
                           pid,
                           starts_with("e.hat."),
                           all_of(outcome.selection),
                           fc_al3)
     
    output$mu1 <- predictions.1
    output$mu0 <- predictions.0
    
  }
  
  
  # Random Forest (Only 500 trees for computation time)
  else if(method == "forest"){
    n <- nrow(data.study)
    predictions.0 <- rep(NA, n)
    predictions.1 <- rep(NA, n)
    

    XX <- model.matrix.lm(train.fmla, 
                          data.study,
                          na.action = NULL)
    
    outcome.vec <- pull(data.study, 
                           all_of(outcome.selection))
    
    clusters.vec <- as.factor(pull(data.study, 
                                   casenumber))
    
    rf.model.raw <- regression_forest(X = XX,
                                      Y = outcome.vec,
                                      num.trees = 500,
                                      clusters = clusters.vec)
    
    varimp <- variable_importance(rf.model.raw)
    selected.idx <- which(varimp>mean(varimp))
    
    if(length(selected.idx)<2){
      selected.idx <- which(varimp>median(varimp))
    }
    
    rf.model <- regression_forest(X = XX[,selected.idx],
                                  Y = outcome.vec,
                                  num.trees=1000,
                                  clusters = clusters.vec)
    
    XX.select.0 <- XX[,selected.idx]
    XX.select.1 <- XX[,selected.idx]
    
    if("fc_al3" %in% colnames(XX.select.0)){
      XX.select.0[,"fc_al3"] <- 0
      XX.select.1[,"fc_al3"] <- 1
    }
    
    predictions.0 <- predict(rf.model, newdata = XX.select.0)$predictions
    predictions.1 <- predict(rf.model, newdata = XX.select.1)$predictions
    output <- select(data.study,
                     casenumber,
                     pid,
                     starts_with("e.hat."),
                     all_of(outcome.selection),
                     fc_al3)
    
    output$mu1 <- predictions.1
    output$mu0 <- predictions.0
  }
  
  return(output)
}

# AIPW: Estimation Function -------------------------------------------------------

aipw <- function(input.data, e.hat.selection, outcome.selection){
  
  input.data <- ungroup(input.data)
  case.vec <- input.data$casenumber
  mu.hat.1 <- input.data$mu1
  mu.hat.0 <- input.data$mu0
  e.hat <- pull(input.data, all_of(e.hat.selection))
  W <- input.data$fc_al3
  Y <- pull(input.data, all_of(outcome.selection))
  
  aipw.scores <- (mu.hat.1 - mu.hat.0
                  + W / e.hat * (Y -  mu.hat.1)
                  - (1 - W) / (1 - e.hat) * (Y -  mu.hat.0))
  ate.aipw.n <- length(aipw.scores)
  ate.aipw.j <- length(unique(case.vec))
  ate.aipw.est <- mean(aipw.scores, na.rm=TRUE)
  
  # Cluster robust standard errors
  input.data <-  input.data %>%
    mutate(scores = aipw.scores) %>%
    group_by(casenumber) %>%
    summarize(tauj=mean(scores)) %>%
    mutate(tau=mean(tauj)) %>%
    mutate(diff=(tau-tauj)^2) %>%
    ungroup()
                
  n.cases <- length(unique(case.vec))
  
  unique.cases <- as.data.frame(unique(case.vec))
  colnames(unique.cases) <- c("casenumber")
  
  unique.cases <- semi_join(unique.cases,
                            select(input.data,
                                   casenumber,
                                   diff),
                            by = "casenumber")
  
  ate.aipw.se <- sqrt(1/(n.cases*(n.cases-1)) * sum(unique.cases$diff, na.rm = TRUE))
  ate.aipw.var <- ate.aipw.se^2
  ate.aipw.upper <- ate.aipw.est + 1.96 * ate.aipw.se
  ate.aipw.lower <- ate.aipw.est - 1.96 * ate.aipw.se
  ate.aipw.tstat <- ate.aipw.est / ate.aipw.se
  ate.aipw.pvalue <- 2*(1 - pnorm(abs(ate.aipw.tstat)))
  ate.aipw.results <- c(n=ate.aipw.n,
                        j=ate.aipw.j,
                        estimate=ate.aipw.est, 
                        var.hat=ate.aipw.var,
                        std.error=ate.aipw.se,
                        ci.upper=ate.aipw.upper,
                        ci.lower=ate.aipw.lower,
                        t.stat=ate.aipw.tstat, 
                        pvalue=ate.aipw.pvalue)
  
  return(ate.aipw.results)
}

# AIPW: Pooled ATE Function & Results -------------------------------------------------------

pooled_aipw_att <- function(data.prefix, event.year, outcome.selection){
  
  data.2009 <- get(paste0(data.prefix, ".2009"))
  data.2010 <- get(paste0(data.prefix, ".2010"))
  data.2011 <- get(paste0(data.prefix, ".2011"))
  data.2012 <- get(paste0(data.prefix, ".2012"))
  
  # Lasso
  lasso.2009 <- aipw(mu_est(data.2009, outcome.selection, event.year, "lasso"), 
                     "e.hat.postlasso", outcome.selection)
  lasso.2010 <- aipw(mu_est(data.2010, outcome.selection, event.year, "lasso"), 
                     "e.hat.postlasso", outcome.selection)
  lasso.2011 <- aipw(mu_est(data.2011, outcome.selection, event.year, "lasso"), 
                     "e.hat.postlasso", outcome.selection)
  lasso.2012 <- aipw(mu_est(data.2012, outcome.selection, event.year, "lasso"), 
                     "e.hat.postlasso", outcome.selection)
  
  lasso.n <- lasso.2009[[1]] + lasso.2010[[1]] + lasso.2011[[1]] + lasso.2012[[1]]
  lasso.j <- lasso.2009[[2]] + lasso.2010[[2]] + lasso.2011[[2]] + lasso.2012[[2]]
  lasso.att <- ((lasso.2009[[1]]/lasso.n)*lasso.2009[[3]]) + 
               ((lasso.2010[[1]]/lasso.n)*lasso.2010[[3]]) +
               ((lasso.2011[[1]]/lasso.n)*lasso.2011[[3]]) +
               ((lasso.2012[[1]]/lasso.n)*lasso.2012[[3]])
  lasso.var <- ((lasso.2009[[2]]/lasso.j)^2*lasso.2009[[4]]) + 
               ((lasso.2010[[2]]/lasso.j)^2*lasso.2010[[4]]) +
               ((lasso.2011[[2]]/lasso.j)^2*lasso.2011[[4]]) +
               ((lasso.2012[[2]]/lasso.j)^2*lasso.2012[[4]])
  lasso.se <- sqrt(lasso.var)
  lasso.t <- lasso.att/lasso.se
  lasso.p <- 1.96*(1 - pnorm(abs(lasso.t)))
  
  # Forest
  rf.2009 <- aipw(mu_est(data.2009, outcome.selection, event.year, "forest"), 
                  "e.hat.rf", outcome.selection)
  rf.2010 <- aipw(mu_est(data.2010, outcome.selection, event.year, "forest"), 
                  "e.hat.rf", outcome.selection)
  rf.2011 <- aipw(mu_est(data.2011, outcome.selection, event.year, "forest"), 
                  "e.hat.rf", outcome.selection)
  rf.2012 <- aipw(mu_est(data.2012, outcome.selection, event.year, "forest"), 
                  "e.hat.rf", outcome.selection)
  
  rf.n <- rf.2009[[1]] + rf.2010[[1]] + rf.2011[[1]] + rf.2012[[1]]
  rf.j <- rf.2009[[2]] + rf.2010[[2]] + rf.2011[[2]] + rf.2012[[2]]
  rf.att <- ((rf.2009[[1]]/rf.n)*rf.2009[[3]]) + 
            ((rf.2010[[1]]/rf.n)*rf.2010[[3]]) +
            ((rf.2011[[1]]/rf.n)*rf.2011[[3]]) +
            ((rf.2012[[1]]/rf.n)*rf.2012[[3]])
  rf.var <- ((rf.2009[[2]]/rf.j)^2*rf.2009[[4]]) + 
            ((rf.2010[[2]]/rf.j)^2*rf.2010[[4]]) +
            ((rf.2011[[2]]/rf.j)^2*rf.2011[[4]]) +
            ((rf.2012[[2]]/rf.j)^2*rf.2012[[4]])
  rf.se <- sqrt(rf.var)
  rf.t <- rf.att/rf.se
  rf.p <- 1.96*(1 - pnorm(abs(rf.t)))
  
  # Tree
  tree.2009 <- aipw(mu_est(data.2009, outcome.selection, event.year, "tree"), 
                    "e.hat.tree", outcome.selection)
  tree.2010 <- aipw(mu_est(data.2010, outcome.selection, event.year, "tree"), 
                    "e.hat.tree", outcome.selection)
  tree.2011 <- aipw(mu_est(data.2011, outcome.selection, event.year, "tree"), 
                    "e.hat.tree", outcome.selection)
  tree.2012 <- aipw(mu_est(data.2012, outcome.selection, event.year, "tree"), 
                    "e.hat.tree", outcome.selection)
  
  tree.n <- tree.2009[[1]] + tree.2010[[1]] + tree.2011[[1]] + tree.2012[[1]]
  tree.j <- tree.2009[[2]] + tree.2010[[2]] + tree.2011[[2]] + tree.2012[[2]]
  tree.att <- ((tree.2009[[1]]/tree.n)*tree.2009[[3]]) + 
              ((tree.2010[[1]]/tree.n)*tree.2010[[3]]) +
              ((tree.2011[[1]]/tree.n)*tree.2011[[3]]) +
              ((tree.2012[[1]]/tree.n)*tree.2012[[3]])
  tree.var <- ((tree.2009[[2]]/tree.j)^2*tree.2009[[4]]) + 
              ((tree.2010[[2]]/tree.j)^2*tree.2010[[4]]) +
              ((tree.2011[[2]]/tree.j)^2*tree.2011[[4]]) +
              ((tree.2012[[2]]/tree.j)^2*tree.2012[[4]])
  tree.se <- sqrt(tree.var)
  tree.t <- tree.att/tree.se
  tree.p <- 1.96*(1 - pnorm(abs(tree.t)))
  
  # Output
  rf.results <- c(n=rf.n,
                  j=rf.j,
                  att=rf.att,
                  var=rf.var,
                  se=rf.se,
                  t=rf.t,
                  p=rf.p)
  
  lasso.results <- c(n=lasso.n,
                     j=lasso.j,
                     att=lasso.att,
                     var=lasso.var,
                     se=lasso.se,
                     t=lasso.t,
                     p=lasso.p)
  
  tree.results <- c(n=tree.n,
                    j=tree.j,
                    att=tree.att,
                    var=tree.var,
                    se=tree.se,
                    t=tree.t,
                    p=tree.p)

  returnme <- list(forest=rf.results, 
                   lasso=lasso.results, 
                   tree=tree.results)
  return(returnme)
}

# ATT Estimates
att.aipw.moved <- pooled_aipw_att("moved", 4, "moved")
att.aipw.school <- pooled_aipw_att("school", 4, "schoolindex_high")
att.aipw.divorce <- pooled_aipw_att("divorce", 4, "cumulnum_divorce")
att.aipw.unpaid <- pooled_aipw_att("unpaid", 4, "num00_collection_unpaid")

# Observe all side-by-side
att.aipw.moved
att.strat.moved
att.aipw.school
att.strat.school
att.aipw.divorce
att.strat.divorce
att.aipw.unpaid
att.strat.unpaid
# AIPW: Event Studies -----------------------------------------------------

## Advanced warning: this is a crude implementation
evt_study_plot_data <- function(prefix.selection, outcome.selection, method.selection){
  
  if(method.selection == "forest"){
    e.hat.select <- "e.hat.rf"
  }
  
  if(method.selection == "tree"){
    e.hat.select <- "e.hat.tree" 
  }
  
  if(method.selection == "lasso"){
    e.hat.select <- "e.hat.postlasso"
  }
  
  data.2009 <- get(paste0(prefix.selection,".2009"))
  data.2010 <- get(paste0(prefix.selection,".2010"))
  data.2011 <- get(paste0(prefix.selection,".2011"))
  data.2012 <- get(paste0(prefix.selection,".2012"))
  
  # I hate doing it like this but I am running out of time
  y09.m4 <- aipw(mu_est(data.2009, outcome.selection, -4, method.selection),e.hat.select, outcome.selection)
  y09.m3 <- aipw(mu_est(data.2009, outcome.selection, -3, method.selection),e.hat.select, outcome.selection)
  y09.m2 <- aipw(mu_est(data.2009, outcome.selection, -2, method.selection),e.hat.select, outcome.selection)
  y09.m1 <- aipw(mu_est(data.2009, outcome.selection, -1, method.selection),e.hat.select, outcome.selection)
  y09.0 <- aipw(mu_est(data.2009, outcome.selection, 0, method.selection),e.hat.select, outcome.selection)
  y09.p1 <- aipw(mu_est(data.2009, outcome.selection, 1, method.selection),e.hat.select, outcome.selection)
  y09.p2 <- aipw(mu_est(data.2009, outcome.selection, 2, method.selection),e.hat.select, outcome.selection)
  y09.p3 <- aipw(mu_est(data.2009, outcome.selection, 3, method.selection),e.hat.select, outcome.selection)
  y09.p4 <- aipw(mu_est(data.2009, outcome.selection, 4, method.selection),e.hat.select, outcome.selection)
  
  y10.m4 <- aipw(mu_est(data.2010, outcome.selection, -4, method.selection),e.hat.select, outcome.selection)
  y10.m3 <- aipw(mu_est(data.2010, outcome.selection, -3, method.selection),e.hat.select, outcome.selection)
  y10.m2 <- aipw(mu_est(data.2010, outcome.selection, -2, method.selection),e.hat.select, outcome.selection)
  y10.m1 <- aipw(mu_est(data.2010, outcome.selection, -1, method.selection),e.hat.select, outcome.selection)
  y10.0 <- aipw(mu_est(data.2010, outcome.selection, 0, method.selection),e.hat.select, outcome.selection)
  y10.p1 <- aipw(mu_est(data.2010, outcome.selection, 1, method.selection),e.hat.select, outcome.selection)
  y10.p2 <- aipw(mu_est(data.2010, outcome.selection, 2, method.selection),e.hat.select, outcome.selection)
  y10.p3 <- aipw(mu_est(data.2010, outcome.selection, 3, method.selection),e.hat.select, outcome.selection)
  y10.p4 <- aipw(mu_est(data.2010, outcome.selection, 4, method.selection),e.hat.select, outcome.selection)
  
  y11.m4 <- aipw(mu_est(data.2011, outcome.selection, -4, method.selection),e.hat.select, outcome.selection)
  y11.m3 <- aipw(mu_est(data.2011, outcome.selection, -3, method.selection),e.hat.select, outcome.selection)
  y11.m2 <- aipw(mu_est(data.2011, outcome.selection, -2, method.selection),e.hat.select, outcome.selection)
  y11.m1 <- aipw(mu_est(data.2011, outcome.selection, -1, method.selection),e.hat.select, outcome.selection)
  y11.0 <- aipw(mu_est(data.2011, outcome.selection, 0, method.selection),e.hat.select, outcome.selection)
  y11.p1 <- aipw(mu_est(data.2011, outcome.selection, 1, method.selection),e.hat.select, outcome.selection)
  y11.p2 <- aipw(mu_est(data.2011, outcome.selection, 2, method.selection),e.hat.select, outcome.selection)
  y11.p3 <- aipw(mu_est(data.2011, outcome.selection, 3, method.selection),e.hat.select, outcome.selection)
  y11.p4 <- aipw(mu_est(data.2011, outcome.selection, 4, method.selection),e.hat.select, outcome.selection)
  
  y12.m4 <- aipw(mu_est(data.2012, outcome.selection, -4, method.selection),e.hat.select, outcome.selection)
  y12.m3 <- aipw(mu_est(data.2012, outcome.selection, -3, method.selection),e.hat.select, outcome.selection)
  y12.m2 <- aipw(mu_est(data.2012, outcome.selection, -2, method.selection),e.hat.select, outcome.selection)
  y12.m1 <- aipw(mu_est(data.2012, outcome.selection, -1, method.selection),e.hat.select, outcome.selection)
  y12.0 <- aipw(mu_est(data.2012, outcome.selection, 0, method.selection),e.hat.select, outcome.selection)
  y12.p1 <- aipw(mu_est(data.2012, outcome.selection, 1, method.selection),e.hat.select, outcome.selection)
  y12.p2 <- aipw(mu_est(data.2012, outcome.selection, 2, method.selection),e.hat.select, outcome.selection)
  y12.p3 <- aipw(mu_est(data.2012, outcome.selection, 3, method.selection),e.hat.select, outcome.selection)
  y12.p4 <- aipw(mu_est(data.2012, outcome.selection, 4, method.selection),e.hat.select, outcome.selection)
  
  fixes <- c("m4", "m3", "m2", "m1", "0", "p1", "p2", "p3", "p4")
  att.2009 <- rep(NA, 9)
  ciu.2009 <- rep(NA, 9)
  cil.2009 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2009[i] <- get(paste0("y09.", checkme))[3]
    ciu.2009[i] <- get(paste0("y09.", checkme))[6]
    cil.2009[i] <- get(paste0("y09.", checkme))[7]
  }
  
  att.2010 <- rep(NA, 9)
  ciu.2010 <- rep(NA, 9)
  cil.2010 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2010[i] <- get(paste0("y10.", checkme))[3]
    ciu.2010[i] <- get(paste0("y10.", checkme))[6]
    cil.2010[i] <- get(paste0("y10.", checkme))[7]
  }
  
  att.2011 <- rep(NA, 9)
  ciu.2011 <- rep(NA, 9)
  cil.2011 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2011[i] <- get(paste0("y11.", checkme))[3]
    ciu.2011[i] <- get(paste0("y11.", checkme))[6]
    cil.2011[i] <- get(paste0("y11.", checkme))[7]
  }
  
  att.2012 <- rep(NA, 9)
  ciu.2012 <- rep(NA, 9)
  cil.2012 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2012[i] <- get(paste0("y12.", checkme))[3]
    ciu.2012[i] <- get(paste0("y12.", checkme))[6]
    cil.2012[i] <- get(paste0("y12.", checkme))[7]
  }
  
  evt.plot.2009 <- as.data.frame(cbind(att.2009, ciu.2009, cil.2009))
  evt.plot.2009$Cohort = rep("2009", 9)
  colnames(evt.plot.2009) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2009$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot.2010<- as.data.frame(cbind(att.2010, ciu.2010, cil.2010))
  evt.plot.2010$Cohort = rep("2010", 9)
  colnames(evt.plot.2010) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2010$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot.2011 <- as.data.frame(cbind(att.2011, ciu.2011, cil.2011))
  evt.plot.2011$Cohort = rep("2011", 9)
  colnames(evt.plot.2011) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2011$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot.2012 <- as.data.frame(cbind(att.2012, ciu.2012, cil.2012))
  evt.plot.2012$Cohort = rep("2012", 9)
  colnames(evt.plot.2012) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2012$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot <- rbind(evt.plot.2009, evt.plot.2010, evt.plot.2011, evt.plot.2012)
  
  return(evt.plot)
  
}

evt_study_plot_gg <- function(plot.data, title.info, subtitle){
  
  ggplot(data=as.data.frame(plot.data), aes(x=Evt, y=ATT, colour=Cohort)) + 
    geom_point() + 
    geom_errorbar(aes(x=Evt, ymin=CIL, ymax=CIU, colour=Cohort)) + 
    facet_wrap(~Cohort) + 
    ylab("Cross-Sectional ATE Estimate") + 
    geom_vline(aes(group=Cohort, xintercept=-0.5), linetype="dashed", alpha=1, color="black") + 
    geom_hline(aes(yintercept=0), linetype="dashed", alpha=0.25, color="black") + 
    theme_bw() + 
    ggtitle(paste("Cohort Event Study:", title.info), 
            subtitle=paste(subtitle)) + 
    scale_x_continuous(name="Event Year", 
                       limits=c(-4.5, 4.5),
                       breaks=c(-4, -3, -2, -1, 0,
                                1, 2, 3, 4))
  
}

plot.event.forest.school <- evt_study_plot_data("school", "schoolindex_high", "forest")
plot.event.forest.moved <- evt_study_plot_data("moved", "moved", "forest")
plot.event.forest.divorce <- evt_study_plot_data("divorce", "cumulnum_divorce", "forest")
plot.event.forest.unpaid <- evt_study_plot_data("unpaid", "num00_collection_unpaid", "forest")

# Plots
pdf("./output/evt_study_aipw_forest_school.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.forest.school,
                  "High School Test Index",
                  "AIPW w/ Random Forest Nuisance Parameter Estimation")
dev.off()

pdf("./output/evt_study_aipw_forest_moved.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.forest.moved,
                  "Moved",
                  "AIPW w/ Random Forest Nuisance Parameter Estimation")
dev.off()

pdf("./output/evt_study_aipw_forest_divorce.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.forest.divorce,
                  "Cumulative Number of Divorces",
                  "AIPW w/ Random Forest Nuisance Parameter Estimation")
dev.off()

pdf("./output/evt_study_aipw_forest_unpaid.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.forest.unpaid,
                  "Number of Unpaid Collections",
                  "AIPW w/ Random Forest Nuisance Parameter Estimation")
dev.off()




# Stratification: Event Studies -------------------------------------------
strat_evt_study_plot_data <- function(prefix.selection, outcome.selection, method.selection){
  
  if(method.selection == "forest"){
    e.hat.select <- "e.hat.rf"
  }
  
  if(method.selection == "tree"){
    e.hat.select <- "e.hat.tree" 
  }
  
  if(method.selection == "lasso"){
    e.hat.select <- "e.hat.postlasso"
  }
  
  data.2009 <- get(paste0(prefix.selection,".2009"))
  data.2010 <- get(paste0(prefix.selection,".2010"))
  data.2011 <- get(paste0(prefix.selection,".2011"))
  data.2012 <- get(paste0(prefix.selection,".2012"))
  
  # I hate doing it like this but I am running out of time
  y09.m4 <- strat_est(data.2009, -4, outcome.selection ,e.hat.select)
  y09.m3 <- strat_est(data.2009, -3, outcome.selection ,e.hat.select)
  y09.m2 <- strat_est(data.2009, -2, outcome.selection ,e.hat.select)
  y09.m1 <- strat_est(data.2009, -1, outcome.selection ,e.hat.select)
  y09.0 <- strat_est(data.2009, 0, outcome.selection ,e.hat.select)
  y09.p1 <- strat_est(data.2009, 1, outcome.selection ,e.hat.select)
  y09.p2 <- strat_est(data.2009, 2, outcome.selection ,e.hat.select)
  y09.p3 <- strat_est(data.2009, 3, outcome.selection ,e.hat.select)
  y09.p4 <- strat_est(data.2009, 4, outcome.selection ,e.hat.select)
  
  y10.m4 <- strat_est(data.2010, -4, outcome.selection ,e.hat.select)
  y10.m3 <- strat_est(data.2010, -3, outcome.selection ,e.hat.select)
  y10.m2 <- strat_est(data.2010, -2, outcome.selection ,e.hat.select)
  y10.m1 <- strat_est(data.2010, -1, outcome.selection ,e.hat.select)
  y10.0 <- strat_est(data.2010, 0, outcome.selection ,e.hat.select)
  y10.p1 <- strat_est(data.2010, 1, outcome.selection ,e.hat.select)
  y10.p2 <- strat_est(data.2010, 2, outcome.selection ,e.hat.select)
  y10.p3 <- strat_est(data.2010, 3, outcome.selection ,e.hat.select)
  y10.p4 <- strat_est(data.2010, 4, outcome.selection ,e.hat.select)
  
  y11.m4 <- strat_est(data.2011, -4, outcome.selection ,e.hat.select)
  y11.m3 <- strat_est(data.2011, -3, outcome.selection ,e.hat.select)
  y11.m2 <- strat_est(data.2011, -2, outcome.selection ,e.hat.select)
  y11.m1 <- strat_est(data.2011, -1, outcome.selection ,e.hat.select)
  y11.0 <- strat_est(data.2011, 0, outcome.selection ,e.hat.select)
  y11.p1 <- strat_est(data.2011, 1, outcome.selection ,e.hat.select)
  y11.p2 <- strat_est(data.2011, 2, outcome.selection ,e.hat.select)
  y11.p3 <- strat_est(data.2011, 3, outcome.selection ,e.hat.select)
  y11.p4 <- strat_est(data.2011, 4, outcome.selection ,e.hat.select)
  
  y12.m4 <- strat_est(data.2012, -4, outcome.selection ,e.hat.select)
  y12.m3 <- strat_est(data.2012, -3, outcome.selection ,e.hat.select)
  y12.m2 <- strat_est(data.2012, -2, outcome.selection ,e.hat.select)
  y12.m1 <- strat_est(data.2012, -1, outcome.selection ,e.hat.select)
  y12.0 <- strat_est(data.2012, 0, outcome.selection ,e.hat.select)
  y12.p1 <- strat_est(data.2012, 1, outcome.selection ,e.hat.select)
  y12.p2 <- strat_est(data.2012, 2, outcome.selection ,e.hat.select)
  y12.p3 <- strat_est(data.2012, 3, outcome.selection ,e.hat.select)
  y12.p4 <- strat_est(data.2012, 4, outcome.selection ,e.hat.select)
  
  fixes <- c("m4", "m3", "m2", "m1", "0", "p1", "p2", "p3", "p4")
  att.2009 <- rep(NA, 9)
  ciu.2009 <- rep(NA, 9)
  cil.2009 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2009[i] <- get(paste0("y09.", checkme))[3]
    ciu.2009[i] <- get(paste0("y09.", checkme))[5]
    cil.2009[i] <- get(paste0("y09.", checkme))[6]
  }
  
  att.2010 <- rep(NA, 9)
  ciu.2010 <- rep(NA, 9)
  cil.2010 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2010[i] <- get(paste0("y10.", checkme))[3]
    ciu.2010[i] <- get(paste0("y10.", checkme))[5]
    cil.2010[i] <- get(paste0("y10.", checkme))[6]
  }
  
  att.2011 <- rep(NA, 9)
  ciu.2011 <- rep(NA, 9)
  cil.2011 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2011[i] <- get(paste0("y11.", checkme))[3]
    ciu.2011[i] <- get(paste0("y11.", checkme))[5]
    cil.2011[i] <- get(paste0("y11.", checkme))[6]
  }
  
  att.2012 <- rep(NA, 9)
  ciu.2012 <- rep(NA, 9)
  cil.2012 <- rep(NA, 9)
  for(i in 1:9){
    checkme <- fixes[i]
    att.2012[i] <- get(paste0("y12.", checkme))[3]
    ciu.2012[i] <- get(paste0("y12.", checkme))[5]
    cil.2012[i] <- get(paste0("y12.", checkme))[6]
  }
  
  evt.plot.2009 <- as.data.frame(cbind(att.2009, ciu.2009, cil.2009))
  evt.plot.2009$Cohort = rep("2009", 9)
  colnames(evt.plot.2009) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2009$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot.2010<- as.data.frame(cbind(att.2010, ciu.2010, cil.2010))
  evt.plot.2010$Cohort = rep("2010", 9)
  colnames(evt.plot.2010) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2010$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot.2011 <- as.data.frame(cbind(att.2011, ciu.2011, cil.2011))
  evt.plot.2011$Cohort = rep("2011", 9)
  colnames(evt.plot.2011) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2011$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot.2012 <- as.data.frame(cbind(att.2012, ciu.2012, cil.2012))
  evt.plot.2012$Cohort = rep("2012", 9)
  colnames(evt.plot.2012) <- c("ATT", "CIU", "CIL", "Cohort")
  evt.plot.2012$Evt <- c(-4, -3, -2, -1, 0, 1, 2, 3, 4)
  
  evt.plot <- rbind.data.frame(evt.plot.2009, evt.plot.2010, evt.plot.2011, evt.plot.2012)
  
  return(evt.plot)
  
}

plot.event.strat.school <- strat_evt_study_plot_data("school", "schoolindex_high", "forest")
plot.event.strat.school <- tidyr::unnest(plot.event.strat.school,cols = c(ATT, CIU, CIL) )
plot.event.strat.moved <- strat_evt_study_plot_data("moved", "moved", "forest")
plot.event.strat.moved <- tidyr::unnest(plot.event.strat.moved,cols = c(ATT, CIU, CIL)) 
plot.event.strat.divorce <- strat_evt_study_plot_data("divorce", "cumulnum_divorce", "forest")
plot.event.strat.divorce <- tidyr::unnest(plot.event.strat.divorce,cols = c(ATT, CIU, CIL)) 
plot.event.strat.unpaid <- strat_evt_study_plot_data("unpaid", "num00_collection_unpaid", "forest")
plot.event.strat.unpaid <- tidyr::unnest(plot.event.strat.unpaid,cols = c(ATT, CIU, CIL)) 

pdf("./output/evt_study_irs_forest_school.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.strat.school,
                  "High School Test Index",
                  "Imbens-Rubin Subclassification Estimator w/ Random Forest Propensity Score Estimation")
dev.off()

pdf("./output/evt_study_irs_forest_moved.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.strat.moved,
                  "Moved",
                  "Imbens-Rubin Subclassification Estimator w/ Random Forest Propensity Score Estimation")
dev.off()

pdf("./output/evt_study_irs_forest_divorce.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.strat.divorce,
                  "Cumulative Number of Divorces",
                  "Imbens-Rubin Subclassification Estimator w/ Random Forest Propensity Score Estimation")
dev.off()

pdf("./output/evt_study_irs_forest_unpaid.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.strat.unpaid,
                  "Number of Unpaid Collections",
                  "Imbens-Rubin Subclassification Estimator w/ Random Forest Propensity Score Estimation")
dev.off()

# Of interest due to different point estimates
plot.event.strat.school.lasso <- strat_evt_study_plot_data("school", "schoolindex_high", "lasso")
plot.event.strat.school.lasso <- tidyr::unnest(plot.event.strat.school.lasso,cols = c(ATT, CIU, CIL) )
pdf("./output/evt_study_irs_lasso_school.pdf", width=8.67, height=5.75)  
evt_study_plot_gg(plot.event.strat.school.lasso,
                  "High School Test Index",
                  "Imbens-Rubin Subclassification Estimator w/ Post-LASSO Propensity Score Estimation")
dev.off()