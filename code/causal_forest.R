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
library(coeftest)
library(sandwich)
set.seed(080699)

#######
# In places there are commented out code followed by a data read, 
# this is because some code runs on the Sherlock HPC cluster, 
# and I pull the data back onto my local machine for graphics processing.
#######

# Panel-Ify ---------------------------------------------------------------

panel_me <- function(data.read, outcome){
  
  baseline <- subset(data.read, pull(data.read, evt)==-1)
  baseline <- mean(pull(baseline, outcome), na.rm=TRUE)
  
  data.read[, all_of(outcome)] <- pull(data.read, outcome) - baseline
  
  
  return(data.read)
  
}

# Import Balanced Panel Data ----------------------------------------------
moved.balanced.cohorts <- fread("./data/moved_balanced_cohorts.csv")
moved.2009 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2009), "moved")
moved.2010 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2010), "moved")
moved.2011 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2011), "moved")
moved.2012 <- panel_me(subset(moved.balanced.cohorts, year_beg_fs == 2012), "moved")

school.balanced.cohorts <- fread("./data/school_balanced_cohorts.csv")
school.2009 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2009), "schoolindex_high")
school.2010 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2010), "schoolindex_high")
school.2011 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2011), "schoolindex_high")
school.2012 <- panel_me(subset(school.balanced.cohorts, year_beg_fs == 2012), "schoolindex_high")

divorce.balanced.cohorts <- fread("./data/divorce.balanced.cohorts")
divorce.2009 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2009), "cumulnum_divorce")
divorce.2010 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2010), "cumulnum_divorce")
divorce.2011 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2011), "cumulnum_divorce")
divorce.2012 <- panel_me(subset(divorce.balanced.cohorts, year_beg_fs == 2012), "cumulnum_divorce")

unpaid.balanced.cohorts <- fread("./data/unpaid.balanced.cohorts")
unpaid.2009 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2009), "num00_collection_unpaid")
unpaid.2010 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2010), "num00_collection_unpaid")
unpaid.2011 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2011), "num00_collection_unpaid")
unpaid.2012 <- panel_me(subset(unpaid.balanced.cohorts, year_beg_fs == 2012), "num00_collection_unpaid")

# Also main data for merges
data <- fread("./data/analysis_data_owners.csv")

# Function: Causal Forest Cohort-Year -------------------------------------
## Important outputs here: the ATT & associated estimates, as well as the 
## vector of CATE estimates. This will let us build event study plots as 
## well as do the relevant CATE analysis when the time comes in a moment. 

covariates.prediction <- c("moved", "isowner", "log_finishedsqft_dq", "zipinc_log",
                           "schoolindex_elem", "schoolindex_middle", "schoolindex_high",
                           "cumulnum_divorce", "cumulnum_crime", "cumulnum_dui", 
                           "cumulnum_bankruptcy", "vantage", "num00_foreclosure", 
                           "num00_collection_unpaid", "num00_trade_auto", "num12_trade30", 
                           "num00_mortgage_loanmod", "num00_mortgage_open", 
                           "rtot00_mortgage_open_bal", "rtot00_mortgage_monthpay")

cf_estimate <- function(data.input, outcome.selection, event.year){
  
  ###
  #Going to use the ex-ante estimated e.hats from pre-treatment. 
  ###
  
  covariates <- covariates.prediction[covariates.prediction!=all_of(outcome.selection)]
  data.study <- left_join(subset(data.input, evt==event.year), 
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
  
  # Training Formula 
  train.fmla <- formula(paste(outcome.selection, "~",
                              paste(covariates, collapse="+")))
  
  XX <- model.matrix.lm(train.fmla, 
                        data.study,
                        na.action = NULL)
  
  outcome.vec <-  pull(data.study, 
                        all_of(outcome.selection))
  
  treat.vec <- pull(data.study, fc_al3)
  
  clusters.vec <- as.factor(pull(data.study, 
                                 casenumber))
  
  e.hat.vec <- pull(data.study,
                    e.hat.rf)
  
  cf.raw <- causal_forest(XX,
                          outcome.vec,
                          treat.vec,
                          clusters = clusters.vec,
                          W.hat = e.hat.vec,
                          num.trees = 500
                          )
  
  varimp <- variable_importance(cf.raw)
  selected.idx <- which(varimp>mean(varimp))
  
  if(length(selected.idx)<2){
    selected.idx <- which(varimp>median(varimp))
  }
  
  cf.model <- causal_forest(XX[,selected.idx],
                            outcome.vec,
                            treat.vec, 
                            num.trees=1000,
                            W.hat = e.hat.vec,
                            clusters = clusters.vec)
  
  tau.hat <- predict(cf.model)$predictions
  ate.hold <- average_treatment_effect(cf.model)
  
  # N should be labelled J but I made the mistake 2000 lines of 
  # code ago and am following through, sorry
  ate.aipw.n <- nrow(XX)
  ate.aipw.j <- length(unique(clusters.vec))
  ate.aipw.est <- ate.hold[1]
  ate.aipw.se <- ate.hold[2]
  ate.aipw.var <- ate.aipw.se^2
  ate.aipw.upper <- ate.aipw.est + 1.96 * ate.aipw.se
  ate.aipw.lower <- ate.aipw.est - 1.96 * ate.aipw.se
  ate.aipw.tstat <- ate.aipw.est / ate.aipw.se
  ate.aipw.pvalue <- 2*(pnorm(1 - abs(ate.aipw.tstat)))
  ate.aipw.results <- c(n=ate.aipw.n,
                        j=ate.aipw.j,
                        estimate=ate.aipw.est, 
                        var.hat=ate.aipw.var,
                        std.error=ate.aipw.se,
                        ci.upper=ate.aipw.upper,
                        ci.lower=ate.aipw.lower,
                        t.stat=ate.aipw.tstat, 
                        pvalue=ate.aipw.pvalue)
  
  data.output <- subset(data.input, evt==event.year)
  Y.hat <- cf.model$Y.hat
  W.hat <- cf.model$W.hat
  data.output[,"tau.hat.cf"] <- tau.hat
  data.output[,"y.hat.cf"] <- Y.hat
  data.output[,"w.hat.cf"] <- W.hat
  
  output.me <- list(ate.aipw.results,
                    data.output,
                    cf.model)
  
  return(output.me)
}

# Function: Pooled Event Studies ------------------------------------------
cf_pooled <- function(data.prefix, outcome.selection, event.year){
  
  data.2009 <- get(paste0(data.prefix, ".2009"))
  data.2010 <- get(paste0(data.prefix, ".2010"))
  data.2011 <- get(paste0(data.prefix, ".2011"))
  data.2012 <- get(paste0(data.prefix, ".2012"))
  
  cf.2009 <- cf_estimate(data.2009, outcome.selection, event.year)
  cf.2010 <- cf_estimate(data.2009, outcome.selection, event.year)
  cf.2011 <- cf_estimate(data.2009, outcome.selection, event.year)
  cf.2012 <- cf_estimate(data.2009, outcome.selection, event.year)
  
  ## CATE estimates
  cf.2009.cate <- select(cf.2009[[2]],
                         casenumber,
                         tau.hat.cf,
                         fc_al3,
                         all_of(outcome.selection),
                         y.hat.cf,
                         w.hat.cf)
  
  cf.2010.cate <- select(cf.2010[[2]],
                         casenumber,
                         tau.hat.cf,
                         fc_al3,
                         all_of(outcome.selection),
                         y.hat.cf,
                         w.hat.cf)
  
  cf.2011.cate <- select(cf.2011[[2]],
                         casenumber,
                         tau.hat.cf,
                         fc_al3,
                         all_of(outcome.selection),
                         y.hat.cf,
                         w.hat.cf)
  
  cf.2012.cate <- select(cf.2012[[2]],
                         casenumber,
                         tau.hat.cf,
                         fc_al3,
                         all_of(outcome.selection),
                         y.hat.cf,
                         w.hat.cf)
  
  cf.cate <- rbind.data.frame(cf.2009.cate,
                              cf.2010.cate,
                              cf.2011.cate,
                              cf.2012.cate)
  
  ## ATE
  cf.n <- cf.2009[[1]][1] + 
          cf.2010[[1]][1] + 
          cf.2011[[1]][1] + 
          cf.2012[[1]][1]
  
  cf.j <- cf.2009[[1]][2] + 
          cf.2010[[1]][2] + 
          cf.2011[[1]][2] + 
          cf.2012[[1]][2]
  
  cf.att <- ((cf.2009[[1]][1]/cf.n)*cf.2009[[1]][3]) + 
            ((cf.2010[[1]][1]/cf.n)*cf.2010[[1]][3]) + 
            ((cf.2011[[1]][1]/cf.n)*cf.2011[[1]][3]) + 
            ((cf.2012[[1]][1]/cf.n)*cf.2012[[1]][3])

  cf.var <- ((cf.2009[[1]][2]/cf.j)^2*cf.2009[[1]][4]) + 
            ((cf.2010[[1]][2]/cf.j)^2*cf.2010[[1]][4]) + 
            ((cf.2011[[1]][2]/cf.j)^2*cf.2011[[1]][4]) + 
            ((cf.2012[[1]][2]/cf.j)^2*cf.2012[[1]][4])

  cf.se <- sqrt(cf.var)
  cf.t <- cf.att/cf.se
  cf.p <- 1.96*(1-pnorm(abs(cf.t)))
  
  cf.results <- c(n=cf.n, 
                  j=cf.j,
                  att=cf.att, 
                  var=cf.var,
                  se=cf.se,
                  t=cf.t,
                  p=cf.p)

  # Return
  returnme <- list(cf.results, cf.cate)
  
  return(returnme)
}

# Causal Forest: ATTs  -----------------------------------------------------

moved.cf <- cf_pooled("moved", "moved", 4)
divorce.cf <- cf_pooled("divorce", "cumulnum_divorce", 4)
school.cf <- cf_pooled("school", "schoolindex_high", 4)
unpaid.cf <- cf_pooled("unpaid", "num00_collection_unpaid", 4)

moved.cf[[1]]
divorce.cf[[1]]
school.cf[[1]]
unpaid.cf[[1]]

# Function: Direct Calibration Test Implementation ------------------------

pooled_test_calibration <- function(pooled.output.data, outcome.selection){
  
  # Note that all predictions in here are cross fit
  
  Y <- pull(pooled.output.data,
            outcome.selection)
  
  W <- pull(pooled.output.data,
            fc_al3)
  
  m.hat <- pull(pooled.output.data,
                y.hat.cf)
  
  e.hat <- pull(pooled.output.data,
                w.hat.cf)
  
  tau.hat <- pull(pooled.output.data,
                  tau.hat.cf)
  
  w.part <- W - e.hat
  t.part <- tau.hat - mean(tau.hat)
  
  y.reg <- Y - m.hat
  w.reg <- mean(tau.hat)*(w.part)
  t.reg <- t.part*w.part
  
  reg <- lm(y.reg ~ 0 + w.reg + t.reg)
  output.model <- summary(reg)
  output.HC3 <- coeftest(reg, vcov=vcovHC(reg))
  
  output <- list(model = output.model,
                 HC3 = output.HC3)
  
  return(output)
  
}

pooled_test_calibration(moved.cf[[2]], "moved")
pooled_test_calibration(divorce.cf[[2]], "cumulnum_divorce")
pooled_test_calibration(school.cf[[2]], "schoolindex_high")
pooled_test_calibration(unpaid.cf[[2]], "num00_collection_unpaid")

pdf("./output/cate_moved_hist.pdf", width=8.67, height=5.75)  
hist(moved.cf[[2]]$tau.hat.cf)
dev.off()

pdf("./output/cate_divorce_hist.pdf", width=8.67, height=5.75)  
hist(divorce.cf[[2]]$tau.hat.cf)
dev.off()

pdf("./output/cate_school_hist.pdf", width=8.67, height=5.75)  
hist(school.cf[[2]]$tau.hat.cf)
dev.off()

pdf("./output/cate_unpaid_hist.pdf", width=8.67, height=5.75)  
hist(unpaid.cf[[2]]$tau.hat.cf)
dev.off()

# Data-Driven Subgroups ---------------------------------------------------

## Implementing this for clustered data is a SERIOUS pain, since you
## can't just use the cluster argument in the causal_forest function. 
## So we need to do our own K-folding (cluster-robust), and then 
## manually do things like calculate within-quintile means and standard
## deviations. To account for the number of times this runs, I cut the 
## number of forest runs by 10 (for 10-fold); I think this emulates how
## the cluster argument would function. There is also the matter that 
## these are pooled event studies, so we need to do this mulitple times. 
## Note that I am making a serious choice to calculate quintiles within-
## cohort and then aggregate. This is to avoid just recovering that one 
## cohort had a weak TE, though of course that is potentially interesting. 

subgroup_cohort <- function(data.input, outcome.selection, event.year){
  
  covariates <- covariates.prediction[covariates.prediction!=all_of(outcome.selection)]
  data.study <- left_join(subset(data.input, evt==event.year), 
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

  # Training Formula 
  train.fmla <- formula(paste(outcome.selection, "~",
                              paste(covariates, collapse="+")))
  
  # Arguments
  num.folds <- 5
  data.study <- ungroup(data.study)
  data.study <- mutate(data.study,
                       casenumber = as.factor(casenumber))
  print(data.study)
  print(sum(is.na(data.study$casenumber)))
  data.study <- fold(data.study,
                     k=num.folds,
                     id_col="casenumber")
  
  data.study <- ungroup(data.study)
  
  Y <- as.matrix(data.study[complete.cases(select(data.study,
                                                  covariates)),"outcome"])
  W <- as.matrix(data.study[complete.cases(select(data.study,
                                                  covariates)),"treatment"])
  folds <- as.matrix(data.study[complete.cases(select(data.study,
                                            covariates)),".folds"])
  XX <- model.matrix.lm(train.fmla, data.study)
  forest <- causal_forest(XX,
                          Y,
                          W,
                          clusters = folds)
  
  tau.hat <- predict(forest)$predictions
  
  ranking <- rep(NA, n)
  for (fold in seq(num.folds)) {
    tau.hat.quantiles <- quantile(tau.hat[folds == fold], probs = seq(0, 1, by=1/num.rankings))
    ranking[folds == fold] <- cut(tau.hat[folds == fold], tau.hat.quantiles, include.lowest=TRUE,labels=seq(num.rankings))
  }
  
  return(ranking)
  
  # Prepare for below loop
  #n <- nrow(data.study)
  #ranking <- rep(NA, n)
  #tau.hat.out <- rep(NA,n)
  #W.hat.out <- rep(NA,n)
  #e.hat.out <- rep(NA,n)
  #num.rankings <- 5  
  #
  # Looping over folds
  #for(fold in 1:num.folds){
  #  
  #  # Within-fold data
  #  X.predict <- model.matrix.lm(train.fmla, data.study[data.study$.folds==fold,], na.action=NULL)
  #  
  #  # Out-of-fold data
  #  X.train <- model.matrix.lm(train.fmla, data.study[data.study$.folds!=fold,], na.action=NULL)
  #  print(nrow(X.train))
  #  Y.train <- as.matrix(data.study[data.study$.folds!=fold,"outcome"])
  #  print(nrow(Y.train))
  #  W.train <- as.matrix(data.study[data.study$.folds!=fold,"treatment"])
  #  print(nrow(W.train))
  #  
  #  cluster.train <- as.matrix((data.study[data.study$.folds!=fold,"casenumber"]))
  #  
  #  # Forest on out-of-fold data, clustering on cluster variable
  #  forest <- causal_forest(X.train, 
  #                          Y.train, 
  #                          W.train, 
  #                          clusters = cluster.train,
  #                          num.trees = 200)
  #  
  #  # Predict on within-fold data, rank within-fold
  #  tau.hat <- predict(forest, X.predict)$predictions
  #  tau.hat.quantiles <- quantile(tau.hat, probs = seq(0, 1, by=1/num.rankings))
  #  ranking[data.study$.folds == fold] <- cut(tau.hat, 
  #                                           tau.hat.quantiles, 
  #                                           include.lowest=TRUE,
  #                                           labels=seq(num.rankings))
  #  
  #  # Gather for output
  #  tau.hat.out[data.study$.folds == fold] <- predict(forest, X.predict)$predictions
  #  W.hat.out[data.study$.folds == fold] <- predict(forest, X.predict)$Y.hat
  #  e.hat.out[data.study$.folds == fold] <- predict(forest, X.predict)$W.hat
  #}
  
  
}

#subgroup_cohort(moved.2009, "moved", 4)


#data_subgroup <- function(data.prefix, outcome.selection, event.year){
 
  data.2009 <- get(paste0(data.prefix, ".2009"))
  data.2010 <- get(paste0(data.prefix, ".2010"))
  data.2011 <- get(paste0(data.prefix, ".2011"))
  data.2012 <- get(paste0(data.prefix, ".2012"))
  
  for(year in 2009:2012){
    
    covariates <- covariates.prediction[covariates.prediction!=all_of(outcome.selection)]
    data.study <- left_join(subset(data.input, evt==event.year), 
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
    
    # Training Formula 
    train.fmla <- formula(paste(outcome.selection, "~",
                                paste(covariates, collapse="+")))
    
    XX <- model.matrix.lm(train.fmla, 
                          data.study,
                          na.action = NULL)
    
    outcome.vec <-  pull(data.study, 
                         all_of(outcome.selection))
    
    treat.vec <- pull(data.study, fc_al3)
    
    clusters.vec <- as.factor(pull(data.study, 
                                   casenumber))
    
    e.hat.vec <- pull(data.study,
                      e.hat.rf)
    
  }
  
}