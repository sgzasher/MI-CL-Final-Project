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
library(lmtest)
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
## I am aware that this tau_hat averaging is potentially not valid, so
## refrain from inference. 

subgroup_cohort <- function(data.input, outcome.selection, event.year){
  
  covariates <- covariates.prediction[covariates.prediction!=all_of(outcome.selection)]
  data.study <- left_join(subset(data.input, evt==event.year), 
                          select(data, 
                                 casenumber,
                                 pid, 
                                 year_orig,
                                 all_of(covariates)),
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

  
  data.study <- fold(data.study,
                     k=num.folds,
                     id_col="casenumber")
  
  data.study <- ungroup(data.study)
  
  data.study$outcome <- data.study[,all_of(outcome.selection)]
  data.study$treatment <- data.study[,"fc_al3"]
  data.study <- data.study %>%
                mutate(casenumber2 = as.numeric(as.factor(casenumber)))
  
  
  # Prepare for below loop
  n <- nrow(data.study)
  ranking <- rep(NA, n)
  tau.hat.out <- rep(NA,n)
  num.rankings <- 5  
  
  # Looping over folds
  for(fold in 1:num.folds){
    
    # Within-fold data
    X.predict <- model.matrix.lm(train.fmla, data.study[data.study$.folds==fold,], na.action=NULL)
    
    # Out-of-fold data
    X.train <- model.matrix.lm(train.fmla, data.study[data.study$.folds!=fold,], na.action=NULL)
    Y.train <- as.matrix(data.study[data.study$.folds!=fold,"outcome"])
    W.train <- as.matrix(data.study[data.study$.folds!=fold,"treatment"])

    cluster.train <- as.matrix(data.study[data.study$.folds!=fold,"casenumber2"])
    
    
    # Forest on out-of-fold data, clustering on cluster variable
    forest <- causal_forest(X.train, 
                            Y.train, 
                            W.train, 
                            clusters = cluster.train,
                            num.trees = 500)
    
    # Predict on within-fold data, rank within-fold
    tau.hat <- predict(forest, X.predict)$predictions
    tau.hat.quantiles <- quantile(tau.hat, probs = seq(0, 1, by=1/num.rankings))
    ranking[data.study$.folds == fold] <- cut(tau.hat, 
                                             tau.hat.quantiles, 
                                             include.lowest=TRUE,
                                             labels=seq(num.rankings))
    
    # Gather for output
    tau.hat.out[data.study$.folds == fold] <- predict(forest, X.predict)$predictions
  }
  
  data.out <- cbind.data.frame(data.study$casenumber,
                               data.study$pid,
                               data.study$year_orig,
                               tau.hat.out,
                               ranking)
  colnames(data.out) <- c("casenumber",
                          "pid",
                          "year_orig",
                          "tau.hat.out",
                          "ranking")
  
  return(data.out)
}

subgroup_plot <- function(data.prefix, outcome.selection, event.year){
  
  data.2009 <- get(paste0(data.prefix, ".2009"))
  data.2010 <- get(paste0(data.prefix, ".2010"))
  data.2011 <- get(paste0(data.prefix, ".2011"))
  data.2012 <- get(paste0(data.prefix, ".2012"))
  
  cf.2009 <- subgroup_cohort(data.2009, outcome.selection, event.year)
  cf.2009 <- left_join(cf.2009,
                       data.2009,
                       by=c("casenumber", "pid", "year_orig"))
  cf.2010 <- subgroup_cohort(data.2010, outcome.selection, event.year)
  cf.2010 <- left_join(cf.2010,
                       data.2010,
                       by=c("casenumber", "pid", "year_orig"))
  cf.2011 <- subgroup_cohort(data.2011, outcome.selection, event.year)
  cf.2011 <- left_join(cf.2011,
                       data.2011,
                       by=c("casenumber", "pid", "year_orig"))
  cf.2012 <- subgroup_cohort(data.2012, outcome.selection, event.year)
  cf.2012 <- left_join(cf.2012,
                       data.2012,
                       by=c("casenumber", "pid", "year_orig"))

  plotdata <- rbind.data.frame(cf.2009,
                               cf.2010,
                               cf.2011,
                               cf.2012)
  
  
  att <- rep(NA, 5)
  se <- rep(NA, 5)
  for(i in 1:5){
    
    current.2009 <- plotdata[plotdata$ranking==i & plotdata$year_beg_fs==2009,]
    current.2010 <- plotdata[plotdata$ranking==i & plotdata$year_beg_fs==2009,]
    current.2011 <- plotdata[plotdata$ranking==i & plotdata$year_beg_fs==2009,]
    current.2012 <- plotdata[plotdata$ranking==i & plotdata$year_beg_fs==2009,]
    
    est.2009 <- aipw(mu_est(current.2009, outcome.selection, event.year, "forest"), 
                       "e.hat.rf", outcome.selection)
    
    est.2010 <- aipw(mu_est(current.2010, outcome.selection, event.year, "forest"), 
                     "e.hat.rf", outcome.selection)
    
    est.2011 <- aipw(mu_est(current.2011, outcome.selection, event.year, "forest"), 
                     "e.hat.rf", outcome.selection)
    
    est.2012 <- aipw(mu_est(current.2012, outcome.selection, event.year, "forest"), 
                     "e.hat.rf", outcome.selection)
    
    rf.n <- est.2009[[1]] + est.2010[[1]] + est.2011[[1]] + est.2012[[1]]
    rf.j <- est.2009[[2]] + est.2010[[2]] + est.2011[[2]] + est.2012[[2]]
    
    att[i] <- ((est.2009[[1]]/rf.n)*est.2009[[3]]) + 
      ((est.2010[[1]]/rf.n)*est.2010[[3]]) +
      ((est.2011[[1]]/rf.n)*est.2011[[3]]) +
      ((est.2012[[1]]/rf.n)*est.2012[[3]])
    rf.var <- ((est.2009[[2]]/rf.j)^2*est.2009[[4]]) + 
      ((est.2010[[2]]/rf.j)^2*est.2010[[4]]) +
      ((est.2011[[2]]/rf.j)^2*est.2011[[4]]) +
      ((est.2012[[2]]/rf.j)^2*est.2012[[4]])
    se[i] <- sqrt(rf.var)

    #n.cases <- length(unique(data.now$casenumber))

    #se[i] <-  unique((data.now %>%
    #              mutate(scores = tau.hat.out) %>%
    #              group_by(casenumber) %>%
    #              dplyr::summarize(tauj=mean(scores)) %>%
    #              mutate(tau=mean(tauj)) %>%
    #              mutate(diff=(tau-tauj)^2) %>%
    #              mutate (se = sqrt(1/(n.cases*(n.cases-1)) * sum(diff))))$se)
    
    #att[i] <- mean(data.now$tau.hat.out)
  }
  
  output <- list(att = att,
                 se = se,
                 data = plotdata)
  
  return(output)
}

plot_subgroups_gg <- function(subgroup.output, title.text){
  
  att <- subgroup.output[[1]]
  print(att)
  se <- subgroup.output[[2]]
  print(se)
  
  CIL <- att - 1.96*se
  print(CIL)
  CIU <- att + 1.96*se
  print(CIU)
  
  ranking <- c(1, 2, 3, 4, 5)
  
  plot <- cbind.data.frame(ATT=att,
                           SE=se,
                           CIL=CIL,
                           CIU=CIU,
                           Rank=ranking)
  
  print(plot)
  
  ggplot(data=as.data.frame(plot), aes(x=Rank, y=ATT)) + 
    geom_point() + 
    geom_errorbar(aes(x=Rank, ymin=CIL, ymax=CIU)) + 
    theme_bw() + 
    ggtitle(paste("Data-Driven Subgroup Discovery:", title.text), 
            subtitle=paste("Cross-fit CATE Quintiles with Forest-AIPW Pooled Event Study Estimates")) 
}

moved.plot <- subgroup_plot("moved", "moved", 4)
divorce.plot <- subgroup_plot("divorce", "cumulnum_divorce", 4)
school.plot <- subgroup_plot("school", "schoolindex_high", 4)
unpaid.plot <- subgroup_plot("unpaid", "num00_collection_unpaid", 4)

pdf("./output/dd_groups_moved.pdf", width=8.67, height=5.75) 
plot_subgroups_gg(moved.plot, "Moved Indicator")
dev.off()


pdf("./output/dd_groups_divorce.pdf", width=8.67, height=5.75) 
plot_subgroups_gg(divorce.plot, "Cumulative Number of Divorces")
dev.off()

pdf("./output/dd_groups_school.pdf", width=8.67, height=5.75) 
plot_subgroups_gg(school.plot, "High School Test Score Index")
dev.off()

pdf("./output/dd_groups_unpaid.pdf", width=8.67, height=5.75) 
plot_subgroups_gg(unpaid.plot, "Number of Unpaid Collections")
dev.off()