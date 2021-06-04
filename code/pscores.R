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

# Raw Data ----------------------------------------------------------------
data <- fread("./data/analysis_data_owners.csv")

# 3 years previous to treatment, event-time
data.m3 <- data[data$evt==-3]

# Pscores: DGT Spec Replication -------------------------------------
covariates.unconfoundedness <- c("moved", "isowner", "log_finishedsqft_dq", "zipinc_log",
                                 "schoolindex_elem", "schoolindex_middle", "schoolindex_high",
                                 "cumulnum_divorce", "cumulnum_crime", "cumulnum_dui", 
                                 "cumulnum_bankruptcy", "vantage", "num00_foreclosure", 
                                 "num00_collection_unpaid", "num00_trade_auto", "num12_trade30", 
                                 "num00_mortgage_loanmod", "num00_mortgage_open", 
                                 "rtot00_mortgage_open_bal", "rtot00_mortgage_monthpay")
treatment <- "fc_al3"
fmla.pscore.DGT.selection <- formula(paste(treatment, "~", 
                                 paste(covariates.unconfoundedness, collapse="+"),
                                 " | zip_year+date_beg_fs"))

# Weights
data.m3 <- data.m3 %>%
  add_count(casenumber) %>%
  mutate(weight=1/n)

# Normalize
for(i in covariates.unconfoundedness){
  vector <- t(as.vector(dplyr::select(data.m3, all_of(i))))
  mean <- wtd.mean(vector, weight=data.m3$weight, na.rm=TRUE, normwt=FALSE)
  sdev <- sqrt(wtd.var(vector, weight=data.m3$weight, na.rm=TRUE, normwt=FALSE))
  vector <- vector-mean
  vector <- vector/sdev
  vector[is.na(vector)] <- 0
  data.m3[,i] <- vector
}

# Variable Selection, DGT
DGT.spec.pscore.pre <- feols(fmla.pscore.DGT.selection, 
                                 data.m3, 
                                 weight=data.m3$weight, 
                                 cluster=~casenumber)

summary(DGT.spec.pscore.pre)

# Variables selected (5 highest by t-stat)
covariates.pscore.DGT <- c("isowner", "log_finishedsqft_dq", "vantage", 
                            "num12_trade30", "num00_mortgage_open")

# FE regression
fmla.pscore.DGT <- formula(paste(treatment, "~", 
                                 paste(covariates.pscore.DGT, collapse="+"),
                                 " | zip_year+date_beg_fs"))

DGT.spec <- feols(fmla.pscore.DGT, 
                  data.m3, 
                  weight=data.m3$weight, 
                  cluster=~casenumber)

# Results
summary(DGT.spec)

# PScores
e.hat.DGT <- predict(DGT.spec, data.m3)
e.hat.DGT[e.hat.DGT>1] <- 1
e.hat.DGT[e.hat.DGT<0] <- 0

# PScores: Post-LASSO -----------------------------------------------------

# Training formula
train.fmla <- formula(paste(treatment, "~",
                            paste(covariates.unconfoundedness, collapse="+")))

# First stage lasso
# No intercept
XX <- model.matrix(train.fmla, 
                   data.m3)[,-1]
Y <- data.m3$fc_al3

lasso <- cv.glmnet(x=XX, 
                   y=Y, 
                   weights=data.m3$weight,  
                   alpha=1.) 

# Results
lambda.grid <- c(0, sort(lasso$lambda))
lasso.coefs <- as.matrix(coef(lasso, s=lambda.grid))
plot(lasso)
plot(lasso$glmnet.fit, xvar="lambda")

# With lambda
selected.lambda <- lasso$lambda.1se
n.folds <- 10
n <- nrow(data.m3)
foldid <- (seq(n) %% n.folds) + 1
coefs <- sapply(seq(n.folds), function(k) {
  lasso.fold <- glmnet(XX[foldid == k,], 
                       Y[foldid == k],
                       weights=data.m3[foldid == k,]$weight,
                       alpha=1.)
  as.matrix(coef(lasso.fold, s=selected.lambda))
})
heatmap(1*(coefs != 0), 
        Rowv = NA, 
        Colv = NA, 
        cexCol = 1, 
        scale="none", 
        col=gray(c(1,0)), 
        margins = c(3, 1), 
        xlab="Fold", 
        labRow=c("Intercept", covariates.unconfoundedness), 
        main="Non-zero coefficient estimates")

# Nonzero Coefs
## Returns nonzero coefs, not including intercept
covariates.postlasso <- rownames(coef(lasso, 
                                      s = 'lambda.1se'))[coef(lasso, s = 'lambda.1se')[,1]!= 0][-1] 


# DGT Regression Model w/ Post-Lasso Coefficients
fmla.pscore.postlasso <- formula(paste(treatment, "~", 
                                        paste(covariates.postlasso, collapse="+"),
                                        " | zip_year+date_beg_fs"))

DGT.spec.postlasso <- feols(fmla.pscore.postlasso, 
                            data.m3, 
                            weight=data.m3$weight, 
                            cluster=~casenumber)

summary(DGT.spec.postlasso)

# PScores
e.hat.postlasso <- predict(DGT.spec.postlasso, data.m3)
e.hat.postlasso[e.hat.postlasso>1] <- 1
e.hat.postlasso[e.hat.postlasso<0] <- 0


# PScores: Regression Forest w/ Clustering -------------------------------

# Training forest with cross fitting case-wise, reweighting
#train.fmla <- formula(paste(treatment, "~",
#                            paste(covariates.unconfoundedness, collapse="+"),
#                            "+ zip_year + date_beg_fs"))
#
#data.m3$zip_year <- as.factor(data.m3$zip_year)
#data.m3$date_beg_fs <- as.factor(data.m3$date_beg_fs)
#
#XX <- model.matrix(train.fmla, 
#                   data.m3)
#
#treat.vec <- data.m3$fc_al3
#clusters.vec <- as.factor(data.m3$casenumber)
#
#rf.model.raw <- regression_forest(X = XX,
#                                  Y = treat.vec,
#                                  num.trees = 500,
#                                  clusters = clusters.vec)
#
#varimp <- variable_importance(rf.model.raw)
#selected.idx <- which (varimp>mean (varimp))
#
#rf.model <- regression_forest(X = XX[,selected.idx],
#                              Y = treat.vec,
#                              num.trees=500,
#                              clusters = clusters.vec)
#
#e.hat.rf <- rf.model$predictions
#
#write.csv(e.hat.rf, "./rf_pscores.csv")

## Note, model chose to ignore ALL the fixed effects factors and focus on
## all of covariates.unconfoundedness. Command to check that was:
## colnames(XX[,selected.idx])

e.hat.rf <- fread("./data/rf_pscores.csv")
e.hat.rf <- e.hat.rf[,2]
e.hat.rf <- e.hat.rf$V1

# PScores: Random Tree w/ Clustering -------------------------------

# See Athey & Wager 2019, Abadie et al. 2017, my report Section 3.1
# Pscores on pre-period outcomes to match unconfoundedness condition
# Folding case-wise to account for correlations within clusters

# Folding case-wise (cases unique to folds)
#data.m3$casenumber <- as.factor(data.m3$casenumber)
#data.m3 <- ungroup(data.m3)
#data.m3 <- fold(data.m3,
#                k = 5,
#                id_col="casenumber")
#
#n <- nrow(data.m3)
#e.hat.tree <- rep(NA, n)
#
#train.fmla <- formula(paste(treatment, "~",
#                            paste(covariates.unconfoundedness, collapse="+")))
#
## Cross-fit p-scores, folding by case 
#for (idx in 1:5){
#  unpruned.tree <- rpart(train.fmla, 
#                         data = data.m3[data.m3$.folds!=idx,],
#                         cp=0,
#                         method="anova")
#  
#  cp.idx <- which.min(unpruned.tree$cptable[,"xerror"])
#  cp.best <- unpruned.tree$cptable[cp.idx, "CP"]
#  propensity.model <- prune(unpruned.tree, cp=cp.best)
#  e.hat.tree[data.m3$.folds==idx] <- predict(propensity.model, 
#                                             newdata=data.m3[data.m3$.folds==idx,], 
#                                             type="vector")
#}

e.hat.tree <- fread("./data/tree_pscores.csv")
e.hat.tree <- e.hat.tree$x

# PScores: Comparison (Histogram) -----------------------------------------------------
# Plot Comparison
pscore.type <- rep("Random Forest", 248363)
rf <- cbind(as.numeric(e.hat.rf), 
            pscore.type)

pscore.type <- rep("Draft Spec", 248363)
lm <- cbind(as.numeric(e.hat.DGT), 
            pscore.type)

pscore.type <- rep("Post-LASSO", 248363)
pl <- cbind(as.numeric(e.hat.postlasso), 
            pscore.type)

pscore.type <- rep("Random Tree", 248363)
tree <- cbind(as.numeric(e.hat.tree), 
                pscore.type)

# Bind
pscores <- as.data.frame(rbind(pl, 
                               lm, 
                               rf,
                               tree))

pscores$V1 <- as.numeric(pscores$V1)

pscores$Treatment <- as.character(data.m3$fc_al3)
pscores$Treatment <- ifelse(pscores$Treatment=="1", 
                            "Treated", 
                            "Not Treated")

pdf("./output/pscores_hist_comparison.pdf", width=8.67, height=5.75)  
ggplot(pscores, aes(x=V1, fill=Treatment, color=Treatment)) + 
  geom_histogram(position = "identity", alpha = 0.4) +
  facet_wrap(~pscore.type) +
  xlab("Propensity Score") +
  ylab("Count") + 
  theme_bw() + 
  ggtitle("Estimated Propensity Scores by Treatment Status")
dev.off()


# PScores: Comparison (Balance) --------------------------------------------

# PScore data with some event time = -1 data
pscore.data <- select(data.m3, 
                      pid, 
                      casenumber)
pscore.data$e.hat.rf <- e.hat.rf
pscore.data$e.hat.DGT <- e.hat.DGT
pscore.data$e.hat.postlasso <- e.hat.postlasso
pscore.data$e.hat.tree <- e.hat.tree

data.m1 <- data[data$evt==-1]
data.m1 <- left_join(data.m1, 
                     pscore.data,
                     by=c("pid", "casenumber"))

covariates.balance <- c("vantage", "age", "log_finishedsqft_dq", "zipinc_log", "schoolindex_high", "moved")
pscore.names <- c("e.hat.rf", "e.hat.DGT", "e.hat.postlasso", "e.hat.tree")

data.m1 <- select(data.m1,
                  all_of(covariates.balance),
                  all_of(pscore.names),
                  pid,
                  casenumber,
                  fc_al3)

# See Crump et al. 2009
data.m1 <- subset(data.m1, 
                  e.hat.rf >= 0.1 & e.hat.rf <= 0.9)
data.m1 <- subset(data.m1, 
                  e.hat.DGT >= 0.1 & e.hat.DGT <= 0.9)
data.m1 <- subset(data.m1, 
                  e.hat.postlasso >= 0.1 & e.hat.postlasso <= 0.9)
data.m1 <- subset(data.m1, 
                  e.hat.tree >= 0.1 & e.hat.tree <= 0.9)

# Formula
fmla.balance <- formula(paste("fc_al3 ~ 0 + ",
                              paste(covariates.balance, 
                                    collapse="+")))

XX <- model.matrix.lm(fmla.balance, data.m1, na.action=NULL)
W <- data.m1[,"fc_al3"]
pp <- ncol(XX)

# Unadjusted covariate means, variances and standardized abs mean differences
means.treat <- apply(XX[W == 1,], 2, mean, na.rm=TRUE)
means.ctrl <- apply(XX[W == 0,], 2, mean, na.rm=TRUE)
abs.mean.diff <- abs(means.treat - means.ctrl)

var.treat <- apply(XX[W == 1,], 2, var, na.rm=TRUE)
var.ctrl <- apply(XX[W == 0,], 2, var, na.rm=TRUE)
std <- sqrt(var.treat + var.ctrl)

# Adjusted covariate means, variances and standardized abs mean differences

### DGT ###
XX.treat <- as.matrix(as.data.frame(apply(XX,2, function(x) x*W/data.m1$e.hat.DGT)))
XX.control <- as.matrix(as.data.frame(apply(XX,2, function(x) x*(1-W)/(1-data.m1$e.hat.DGT))))
means.treat.adj <- apply(XX.treat, 2, mean, na.rm=TRUE)
means.ctrl.adj <- apply(XX.control, 2, mean, na.rm=TRUE)
abs.mean.diff.adj <- abs(means.treat.adj - means.ctrl.adj)

var.treat.adj <- apply(XX.treat, 2, var, na.rm=TRUE)
var.ctrl.adj <- apply(XX.control, 2, var, na.rm=TRUE)
std.adj <- sqrt(var.treat.adj + var.ctrl.adj)

### Post-LASSO ###
XX.treat2 <- as.matrix(as.data.frame(apply(XX,2, function(x) x*W/data.m1$e.hat.postlasso)))
XX.control2 <- as.matrix(as.data.frame(apply(XX,2, function(x) x*(1-W)/(1-data.m1$e.hat.postlasso))))
means.treat.adj2 <- apply(XX.treat2, 2, mean, na.rm=TRUE)
means.ctrl.adj2 <- apply(XX.control2, 2, mean, na.rm=TRUE)
abs.mean.diff.adj2 <- abs(means.treat.adj2 - means.ctrl.adj2)

var.treat.adj2 <- apply(XX.treat2, 2, var, na.rm=TRUE)
var.ctrl.adj2 <- apply(XX.control2, 2, var, na.rm=TRUE)
std.adj2 <- sqrt(var.treat.adj2 + var.ctrl.adj2)

### Neural Net ###
XX.treat3 <- as.matrix(as.data.frame(apply(XX,2, function(x) x*W/data.m1$e.hat.tree)))
XX.control3 <- as.matrix(as.data.frame(apply(XX,2, function(x) x*(1-W)/(1-data.m1$e.hat.tree))))
means.treat.adj3 <- apply(XX.treat3, 2, mean, na.rm=TRUE)
means.ctrl.adj3 <- apply(XX.control3, 2, mean, na.rm=TRUE)
abs.mean.diff.adj3 <- abs(means.treat.adj3 - means.ctrl.adj3)

var.treat.adj3 <- apply(XX.treat3, 2, var, na.rm=TRUE)
var.ctrl.adj3 <- apply(XX.control3, 2, var, na.rm=TRUE)
std.adj3 <- sqrt(var.treat.adj3 + var.ctrl.adj3)

### Random Forest ###
XX.treat4 <- as.matrix(as.data.frame(apply(XX,2, function(x) x*W/data.m1$e.hat.rf)))
XX.control4 <- as.matrix(as.data.frame(apply(XX,2, function(x) x*(1-W)/(1-data.m1$e.hat.rf))))
means.treat.adj4 <- apply(XX.treat4, 2, mean, na.rm=TRUE)
means.ctrl.adj4 <- apply(XX.control4, 2, mean, na.rm=TRUE)
abs.mean.diff.adj4 <- abs(means.treat.adj4 - means.ctrl.adj4)

var.treat.adj4 <- apply(XX.treat4, 2, var, na.rm=TRUE)
var.ctrl.adj4 <- apply(XX.control4, 2, var, na.rm=TRUE)
std.adj4 <- sqrt(var.treat.adj4 + var.ctrl.adj4)

# Plot
pdf("./output/pscores_balance_comparison.pdf", width=8.67, height=5.75)  
par(oma=c(0,4,0,0))
plot(-2, xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(-.01, 0.1), ylim=c(0, pp+1), 
     main="Propensity Score Model Comparison",
     sub="Standardized Absolute Mean Differences One Year Pre-Treatment")
axis(side=1, at=c(-1, 0, 0.05, 0.1, 0.15), las=1)
lines(abs.mean.diff / std, seq(1, pp), type="p", col="blue", pch=20)
lines(abs.mean.diff.adj / std.adj, seq(1, pp), type="p", col="purple", pch=19)
lines(abs.mean.diff.adj2 / std.adj2, seq(1, pp), type="p", col="green", pch=20)
lines(abs.mean.diff.adj3 / std.adj3, seq(1, pp), type="p", col="orange", pch=20)
lines(abs.mean.diff.adj4 / std.adj3, seq(1, pp), type="p", col="red", pch=20)
legend("topright", c("Unadjusted", "DGT", "Post-Lasso", "Random Tree", "Random Forest"), col=c("blue", "purple", "green", "orange", "red"), pch=19, ncol=3, bty='n', cex=0.8)
abline(v = seq(0, 1, by=.25), lty = 2, col = "grey", lwd=.5)
abline(h = 1:pp,  lty = 2, col = "grey", lwd=.5)
mtext(colnames(XX), side=2, cex=0.7, at=1:pp, padj=.4, adj=1, col="black", las=1, line=.3)
abline(v = 0)
dev.off()

# PScores: Comparison (Correlation) -------------------------------------------------------

# Correlation matrix (6 entries!)
# Scatter Plots 
# Decile mean comparisons

data.m1 <- data[data$evt==-1]
data.m1 <- left_join(data.m1, 
                     pscore.data,
                     by=c("pid", "casenumber"))

# See Crump et al. 2009
data.m1 <- data.m1 %>%
           mutate(tag = ifelse(e.hat.rf >= 0.1 & e.hat.rf <= 0.9 &
                                 e.hat.DGT >= 0.1 & e.hat.DGT <= 0.9 & 
                                 e.hat.postlasso >= 0.1 & e.hat.postlasso <= 0.9 & 
                                 e.hat.tree >= 0.1 & e.hat.tree <= 0.9, 0, 1))

tags <- data.m1$tag
data.m1 <- subset(data.m1, tag==0)

covariates.balance <- c("vantage", "age", "log_finishedsqft_dq", "zipinc_log", "schoolindex_high", "moved")
pscores.all <- cbind(e.hat.DGT, e.hat.tree, e.hat.postlasso, e.hat.rf)

# Correlation Analysis
#pdf("./output/pscores_corr_comparison.pdf", width=8.67, height=5.75) 
#pairs.panels(pscores.all,
#             show.points=FALSE,
#             method="spearman",
#             )
#dev.off()


# Imbens-Rubin Recursive Stratification -----------------------------------

# See Imbens & Rubun 2015, Chapter 17

imbens_recurse <- function(data, scores, treatment) {
  # Get scores into data
  data[, "pscore"] <- scores
  data[, "treatment"] <- select(data, all_of(treatment))
  
  # Decide on Propensity Score Splits, starting from only one bin
  j = 1
  p.bins = c(0,1)
  n.strata = 1
  check.me = 0
  
  # Loop over a condition for when all bins fit criterion
  while(j==1){
    idx <- seq(1, n.strata)
    bin.avg.prop <- rep(NA, n.strata)
    bin.t <- rep(NA, n.strata)
    bin.n <- rep(NA, n.strata)
    
    # Looping over current stratum definition
    for(i in idx){
      
      # Find details of this stratum
      stratum.min <- p.bins[i]
      stratum.max <- p.bins[1 + i]
      stratum.contents <- subset(data, stratum.min<data$pscore & data$pscore<=stratum.max)
      stratum.treat <- subset(stratum.contents, stratum.contents$treatment==1)
      stratum.control <- subset(stratum.contents, stratum.contents$treatment==0)
      
      # Calculate t-stat for the linear difference in propensity score means
      # between treatment and control groups inside the bin
      
      if(max(stratum.contents$pscore)-min(stratum.contents$pscore) > 0.01){
        stratum.t <- t.test(pscore ~ fc_al3, data=stratum.contents)
      }
      else{
        stratum.t <- 0
      }
      
      
      # Note some summary stats for the bin
      bin.avg.prop[i] <- mean(stratum.contents$pscore)
      bin.t[i] <- stratum.t[[1]]
      bin.n[i] <- nrow(stratum.contents)
      
      # If difference in scores is too high, and bin is large enough, 
      # split and loop over all bins again
      if(abs(stratum.t[[1]]) > 1.96 & 
         nrow(subset(stratum.treat, stratum.treat$pscore<=median(stratum.contents$pscore)))>3 &
         nrow(subset(stratum.treat, stratum.treat$pscore>median(stratum.contents$pscore)))>3 &
         nrow(subset(stratum.control, stratum.control$pscore<=median(stratum.contents$pscore)))>3 &
         nrow(subset(stratum.control, stratum.control$pscore>median(stratum.contents$pscore)))>3 & 
         nrow(subset(stratum.contents, stratum.contents$pscore<=median(stratum.contents$pscore)))>23 & 
         nrow(subset(stratum.contents, stratum.contents$pscore>median(stratum.contents$pscore)))>23){
        
        n.strata <- n.strata + 1
        new.val <- median(stratum.contents$pscore)
        p.bins <- append(p.bins, new.val, after = i)
        check.me = 0
        j=1
      }
      
      # If this bin is fine, note that and go to next bin
      else{
        check.me = check.me + 1
      }
      
      # If this bin is fine and all previous bins look fine,
      # end
      if(i == n.strata & check.me >= n.strata){
        j = 0
      }
    }
  }
  
  # Return: number of strata, cut points for the bins, 
  # average propensity score within strata, t-stats for each stratum,
  # and the number of obs. in each stratum. 
  output.imbens <- list(n.strata, p.bins, bin.avg.prop, bin.t, bin.n)
  return(output.imbens)
}


# PScores: Comparison (Boxplots) -------------------------------------------

# Decile Balance Comparisons
## Contruct Deciles
for(p in pscore.names){
  dec.name <- paste0("dec.",p)
  decile <- as.numeric(cut2(get(p)[tags==0], g=10))
  data.m1[,all_of(dec.name)] <- decile
}

## Group by Imbens-Rubin Instead
imbens.rf <- imbens_recurse(data.m1, e.hat.rf[tags==0], "fc_al3")
imbens.postlasso <- imbens_recurse(data.m1, e.hat.postlasso[tags==0], "fc_al3")
imbens.DGT <- imbens_recurse(data.m1, e.hat.DGT[tags==0], "fc_al3")
imbens.tree <- imbens_recurse(data.m1, e.hat.tree[tags==0], "fc_al3")

for(p in pscore.names){
  suffix <- str_extract(p, "(?<=.)[a-zA-Z]*$")
  imb.name <- paste0("imbens.",suffix)
  groups <- as.numeric(cut(get(p)[tags==0], 
                           breaks = c(get(imb.name)[[2]])))
  
  data.m1[,all_of(imb.name)] <- groups
}


# For group-variable, get t-stat for each bin 
var_score_t <- function(data, variable, strata){
  data$vector <- select(data, all_of(strata))
  loopto <- max(data$vector)
  bin.t <- rep(NA, loopto)
  for(q in 1:loopto){
    stratum.contents <- subset(data, data$vector==q) %>%
      select(all_of(variable), fc_al3) %>%
      rename(variable = all_of(variable))
    
    stratum.contents <- subset(stratum.contents, 
                               complete.cases(stratum.contents))
    
    bin.t[q] <- t.test(variable ~ fc_al3, data=stratum.contents)[[1]]
  }
  return(as.data.frame(bin.t))
}

# Above for every group-variable pair of interest
strata_tstat_box <- function(data, variable){
  hold.rf <- var_score_t(data, variable, "imbens.rf")
  hold.tree <- var_score_t(data, variable, "imbens.tree")
  hold.postlasso <- var_score_t(data, variable, "imbens.postlasso")
  hold.DGT <- var_score_t(data, variable, "imbens.DGT")
  
  hold.rf$Bin <- rownames(hold.rf)
  hold.rf$Method <- "Random Forest"
  
  hold.tree$Bin <- rownames(hold.tree)
  hold.tree$Method <- "Random Tree"
  
  hold.postlasso$Bin <- rownames(hold.postlasso)
  hold.postlasso$Method <- "Post-LASSO"
  
  hold.DGT$Bin <- rownames(hold.DGT)
  hold.DGT$Method <- "DGT Spec"
  
  output <- rbind(hold.rf, hold.tree, hold.postlasso, hold.DGT)
  return(output)
}

# Relevant variables
t.moved <- strata_tstat_box(data.m1, "moved")
t.moved$Variable <- "Moved"
t.school_high <- strata_tstat_box(data.m1, "schoolindex_high")
t.school_high$Variable <- "High School Index"
t.zipinc_log <- strata_tstat_box(data.m1, "zipinc_log")
t.zipinc_log$Variable <- "Zip Income"
t.log_finished <- strata_tstat_box(data.m1, "log_finishedsqft_dq")
t.log_finished$Variable <- "Square Feet"
t.age <- strata_tstat_box(data.m1, "age")
t.age$Variable <- "Age"
t.vantage <- strata_tstat_box(data.m1, "vantage")
t.vantage$Variable <- "Vantage Score"

# Bind
data.tboxes <- rbind(t.moved,
                     t.school_high,
                     t.zipinc_log,
                     t.log_finished,
                     t.age,
                     t.vantage)

# MSE
data.tboxes <- data.tboxes %>%
               group_by(Method) %>%
               mutate(mse = mean(bin.t^2, na.rm=TRUE)) %>%
               ungroup

MSE.strata <- as.data.frame(unique(data.tboxes$mse))
rownames(MSE.strata) <- c("RF", "Tree", "Post-LASSO", "DGT")
            

# Plot: Strata Plots
pdf("./output/pscores_imbens_boxplot_comparison.pdf", width=8.67, height=5.75)  
ggplot(data.tboxes, aes(Method, bin.t)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, aes(colour=Variable)) + 
  geom_hline(yintercept = 1.96, linetype="dashed", alpha=0.4, color="black") + 
  geom_hline(yintercept = -1.96, linetype="dashed", alpha=0.4, color="black") + 
  xlab("Method") +
  ylab("Stratum T-Stat") + 
  theme_bw() + 
  ggtitle("Stratification Balance Checks", subtitle="Imbens-Rubin Stratification Bins")
dev.off()

decile_tstat_box <- function(variable){
  
  # Takes a variable as input, contructs dataframe with t-stats, plots
  bin.t <- as.data.frame(matrix(, ncol=4, nrow=10))
  colnames(bin.t) <- c("Random Forest", "DGT Spec", "Post-LASSO", "Random Tree")

  for(j in 1:4){
    pscore.now <- pscore.names[j]
    dec.name <- paste0("dec.",pscore.now)

    for(i in 1:10){
      
      # Decile contents
      data.m3 <- as.data.frame(data.m1)
      stratum.contents <- subset(data.m1, get(dec.name)==i) %>%
                          select(all_of(variable), fc_al3) %>%
                          rename(variable = all_of(variable))
      
      stratum.contents <- subset(stratum.contents, 
                                 complete.cases(stratum.contents))
      
      # Return to the matrix
      bin.t[i,j] <- t.test(variable ~ fc_al3, data=stratum.contents)[1]

    }
  }
  
  # Coerce For Plot
  bin.t[,"Decile"] <- rownames(bin.t)
  bin.t[,"Variable"] <- variable
  
  bin.t <- melt(as.data.table(bin.t),
                id.vars = c("Decile", "Variable"),
                variable.name = "e.hat")
  
  return(bin.t)
}

# Variables of Interest
t.moved <- decile_tstat_box("moved")
t.school_high <- decile_tstat_box("schoolindex_high")
t.zipinc_log <- decile_tstat_box("zipinc_log")
t.log_finished <- decile_tstat_box("log_finishedsqft_dq")
t.age <- decile_tstat_box("age")
t.vantage <- decile_tstat_box("vantage")

# Bind
data.tboxes <- rbind(t.moved,
                     t.school_high,
                     t.zipinc_log,
                     t.log_finished,
                     t.age,
                     t.vantage)

data.tboxes <- data.tboxes %>%
               group_by(e.hat) %>%
               mutate(mse = mean(value^2)) %>%
               ungroup()

MSE.decile <- as.data.frame(unique(data.tboxes$mse))
rownames(MSE.decile) <- c("RF", "DGT", "Post-LASSO", "Tree")


# Plot: Decile Cuts
pdf("./output/pscores_decile_boxplot_comparison.pdf", width=8.67, height=5.75)
ggplot(data.tboxes, aes(e.hat, value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, aes(colour=Variable)) + 
  geom_hline(yintercept = 1.96, linetype="dashed", alpha=0.4, color="black") + 
  geom_hline(yintercept = -1.96, linetype="dashed", alpha=0.4, color="black") +
  xlab("Method") +
  ylab("Decile T-Stat") + 
  theme_bw() + 
  ggtitle("Stratification Balance Checks", subtitle="Decile Stratification Bins")
dev.off()

## Some summary stats for inspection and reporting
## Looks like Imbens-Rubin DGT has done well because it can be easily split
MSE.strata
MSE.decile
max(data.m1$imbens.DGT)
max(data.m1$imbens.rf)
max(data.m1$imbens.postlasso)
max(data.m1$imbens.tree)

# DGT Spec Robustness -----------------------------------------------------

# Run rebecca's specification for each propensity score method, on:
# moved, school test rank high, divorce, num00_collection_unpaid
# With both decile cuts and Imbens-Rubin cuts. 

# Main Data
data <- fread("./data/analysis_data_owners.csv")

# Get propensity scores to merge to analysis data
merge.pscores <- select(data.m3,
                        casenumber,
                        pid)

merge.pscores$e.hat.DGT <- e.hat.DGT
merge.pscores$e.hat.rf <- e.hat.rf
merge.pscores$e.hat.tree <- e.hat.tree
merge.pscores$e.hat.postlasso <- e.hat.postlasso

# Deciles and Stratification Bins (No Pscore Slicing)
merge.pscores$dec.e.hat.rf <- as.numeric(cut2(e.hat.rf, g=10))
merge.pscores$dec.e.hat.DGT <- as.numeric(cut2(e.hat.DGT, g=10))
merge.pscores$dec.e.hat.postlasso <- as.numeric(cut2(e.hat.postlasso, g=10))
merge.pscores$dec.e.hat.tree <- as.numeric(cut2(e.hat.tree, g=10))

imbens.rf <- imbens_recurse(data.m3, e.hat.rf, "fc_al3")
imbens.postlasso <- imbens_recurse(data.m3, e.hat.postlasso, "fc_al3")
imbens.DGT <- imbens_recurse(data.m3, e.hat.DGT, "fc_al3")
imbens.tree <- imbens_recurse(data.m3, e.hat.tree, "fc_al3")

merge.pscores$imbens.e.hat.rf <- as.numeric(cut(merge.pscores$e.hat.rf, 
                                                breaks = imbens.rf[[2]],
                                                include.lowest=TRUE))
merge.pscores$imbens.e.hat.postlasso <- as.numeric(cut(merge.pscores$e.hat.postlasso, 
                                                       breaks = imbens.postlasso[[2]],
                                                       include.lowest=TRUE))
merge.pscores$imbens.e.hat.DGT <- as.numeric(cut(merge.pscores$e.hat.DGT, 
                                                 breaks = imbens.DGT[[2]],
                                                 include.lowest=TRUE))
merge.pscores$imbens.e.hat.tree <- as.numeric(cut(merge.pscores$e.hat.tree, 
                                                  breaks = imbens.tree[[2]],
                                                  include.lowest=TRUE))

# Analysis data for pooled spec is pre-period, reference 0 period,
# and periods 3 + 4 post. See DGT Table 7
data.analysis <- subset(data,
                        data$evt<=0 | data$evt==3 | data$evt==4)

data.analysis <- left_join(data.analysis,
                           merge.pscores,
                           by = c("pid", "casenumber"))

## Indicators
# Event pre treatment
data.analysis <- data.analysis %>%
  mutate(fc_evt_pre=ifelse(evt<0, fc_al3, 0))

# Event post treatment
data.analysis <- data.analysis %>%
  mutate(fc_evt_post=ifelse(evt>0, fc_al3, 0))

# Event*Decile*Date
data.analysis$date_beg_fs <- as.factor(data.analysis$date_beg_fs)

data.analysis <- data.analysis %>%
  mutate(evt.date.dec.DGT=interaction(evt,
                                      dec.e.hat.DGT,
                                      date_beg_fs))

data.analysis <- data.analysis %>%
  mutate(evt.date.dec.rf=interaction(evt,
                                     dec.e.hat.rf,
                                     date_beg_fs))

data.analysis <- data.analysis %>%
  mutate(evt.date.dec.postlasso=interaction(evt,
                                            dec.e.hat.postlasso,
                                            date_beg_fs))

data.analysis <- data.analysis %>%
  mutate(evt.date.dec.tree=interaction(evt,
                                       dec.e.hat.tree,
                                       date_beg_fs))

# Event*Decile*Year
data.analysis$zip_year <- as.factor(data.analysis$zip_year)

data.analysis <- data.analysis %>%
  mutate(evt.year.dec.DGT=interaction(evt,
                                      dec.e.hat.DGT,
                                      zip_year))

data.analysis <- data.analysis %>%
  mutate(evt.year.dec.rf=interaction(evt,
                                     dec.e.hat.rf,
                                      zip_year))

data.analysis <- data.analysis %>%
  mutate(evt.year.dec.postlasso=interaction(evt,
                                      dec.e.hat.postlasso,
                                      zip_year))

data.analysis <- data.analysis %>%
  mutate(evt.year.dec.tree=interaction(evt,
                                      dec.e.hat.tree,
                                      zip_year))

# Event*Strat*Date
data.analysis <- data.analysis %>%
  mutate(evt.date.strat.DGT=interaction(evt,
                                      imbens.e.hat.DGT,
                                      date_beg_fs))

data.analysis <- data.analysis %>%
  mutate(evt.date.strat.rf=interaction(evt,
                                     imbens.e.hat.rf,
                                     date_beg_fs))

data.analysis <- data.analysis %>%
  mutate(evt.date.strat.postlasso=interaction(evt,
                                            imbens.e.hat.postlasso,
                                            date_beg_fs))

data.analysis <- data.analysis %>%
  mutate(evt.date.strat.tree=interaction(evt,
                                       imbens.e.hat.tree,
                                       date_beg_fs))

# Event*Strat*Year
data.analysis <- data.analysis %>%
  mutate(evt.year.strat.DGT=interaction(evt,
                                        imbens.e.hat.DGT,
                                      zip_year))

data.analysis <- data.analysis %>%
  mutate(evt.year.strat.rf=interaction(evt,
                                       imbens.e.hat.rf,
                                     zip_year))

data.analysis <- data.analysis %>%
  mutate(evt.year.strat.postlasso=interaction(evt,
                                              imbens.e.hat.postlasso,
                                            zip_year))

data.analysis <- data.analysis %>%
  mutate(evt.year.strat.tree=interaction(evt,
                                         imbens.e.hat.tree,
                                       zip_year))

## Regressions
psm.fmla <- formula(paste("moved ~ fc_evt_pre+fc_evt_post | evt+ fc_al3+ pid+ evt.match.date+ evt.match.year"))

match.date.dec <- c("evt.date.dec.DGT", "evt.date.dec.rf", "evt.date.dec.postlasso", "evt.date.dec.tree")
match.year.dec <- c("evt.year.dec.DGT", "evt.year.dec.rf", "evt.year.dec.postlasso", "evt.year.dec.tree")
match.date.strat <- c("evt.date.strat.DGT", "evt.date.strat.rf", "evt.date.strat.postlasso", "evt.date.strat.tree")
match.year.strat <- c("evt.year.strat.DGT", "evt.year.strat.rf", "evt.year.strat.postlasso", "evt.year.strat.tree")

# Weights
data.analysis <- data.analysis %>%
  add_count(casenumber) %>%
  mutate(weight=1/n)


dgt.robustness <- function(variable){
  
  # Variable as dep var
  psm.fmla <- formula(paste0(variable, " ~ fc_evt_pre+fc_evt_post | evt+ fc_al3+ pid+ evt.match.date+ evt.match.year"))
  
  for(i in 1:4){
    
    # Get the correct fixed effects
    date.dec <- match.date.dec[i]
    year.dec <- match.year.dec[i]
    date.strat <- match.date.strat[i]
    year.strat <- match.year.strat[i]
    
    # Using these ones: deciles
    print(date.dec)
    print(year.dec)
    
    data.analysis$evt.match.date <- select(data.analysis, all_of(date.dec))
    data.analysis$evt.match.year <- select(data.analysis, all_of(year.dec))
    
    print(feols(psm.fmla, data.analysis, weight=data.analysis$weight, cluster=~casenumber))
    
    # Using these ones: Imbens-Rubin
    print(date.strat)
    print(year.strat)
    
    data.analysis$evt.match.date <- select(data.analysis, all_of(date.strat))
    data.analysis$evt.match.year <- select(data.analysis, all_of(year.strat))
    
    print(feols(psm.fmla, data.analysis, weight=data.analysis$weight, cluster=~casenumber))
    
  }
  
}

dgt.robustness("moved")
dgt.robustness("schoolindex_high")
dgt.robustness("cumulnum_divorce")
dgt.robustness("num00_collection_unpaid")
