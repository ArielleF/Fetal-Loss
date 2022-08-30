#!/usr/bin/env Rscript

# Script for running logistic regression models predicting fetal loss

# Load the data in Table S6
data <- read.csv("groom_anonymized.csv", header=T) 

# Load R libraries
library(glmmTMB)
library(dplyr)

#############################################################################################################################
# Model 1. Main logistic regression model predicting fetal loss.
#############################################################################################################################

# run model
m1 <- glmmTMB(, data=data, family="binomial")

# look at model output
summary(m1)
m1_betas <- fixef(m1)$cond # get betas (effect estimates) for intercept and predictor variables 
round(m1_betas, 3) # effect estimates reported in Table 1  


#############################################################################################################################
# Model 2. Logistic regression model predicting fetal loss using the most conservative inclusion criteria.
#############################################################################################################################

# Exclude females in non-wild feeding social groups, females with birthdates estimated with greater than a few days’ error, pregnancies that overlapped the 2009 hydrological year (a severe drought), pregnancies that occurred during several periods of reduced data collection, pregnancies that occurred during social group fissions and fusions, and pregnancies with conception and end dates that were not known to within a few days’ error
tmp <- 

nrow(tmp) # n=675 pregnancies
length(distinct(tmp, female)) # n=142 females

# run model (same as m1 above but on reduced data set)
m_s1 <- glmmTMB(, data=tmp, family="binomial")

# look at model output
summary(m_s1)
m_s1_betas <- fixef(m_s1)$cond # get betas (effect estimates) for intercept and predictor variables 
round(m_s1_betas, 3) # effect estimates reported in Table S1  

rm(tmp)

#############################################################################################################################
# Model 3. Logistic regression model predicting fetal loss using a mother’s proportional rank instead of ordinal rank at the time of conception as a predictor variable.
#############################################################################################################################

# run model (same as m1 but replace ordinal_rank with proportional_rank)
m_s2 <- glmmTMB(, data=data, family="binomial")

# look at model output
summary(m_s2)
m_s2_betas <- fixef(m_s2)$cond # get betas (effect estimates) for intercept and predictor variables 
round(m_s2_betas, 3) # effect estimates reported in Table S3


#############################################################################################################################
# Model 4. Logistic regression model predicting fetal loss excluding fetal losses due to feticide or presumed feticide.
#############################################################################################################################

# Exclude pregnancies that did not overlap the time period evaluated by Zipple et al. 2017 Proc Biol Sci, 284(1847) (doi:10.1098/rspb.2016.2561)
tmp <- 

nrow(tmp) # n=xxx pregnancies
length(distinct(tmp, female)) # n=xx females
summary(tmp$preg_outcome)
# 113 fetal losses total
summary(tmp$feticide)
# 8 feticides

# Exclude feticides or presumed feticides
tmp2 <- tmp[tmp

# run model (same as m1 above but on reduced data set)
m_s3 <- glmmTMB(, data=tmp2, family="binomial")

# look at model output
summary(m_s3)
m_s3_betas <- fixef(m_s3)$cond # get betas (effect estimates) for intercept and predictor variables 
round(m_s3_betas, 3) # effect estimates reported in Table S4

rm(tmp, tmp2)

#############################################################################################################################
# Model 5. Logistic regression model predicting fetal loss substituting hybrid status for genome-wide ancestry as a predictor variable.
#############################################################################################################################

# Exclude females where we could not assign a hybrid status
tmp <- 

nrow(tmp) # n=944 pregnancies
length(distinct(tmp, female)) # n=164 females

# run model (same as m1 above but replace genome_wide_anubis_ancestry and genome_wide_anubis_ancestry^2 with hybrid_status)
m_s4 <- glmmTMB(, data=tmp, family="binomial")

# look at model output
summary(m_s4)
m_s4_betas <- fixef(m_s4)$cond # get betas (effect estimates) for intercept and predictor variables 
round(m_s4_betas, 3) # effect estimates reported in Table S5

