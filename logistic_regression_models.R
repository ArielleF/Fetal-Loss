#!/usr/bin/env Rscript

# Script for running logistic regression models predicting fetal loss

# Load the data in Table S4
data <- read.csv("TableS4.csv", header=T) 

# Load R libraries
library(glmmTMB)
library(dplyr)

#############################################################################################################################
# Model 1. Main logistic regression model predicting fetal loss.
#############################################################################################################################

# code pregnancy outcome (live birth or fetal loss) as a binary variable where fetal loss=1 and live birth=0
data$pregnancy_outcome2 <- as.factor(ifelse(data$pregnancy_outcome=="live_birth", 0, 1))

# run model (note that I(variable^2) codes for the quadratic transformation of the variable in the model)
m1 <- glmmTMB(pregnancy_outcome2 ~ genome_wide_anubis_ancestry + I(genome_wide_anubis_ancestry^2)  +  female_age_conception + I(female_age_conception^2) + previous_fetal_losses + female_ordinal_rank_conception + group_size_females_conception + temp_max_2months_conception + temp_max_2months_endpreg +   rainfall_5months_conception + rainfall_5months_endpreg +  habitat_quality + (1|female_id), family=binomial, data=data)

# look at model output
summary(m1)
write.csv(round(summary(m1)$coefficients$cond, 3), "table_1_results.csv") # results reported in Table 1  


#############################################################################################################################
# Model 2. Logistic regression model predicting fetal loss using the most conservative inclusion criteria.
#############################################################################################################################

# Exclude females in non-wild feeding social groups, females with birthdates estimated with greater than a few days’ error, pregnancies that overlapped the 2009 hydrological year (a severe drought), pregnancies that occurred during several periods of reduced data collection, pregnancies that occurred during social group fissions and fusions, and pregnancies with conception and end dates that were not known to within a few days’ error
# Exclude the following pregnancies:
# females in non-wild feeding social groups
# females with birthdates estimated with greater than a few days’ error
# pregnancies that overlapped the 2009 hydrological year (a severe drought)
# pregnancies that occurred during several periods of reduced data collection
# pregnancies that occurred during social group fissions and fusions
# pregnancies with conception and end dates that were not known to within a few days’ error
tmp <- data[data$social_group_wild_vs_provisioned_conception!="provisioned" & data$social_group_wild_vs_provisioned_endpreg!="provisioned" & data$female_age_uncertainty=="few_days" & data$drought_2009!="drought2009" &
data$reduced_data_collected_conception=="No" &
data$reduced_data_collected_endpreg=="No" &
data$social_group_fission_fusion_conception=="No" & 
data$social_group_fission_fusion_endpreg=="No" &
data$conception_date_uncertainty=="No" & data$pregnancy_date_uncertainty=="No",]

nrow(tmp) # n=676 pregnancies
nrow(distinct(tmp, female_id)) # n=142 unique females

# run model (same as m1 above but on reduced data set)
m_s1 <-  glmmTMB(pregnancy_outcome2 ~ genome_wide_anubis_ancestry + I(genome_wide_anubis_ancestry^2)  +  female_age_conception + I(female_age_conception^2) + previous_fetal_losses + female_ordinal_rank_conception + group_size_females_conception + temp_max_2months_conception + temp_max_2months_endpreg +   rainfall_5months_conception + rainfall_5months_endpreg +  habitat_quality + (1|female_id), family=binomial, data=tmp)

# look at model output
summary(m_s1)
write.csv(round(summary(m_s1)$coefficients$cond, 3), "table_s1_results.csv") # results reported in Table S1  

rm(tmp)

#############################################################################################################################
# Model 3. Logistic regression model predicting fetal loss using a mother’s proportional rank instead of ordinal rank at the time of conception as a predictor variable.
#############################################################################################################################

# run model (same as m1 but replace ordinal_rank with proportional_rank)
m_s2 <- glmmTMB(pregnancy_outcome2 ~ genome_wide_anubis_ancestry + I(genome_wide_anubis_ancestry^2)  +  female_age_conception + I(female_age_conception^2) + previous_fetal_losses + female_proportional_rank_conception + group_size_females_conception + temp_max_2months_conception + temp_max_2months_endpreg +   rainfall_5months_conception + rainfall_5months_endpreg +  habitat_quality + (1|female_id), family=binomial, data=data)

# look at model output
summary(m_s2)
write.csv(round(summary(m_s2)$coefficients$cond, 3), "table_s3_results.csv") # results reported in Table S3  


#############################################################################################################################
# Model 4. Logistic regression model predicting fetal loss excluding fetal losses due to feticide or presumed feticide.
#############################################################################################################################

# Exclude pregnancies that were high confidence or possible feticides based on Zipple et al. 2017 Proc Biol Sci, 284(1847) (doi:10.1098/rspb.2016.2561)
tmp <- data[data$high_confidence_feticide!="Yes" & data$possible_feticide!="Yes",]

nrow(tmp) # n=1,012 pregnancies
nrow(distinct(tmp, female_id)) # n=175 unique females
summary(tmp$pregnancy_outcome)

# run model (same as m1 above but on reduced data set)
m_s3 <- glmmTMB(pregnancy_outcome2 ~ genome_wide_anubis_ancestry + I(genome_wide_anubis_ancestry^2)  +  female_age_conception + I(female_age_conception^2) + previous_fetal_losses + female_ordinal_rank_conception + group_size_females_conception + temp_max_2months_conception + temp_max_2months_endpreg +   rainfall_5months_conception + rainfall_5months_endpreg +  habitat_quality + (1|female_id), family=binomial, data=tmp)

# look at model output
summary(m_s3)
write.csv(round(summary(m_s3)$coefficients$cond, 3), "table_s5_results.csv") # results reported in Table S5  

rm(tmp)

#############################################################################################################################
# Model 5. Logistic regression model predicting fetal loss substituting hybrid status for genome-wide ancestry as a predictor variable.
#############################################################################################################################

# Exclude females where we could not assign a hybrid status
tmp <- data[data$hybrid_status!="unassigned",]

nrow(tmp) # n=944 pregnancies
nrow(distinct(tmp, female_id)) # n=164 unique females

# run model (same as m1 above but replace genome_wide_anubis_ancestry and genome_wide_anubis_ancestry^2 with hybrid_status)
m_s4 <- glmmTMB(pregnancy_outcome2 ~ hybrid_status  +  female_age_conception + I(female_age_conception^2) + previous_fetal_losses + female_ordinal_rank_conception + group_size_females_conception + temp_max_2months_conception + temp_max_2months_endpreg +   rainfall_5months_conception + rainfall_5months_endpreg +  habitat_quality + (1|female_id), family=binomial, data=tmp)

# look at model output
summary(m_s4)
write.csv(round(summary(m_s4)$coefficients$cond, 3), "table_s6_results.csv") # results reported in Table S6 
            
rm(tmp)

