#!/usr/bin/env Rscript

# Scripts for producing Figures 1-3 and S1

# Load the data in Table S6
data <- read.csv("TableS6.csv", header=T) 

# Load R libraries
library(ggplot2)
library(glmmTMB)
library(dplyr)

set.seed(1234) # so figures with any jittering are reproducible - will not completely match the published figures because although this set.seed was used to generate the paper's figures, the figures were produced in a different order below (order below reflects order of figures in the paper)

#############################################################################################################################
# Figure 1. Live birth and fetal loss patterns in Amboseli baboon females.
#############################################################################################################################



# Plot figure 1
ggplot()
ggsave("fig1.png")

#############################################################################################################################
# Figure 2. Records of fetal loss by trimester in the data set.
#############################################################################################################################

# Plot figure 2
ggplot()
ggsave("fig2.png")

#############################################################################################################################
# Figure 3. Records of fetal loss by trimester in the data set.
#############################################################################################################################

# Need to run original logistic regression model for some of the figure making
# code pregnancy outcome (live birth or fetal loss) as a binary variable where fetal loss=1 and live birth=0
data$pregnancy_outcome2 <- as.factor(ifelse(data$pregnancy_outcome=="live_birth", 0, 1))

####################
# Fig. 3A
####################

# Plot figure 3A
ggplot()
ggsave("fig3A.png")

####################
# Fig. 3B
####################
# For the dashed line which shows the predicted relationship between male dominance rank and grooming based on model estimates, assuming average values for all other covariates
# Get the grooming probability for male rank given an average female social partner in an average demographic environment
# Note that we exclude male rank and female id and male id are coded as NA
groomm_r <- data.frame(females_in_group=mean(groom$females_in_group), males_in_group=mean(groom$males_in_group), rank_female=mean(groom$rank_female), genetic_ancestry_female=mean(groom$genetic_ancestry_female), heterozygosity_male=mean(groom$heterozygosity_male), heterozygosity_female=mean(groom$heterozygosity_female), genetic_relatedness=mean(groom$genetic_relatedness), observer_effort=mean(groom$observer_effort), female_age=mean(groom$female_age), female_age_transformed=mean(groom$female_age_transformed), reproductive_state=1, pair_coresidency=mean(groom$pair_coresidency), genetic_ancestry_male=mean(groom$genetic_ancestry_male), assortative_genetic_ancestry_index=mean(groom$assortative_genetic_ancestry_index), female_id=NA, male_id=NA)
# Get the range of empirical values for male rank
m_rank <- as.data.frame(seq(min(groom$rank_male),max(groom$rank_male), by=1))
colnames(m_rank) <- c("rank_male")

tmp <- merge(groomm_r,m_rank) # combine dataframe of averaged variables with the dataframe of the empirical range of male rank

# Get the grooming probability for the range of male rank values assuming average values for all other covariates
tmp$probs <- predict(groom_model, tmp, type="response", re.form=NA)

# For the colored dots which show probabilities based on counts of grooming occurrences (same calculation as in figure 1A but per rank)
raw <- groom %>% group_by(rank_male) %>% mutate(count_groom=dplyr::n()) # count the total number of grooming opportunities (i.e., the total number of co-resident pairings) for each male rank
raw2 <- raw %>% group_by(rank_male, count_groom) %>% mutate(yes_groom=sum(groom_two_month)) # count the total number of grooming occurrences (i.e., 1 for groom_two_month) for each male rank
raw3 <- distinct(raw2, rank_male, count_groom, yes_groom)  # only grab one row per rank
raw3$groom_prox <- raw3$yes_groom/raw3$count_groom # divide the total number of grooming occurrences by the total number of grooming opportunities to get the probability of grooming for each male rank (without adjustment for other covariates)

# Plot figure 3B
ggplot()
ggsave("fig3B.png")

####################
# Fig. 3C
####################

# Plot figure 3C
ggplot()
ggsave("fig3C.png")

####################
# Fig. 3D
####################

# Plot figure 3D
ggplot()
ggsave("fig3D.png")

#############################################################################################################################
# Figure S1. The total number of pregnancies per female and its relationship to the number of different pregnancy outcomes.
#############################################################################################################################

# Plot figure S1
ggplot()
ggsave("figS1.png")
