#!/usr/bin/env Rscript

# Scripts for producing Figures 1-3 and S1

# Load the data in Table S4
data <- read.csv("TableS4.csv", header=T) 

# Load R libraries
library(ggplot2)
library(glmmTMB)
library(data.table)
library(dplyr)
library(patchwork)

set.seed(1234) # so figures with any jittering are reproducible - will not completely match the published figures because although this set.seed was used to generate the paper's figures, the figures were produced in a different order below (order below reflects order of figures in the paper)

#############################################################################################################################
# Figure 1. Live birth and fetal loss patterns in Amboseli baboon females.
#############################################################################################################################

# code pregnancy outcome (live birth or fetal loss) as a binary variable where fetal loss=1 and live birth=0
data$pregnancy_outcome2 <- as.factor(ifelse(data$pregnancy_outcome=="live_birth", 0, 1))

# order female id by her age at death or censorship
data$female_id_ordered <- as.numeric(reorder(data$female_id, data$age_at_death_or_censorship))

# Plot figure 1
ggplot(data=data) + geom_linerange(aes(x=female_id_ordered, ymin=0, ymax=age_at_death_or_censorship, color=genome_wide_anubis_ancestry))  + geom_segment(data=data[data$female_censorship=="dead",], aes(x=female_id_ordered-0.65, y=age_at_death_or_censorship, xend=female_id_ordered+0.65, yend=age_at_death_or_censorship), size=0.25) +  geom_segment(data=data[data$female_censorship=="dead",], aes(x=female_id_ordered, y=age_at_death_or_censorship, xend=female_id_ordered, yend=age_at_death_or_censorship-0.4), size=0.25) + geom_point(aes(x=female_id_ordered, y=female_age_conception, shape=pregnancy_outcome2, fill=pregnancy_outcome2), size=2.5)  +  coord_flip() + scale_fill_manual(breaks = c("0", "1"), values=c("blue", "red"))  + scale_shape_manual(breaks = c("0", "1"), values=c(1,8)) + theme_classic() + theme(text=element_text(size=30), axis.text = element_text(color="black"),  axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") + labs(x="females ordered by lifespan", y="female age (years)")+scale_colour_gradientn(colors=c("#FFED4F", "orange","#009E73"), values=scales::rescale(c(0,0.25,0.45,0.55,1))) + scale_x_continuous(expand = c(0.015, 0.015))

ggsave("fig1.pdf", width=16, height=22)

#############################################################################################################################
# Figure 2. Records of fetal loss by trimester in the data set.
#############################################################################################################################

# Plot figure 2
ggplot(data=data[data$pregnancy_outcome=="fetal_loss",]) + geom_rect(aes(xmin = 0, xmax = 60, ymin = 0, ymax = Inf), fill = "grey90", alpha = 0.1) + geom_rect(aes(xmin = 120, xmax = 201, ymin = 0, ymax = Inf), fill = "grey90", alpha = 0.1) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = Inf), linetype = "dashed", col = "black")+geom_segment(aes(x = 60, y = 0, xend = 60, yend = Inf), linetype = "dashed", col = "black") + geom_segment(aes(x = 120, y = 0, xend = 120, yend = Inf), linetype = "dashed", col = "black") + geom_segment(aes(x = 201, y = 0, xend = 201, yend = Inf), linetype = "dashed", col = "black") + geom_histogram(aes(pregnancy_length), color="black", fill="#99D492", binwidth=3) + scale_x_continuous(name="length of pregnancy that ended in fetal loss (days)") + scale_y_continuous(name="count of fetal losses") + annotate("text",x=30,y=8.5, label="first trimester", size=6) + annotate("text",x=90,y=8.5, label="second trimester", size=6) + annotate("text",x=160,y=8.5, label="third trimester", size=6) + theme_classic()+theme(text=element_text(size=20), axis.text = element_text(color="black"))

ggsave("fig2.pdf", width=8, height=8)

#############################################################################################################################
# Figure 3. Records of fetal loss by trimester in the data set.
#############################################################################################################################

# Need to run main logistic regression model for some panels in figure 3 (code in logistic_regression_models.R but also below)
# Note that we could use I(variable^2) for the quadratic transformation of the variable in the model, however, for figure making, it's easier to create a squared variable and then input that variable into the model directly (get the same results no matter which way you do it)
data$genome_wide_anubis_ancestry_squared <- data$genome_wide_anubis_ancestry^2
data$female_age_conception_squared <- data$female_age_conception^2

m1 <- glmmTMB(pregnancy_outcome2 ~ genome_wide_anubis_ancestry + genome_wide_anubis_ancestry_squared  +  female_age_conception + female_age_conception_squared + previous_fetal_losses + female_ordinal_rank_conception + group_size_females_conception + temp_max_2months_conception + temp_max_2months_endpreg +   rainfall_5months_conception + rainfall_5months_endpreg +  habitat_quality + (1|female_id), family=binomial, data=data)

####################
# Fig. 3A
####################

main_model_results <- as.data.frame(summary(m1)$coefficients$cond) # get effect sizes (i.e., betas), standard errors, z values, and p-values from the model
main_model_results <- setDT(main_model_results, keep.rownames = T) # set predictor variable names/intercept as their own column (not just row names)
names(main_model_results)[1] <- c("predictor")
main_model_results$number <- seq(1, nrow(main_model_results)) # assign numbers to each predictor variable (i.e., row) which will help us not right out each variable name later

# For genome-wide ancestry-related variables and the intercept, divide betas and standard errors by 100 to place them on a similar scale as the other predictors in order to facilitate visualization
main_model_results$Estimate_scaled <- ifelse(main_model_results$predictor %in% c("(Intercept)",                       "genome_wide_anubis_ancestry","genome_wide_anubis_ancestry_squared"), main_model_results$Estimate/100, main_model_results$Estimate)
main_model_results$SE_scaled <- ifelse(main_model_results$predictor %in% c("(Intercept)", "genome_wide_anubis_ancestry","genome_wide_anubis_ancestry_squared"), main_model_results$`Std. Error`/100, main_model_results$`Std. Error`)
# For maternal ancestry squared, multiply beta and standard error by 100 to place them on a similar scale as the other predictors in order to facilitate visualization
main_model_results$Estimate_scaled <- ifelse(main_model_results$predictor %in% c("female_age_conception_squared"), main_model_results$Estimate*100, main_model_results$Estimate_scaled)
main_model_results$SE_scaled <- ifelse(main_model_results$predictor %in% c("female_age_conception_squared"), main_model_results$`Std. Error`*100, main_model_results$SE_scaled)

# Assign predictor variables into broad categories (individual-level effects, social and demographic effects, ecological effects, or intercept)
main_model_results$category <- ifelse(main_model_results$number>1 & main_model_results$number<=6, "individual effects", ifelse(main_model_results$number>6 & main_model_results$number<=8, "social and demographic effects", ifelse(main_model_results$number>8, "ecological effects", "intercept")))

# Assign a color to each predictor variable category
main_model_results$color <- ifelse(main_model_results$number>1 & main_model_results$number<=6, "#1EAE9E", ifelse(main_model_results$number>6 & main_model_results$number<=8, "#8478A5", ifelse(main_model_results$number>8, "#C77657", "grey50")))

# Identify which predictor variables were scaled and which were not
main_model_results$scaled <- ifelse(main_model_results$predictor %in% c("(Intercept)",                       "genome_wide_anubis_ancestry","genome_wide_anubis_ancestry_squared", "female_age_conception_squared"), "yes", "no") 

# Color key
key <- main_model_results$color[order(-main_model_results$number)]

# Plot figure 3A
fig3a <- ggplot(data = main_model_results, aes(x=as.factor(-number), y=Estimate_scaled, ymin=Estimate_scaled-SE_scaled, ymax=Estimate_scaled+SE_scaled, color=category, linetype=scaled)) + geom_hline(yintercept = 0, color="black", linetype="dashed", size=0.75) + geom_pointrange(size=1)  + annotate("text", x = "-5", y = 1.55, label = "*", size=8, vjust=0.75) + annotate("text", x = "-10", y = 0.38, label = "*", size=8, vjust=0.75) + annotate("text", x = "-13", y = 1.42, label = "*", size=8, vjust=0.75) + coord_flip() +  theme_classic() + theme(axis.title.y = element_blank(), text=element_text(size=20), axis.text = element_text(color="black"), axis.text.y = element_text(color=key, size=20.5), legend.position = "none") +  coord_flip() + scale_x_discrete(labels=c("-1"="intercept", "-2"="maternal genetic ancestry", "-3"=bquote("maternal genetic ancestry"^2), "-4"="maternal age", "-5"=bquote("maternal age"^2), "-6"="number of previous fetal losses", "-7"="maternal dominance rank", "-8"="group size", "-9"= "temperature before conception", "-10"= "temperature before end of pregnancy", "-11"= "rainfall before conception", "-12"="rainfall before end of pregnancy", "-13"="habitat quality")) + scale_y_continuous(name="effect size") + scale_color_manual(values = c("intercept"="grey50", "individual effects"="#1EAE9E", "social and demographic effects"="#8478A5", "ecological effects"="#C77657")); fig3a

####################
# Fig. 3B
####################

# For the dashed line which shows the predicted relationship between female age and fetal loss based on model estimates, assuming average values for all other covariates
# Get the fetal loss probability for female age assuming average values for all other covariates (i.e., the female is average with respect to individual, social, demographic, and ecological effects apart from her age which will vary)
# Note that we exclude female age-related variables and female id is coded as NA
tmp <- data.frame(genome_wide_anubis_ancestry=mean(data$genome_wide_anubis_ancestry), genome_wide_anubis_ancestry_squared=mean(data$genome_wide_anubis_ancestry)^2, group_size_females_conception=mean(data$group_size_females_conception), previous_fetal_losses=mean(data$previous_fetal_losses), habitat_quality="high", female_ordinal_rank_conception=mean(data$female_ordinal_rank_conception),temp_max_2months_conception=mean(data$temp_max_2months_conception), temp_max_2months_endpreg=mean(data$temp_max_2months_endpreg), rainfall_5months_endpreg=mean(data$rainfall_5months_endpreg), rainfall_5months_conception=mean(data$rainfall_5months_conception), female_id=NA)

# Get the range of empirical values for female age
age <- as.data.frame(seq(min(data$female_age_conception),max(data$female_age_conception), by=0.5))
colnames(age) <- c("female_age_conception")
age$female_age_conception_squared <- age$female_age_conception^2

tmp2 <- merge(tmp,age) # combine data frame of averaged variables with the data frame of the empirical range of female age

# Get the fetal loss probability (and standard error) for the range of female age values assuming average values for all other covariates
tmp2[c(14:15)] <- predict(m1, tmp2, type="response", re.form=NA, se.fit = TRUE)

# Plot figure 3b
fig3b <- ggplot() + geom_jitter(data=data, aes(female_age_conception, ifelse(pregnancy_outcome2==0, 0, 1), color = ifelse(pregnancy_outcome2==0, "0", "1")), alpha=0.4, height=0.05, size=3) + geom_pointrange(data=tmp2, aes(x=female_age_conception, y=fit, ymin=fit-se.fit, ymax=fit+se.fit), size=0.25)  + scale_x_continuous(name="maternal age at conception (years)") +  theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"), legend.position="none") + scale_y_continuous(name="fetal loss probability") + scale_color_manual(values = c("0" = "skyblue3", "1" = "#99D492")); fig3b

####################
# Fig. 3C
####################

# For the dashed line which shows the predicted relationship between temperature two months before the end of the pregnancy and fetal loss based on model estimates, assuming average values for all other covariates
# Get the fetal loss probability for temperature two months before the end of the pregnancy assuming average values for all other covariates (i.e., the female is average with respect to individual, social, demographic, and ecological effects apart from temperature two months before the end of the pregnancy which will vary)
# Note that we exclude temperature two months before the end of the pregnancy and female id is coded as NA
tmp <- data.frame(genome_wide_anubis_ancestry=mean(data$genome_wide_anubis_ancestry), genome_wide_anubis_ancestry_squared=mean(data$genome_wide_anubis_ancestry)^2, group_size_females_conception=mean(data$group_size_females_conception), habitat_quality="high", previous_fetal_losses=mean(data$previous_fetal_losses), female_ordinal_rank_conception=mean(data$female_ordinal_rank_conception),female_age_conception=mean(data$female_age_conception), female_age_conception_squared=mean(data$female_age_conception)^2,  temp_max_2months_conception=mean(data$temp_max_2months_conception), rainfall_5months_endpreg=mean(data$rainfall_5months_endpreg), rainfall_5months_conception=mean(data$rainfall_5months_conception), female_id=NA)

# Get the empirical range of temperature two months before the end of the pregnancy
temp_endpreg <- as.data.frame(seq(min(data$temp_max_2months_endpreg),max(data$temp_max_2months_endpreg), by=0.25))
colnames(temp_endpreg) <- c("temp_max_2months_endpreg")

tmp2 <- merge(tmp,temp_endpreg) # combine data frame of averaged variables with the data frame of the empirical range of temperature two months before the end of the pregnancy

# Get the fetal loss probability (and standard error) for the range of temperature two months before the end of the pregnancy values assuming average values for all other covariates
tmp2[c(14:15)]  <- predict(m1, tmp2, type="response", re.form=NA, se.fit = TRUE)

# Plot figure 3c
fig3c <- ggplot() + geom_jitter(data=data, aes(temp_max_2months_endpreg, ifelse(pregnancy_outcome2==0, 0, 1), color = ifelse(pregnancy_outcome2==0, "0", "1")), alpha=0.4, height=0.05, size=3) + geom_pointrange(data=tmp2, aes(x=temp_max_2months_endpreg, y=fit, ymin=fit-se.fit, ymax=fit+se.fit), size=0.25)  + scale_x_continuous(name="mean average daily maximum temperature\n2 months before live birth/fetal loss (°C)") +  theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black"), legend.position = "none") + scale_y_continuous(name="fetal loss probability")  + scale_color_manual(values = c("0" = "skyblue3", "1" = "#99D492")); fig3c

####################
# Fig. 3D
####################

# Plot figure 3d
fig3d <- ggplot(data=data)+geom_bar(aes(habitat_quality, fill=relevel(pregnancy_outcome, ref="fetal_loss")), color="black", position="fill", size=1) + theme_classic()+theme(text=element_text(size=20), axis.text = element_text(color="black"), legend.position = "none") + scale_x_discrete( name="habitat quality")+scale_fill_manual(breaks = c("live_birth", "fetal_loss"), values=c("skyblue3", "#99D492"))+scale_y_continuous(name="proportion of all pregnancies") + geom_text(aes(habitat_quality, group=pregnancy_outcome, label=..count..), vjust=1.5, stat="count", position="fill", size=6, color="black"); fig3d

####################
# Fig. 3A-D
####################

# Combine main panel 3A with panels 3B, 3C, 3D
(fig3a | (fig3b + theme(plot.margin = unit(c(0,40,0,0), "pt"))) / (fig3d + theme(plot.margin = unit(c(0,40,0,0), "pt"))) / (fig3b + theme(plot.margin = unit(c(0,40,0,0), "pt"))))  + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1.2, 1))

ggsave("fig3A-D.pdf", height=12, width=16)

#############################################################################################################################
# Figure S1. The total number of pregnancies per female and its relationship to the number of different pregnancy outcomes.
#############################################################################################################################

# What is the relationship between the total number of pregnancies per females and her total number of live births and fetal losses? 

# Get the number of times each female appears in our data set (i.e., how many pregnancies she had)
females <- select(data, female_id) %>% group_by(female_id) %>% tally()
names(females)[2] <- "total_pregs_count"

# Get the total number of pregnancy outcome types (i.e., live birth or fetal loss) per female
tmp <- select(data, female_id, pregnancy_outcome2) %>% group_by(female_id, pregnancy_outcome2) %>% tally()

# Grab just one pregnancy outcome type to combine with the females data frame (can then back-calculate the other pregnancy outcome type based on the total pregnancies)
tmp2 <- tmp[tmp$pregnancy_outcome2==0,]
names(tmp2)[3] <- c("num_live_births")
females <- merge(females, tmp2[c(1,3)], by=c("female_id"))

# Calculate the total number of fetal losses per female
females$num_fetal_losses <- females$total_pregs_count-females$num_live_births

round(summary(lm(num_live_births~total_pregs_count, data=females))$coefficients[2],3) # Total pregs beta

round(summary(lm(num_fetal_losses~total_pregs_count, data=females))$coefficients[2],3) # Total pregs beta

# Plot total pregnancies per female vs. her total number of live births (S1A) and total number of fetal losses (S1B)
# Because multiple females have the same combinations of x-y values (total number of pregnancies-total number of live births or total number of pregnancies-total number of fetal losses), depict multiple females at a given x-y value combination based on the size of the dot
s1a <- ggplot() + geom_count(data=females, aes(total_pregs_count, num_live_births), color="skyblue3") + scale_x_continuous(name="number of pregnancies\nper female", breaks=seq(0,15,by=5)) + scale_y_continuous(name="number of live births per female") + theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black")); s1a
s1b <- ggplot() + geom_count(data=females, aes(total_pregs_count, num_fetal_losses), color="#99D492") + scale_x_continuous(name="number of pregnancies\nper female", breaks=seq(0,15,by=5)) + scale_y_continuous(name="number of fetal losses per female") + theme_classic() + theme(text=element_text(size=20), axis.text = element_text(color="black")); s1b

# Plot S1 panels in a single plot
(s1a | s1b) + plot_annotation(tag_levels = 'A')

ggsave(file="figs1.pdf", width=15, height=7)
