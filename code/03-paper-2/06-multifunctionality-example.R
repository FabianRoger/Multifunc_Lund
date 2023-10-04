#'
#' @title Example of the workflow
#' 
#' @description Reanalyse data from Gamfeldt et al. (2013) to show how we might
#' operationalise some of the steps recommended
#' 

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# load plotting theme
source("code/helper-plotting-theme.R")

# load the data
mf_df <- readr::read_delim("data/gamfeldt_2013_forest_data.txt", delim = "\t")
head(mf_df)

# check the data
View(mf_df)

# subset out the relevant columns
mf_df <-
  mf_df |>
  dplyr::select(year, trakt, plot, S, production, Cstock, bilberries, dw_occ)

# check the initial relationships among the data
pairs(mf_df[,-c(1, 2, 3) ])

# get only complete cases
mf_df <- mf_df[complete.cases(mf_df),]

# randomly select a single plot from each trakt
unique(mf_df$year)

# how many unique trakts are there in each year
mf_df |>
  dplyr::group_by(year) |>
  dplyr::summarise(n_trakts = length(unique(trakt)))

# take only data from 1999 and get one plot per trakt
set.seed(854875)
mf_df <- 
  mf_df |>
  dplyr::filter(year == 1999) |>
  dplyr::group_by(trakt) |>
  dplyr::sample_n(size = 1)

# convert Cstock to 100 g per m2
mf_df <-
  mf_df |>
  dplyr::mutate(Cstock = Cstock/100)

# logit transform bilberries
mf_df <-
  mf_df |>
  dplyr::mutate(bilberries = log( (bilberries+0.001)/(1 - (bilberries+0.001) ) ))

# remove one datapoint with extremely high Cstock
mf_df <- 
  mf_df |>
  dplyr::filter(Cstock < 400)

# check the raw relationships again
pairs(mf_df[,-c(1, 2, 3) ])

# draw the DAG
dag1 <- dagitty::downloadGraph("dagitty.net/mpkrtMC")
plot(dag1)

# test the conditional independencies
dagitty::impliedConditionalIndependencies(dag1)

# test the conditional independencies using likelihood ratio tests

# Blbr _||_ TrSR | TrPr
lm1 <- lm(bilberries ~ production, data = mf_df)
lm2 <- lm(bilberries ~ production + S, data = mf_df)

# get the likelihood ratio test
ci1 <- lmtest::lrtest(lm2, lm1)
p1 <- ci1$`Pr(>Chisq)`[2]

# DW _||_ SolC | TrPr
lm1 <- glm(dw_occ ~ production, family = binomial(link = "logit"), data = mf_df)
lm2 <- glm(dw_occ ~ production + Cstock, family = binomial(link = "logit"), data = mf_df)

# get the likelihood ratio test
ci2 <- lmtest::lrtest(lm2, lm1)
p2 <- ci2$`Pr(>Chisq)`[2]

# DW _||_ TrSR | TrPr
lm1 <- glm(dw_occ ~ production, family = binomial(link = "logit"), data = mf_df)
lm2 <- glm(dw_occ ~ production + S, family = binomial(link = "logit"), data = mf_df)

# get the likelihood ratio test
ci3 <- lmtest::lrtest(lm2, lm1)
p3 <- ci3$`Pr(>Chisq)`[2]

# SolC _||_ TrSR | TrPr
lm1 <- lm(Cstock ~ production, data = mf_df)
lm2 <- lm(Cstock ~ production + S, data = mf_df)

# get the likelihood ratio test
ci4 <- lmtest::lrtest(lm2, lm1)
p4 <- ci4$`Pr(>Chisq)`[2]

# conduct p-value correction
p_vals <- p.adjust(p = c(p1, p2, p3, p4), method = "bonferroni")
print(p_vals)

# get the coefficients for each path
plot(dag1)

# adjustment set required for the effect of TreeSR on TreeProd
dagitty::adjustmentSets(x = dag1, exposure = "TreeSR", outcome = "TreeProd", effect = "direct")

# fit the model
lm1 <- lm(production ~ S, data = mf_df)
summary(lm1)

ggplot(data = mf_df, mapping = aes(x = S, y = production)) +
  geom_jitter(width = 0.1) +
  ylab("Production") +
  xlab("Tree species richness") +
  theme_meta()

# adjustment set required for the effect of TreeProd on SoilC
dagitty::adjustmentSets(x = dag1, exposure = "TreeProd", outcome = "SoilC", effect = "direct")

# fit the model
lm2 <- lm(Cstock ~ production, data = mf_df)
summary(lm2)

ggplot(data = mf_df, mapping = aes(x = production, y = Cstock)) +
  geom_jitter(width = 0.1) +
  ylab("Cstock") +
  xlab("Production") +
  theme_meta()

# adjustment set required for the effect of TreeProd on SoilC
dagitty::adjustmentSets(x = dag1, exposure = "SoilC", outcome = "Bilberry", effect = "direct")

# fit the model
lm3 <- lm(bilberries ~ production + Cstock, data = mf_df)
summary(lm3)

# fit models to get the residual plots
lm3a <- lm(bilberries ~ production, data = mf_df)
summary(lm3a)
lm3b <- lm(Cstock ~ production, data = mf_df)
summary(lm3b)

# make a residual plot
ggplot(data = mf_df, mapping = aes(x = residuals(lm3b), y = residuals(lm3a))) +
  geom_jitter(width = 0.1) +
  ylab("Bilberry production") +
  xlab("Soil carbon") +
  theme_meta()

# adjustment set required for the effect of TreeProd on Bilberry
dagitty::adjustmentSets(x = dag1, exposure = "TreeProd", outcome = "Bilberry", effect = "direct")

# fit the model
lm4 <- lm(bilberries ~ production + dw_occ + Cstock, data = mf_df)
summary(lm4)

# fit models to get the residual plots
lm4a <- lm(bilberries ~ dw_occ + Cstock, data = mf_df)
summary(lm4a)
lm4b <- lm(production ~ dw_occ + Cstock, data = mf_df)
summary(lm4b)

ggplot(data = mf_df, mapping = aes(x = residuals(lm4b), y = residuals(lm4a))) +
  geom_jitter(width = 0.1) +
  ylab("Bilberry production") +
  xlab("Production") +
  theme_meta()

# adjustment set required for the effect of TreeProd on DW
dagitty::adjustmentSets(x = dag1, exposure = "TreeProd", outcome = "DW", effect = "direct")

# fit the model
lm5 <- glm(dw_occ ~ production, family = binomial(link = "logit"), data = mf_df)
summary(lm5)

ggplot(data = mf_df, mapping = aes(x = production, y = dw_occ)) +
  geom_jitter(width = 0.1) +
  ylab("Dead wood") +
  xlab("Production") +
  theme_meta()

# adjustment set required for the effect of DW on Bilberries
dagitty::adjustmentSets(x = dag1, exposure = "DW", outcome = "Bilberry", effect = "direct")

# fit the model
lm6 <- lm(bilberries ~ production + dw_occ, data = mf_df)
summary(lm6)

# fit models to get the residual plots
lm6a <- lm(bilberries ~ production, data = mf_df)
summary(lm4a)
lm6b <- glm(dw_occ ~ production, family = binomial(link = "logit"), data = mf_df)
summary(lm6b)

# inverse logit transform
x <- exp(residuals(lm6b))/(1+exp(residuals(lm6b)))

ggplot(data = mf_df, mapping = aes(x = x, y = residuals(lm6a))) +
  geom_jitter(width = 0.1) +
  ylab("Bilberry production") +
  xlab("Dead wood") +
  theme_meta()




