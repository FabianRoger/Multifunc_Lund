
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Illustrate the single function approach using the Jena data

# load relevant libraries
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(DHARMa)
library(MASS)
library(corrr)
library(tidyr)
library(forcats)

rm(list = ls())

# link to scripts with the relevant functions
source(here("Scripts/function_plotting_theme.R"))

# read in the Jena data
jena.dat <- read_csv(file = here("data/jena_data_Jochum_2020_clean.csv"))

# check the column names
col_names <- names(jena.dat)

# remove the individual species abundances
spp <- col_names[grepl(pattern = "[A-Z][a-z]{2}[.][a-z]{3}", col_names)]
length(spp)
rm(col_names)

jena.dat <- 
  jena.dat %>%
  select(-all_of(spp))

# choose four ecosystem functions to examine in the single function approach
names(jena.dat)
func.names <- c("biomass", "soilorgC", "poll", "rootBM", "plantCN")

jena.dat <- 
  jena.dat %>%
  select(plotcode, block, sowndiv, realised_richness1, all_of(func.names) )

# check the summary statistics and distributions
summary(jena.dat)
length(unique(jena.dat$block))
lapply(jena.dat[, func.names], hist)

# scale the sowndiv as this is our predictor
jena.dat <- 
  jena.dat %>%
  mutate(sowndiv_scale = scale(sowndiv, center = TRUE, scale = TRUE)[,1])

# fit models with each variable as a response and sowndiv as a predictor

# biomass
hist(sqrt(jena.dat$biomass) )
lm.bm <- glm((biomass) ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.bm) # assumptions look good

summary(lm.bm)

# soilC
hist(jena.dat$soilorgC)
lm.sc <- glm(soilorgC ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.sc)

summary(lm.sc)

# pollination
hist(jena.dat$poll)
lm.po <- glm(round(poll, 0) ~ sowndiv_scale, data = jena.dat, family = "poisson")
plot(lm.po)

lm.po.s <- summary(lm.po)
lm.po.s

# ratio should be approximately one: it is seven...
lm.po.s$deviance/lm.po.s$df.residual
lm.po.s$dispersion

# use the DHARMa library
sim_fmp <- simulateResiduals(lm.po, refit=F)
testOverdispersion(sim_fmp)

# what to do about overdispersion? 
# 1. we can fit a dispersion parameter
# 2. we can use the negative binomial which can model different variances
# 3. we can use an observation-level random effect

# https://biometry.github.io/APES/LectureNotes/2016-JAGS/Overdispersion/OverdispersionJAGS.pdf

# 2. use a negative binomial distribution
library(MASS)
lm.po2 <- glm.nb(round(poll, 0) ~ sowndiv_scale, data = jena.dat, link = log)

# check for overdispersion
lm.po2$deviance/lm.po2$df.residual

# graphical analysis of the residuals
plot(lm.po2)

# overdispersion test using the DHARMa library
sim_fmp <- simulateResiduals(lm.po2, refit=F)
testOverdispersion(sim_fmp)
plot(sim_fmp)

# run the summary function
summary(lm.po2)

# rootBN
hist(jena.dat$rootBM)
lm.rb <- glm(sqrt(rootBM) ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.rb)

summary(lm.rb)

# plantCN
hist(jena.dat$plantCN)
lm.cn <- glm(plantCN ~ sowndiv_scale, data = jena.dat, family = "gaussian")
plot(lm.cn)

summary(lm.cn)

# generate a correlation plot between these functions
jena.corr.raw <- jena.dat[, func.names]
names(jena.corr) <- c("Biomass", "Soil Org. C", "Pollination", "Root biomass", "C/N")

jena.corr <-
  jena.corr.raw %>%
  corrr::correlate(diagonal = 1) %>%
  corrr::shave(upper = FALSE)

# calculate the correlation and p-values using Hmisc
p.corr <- Hmisc::rcorr(as.matrix(jena.corr.raw))

# extract a matrix of p-values
p.out <- p.corr$P

# set the diagonals to 1
diag(p.out) <- 1

#replace the lower half with NAs to match the correlation
p.out[lower.tri(p.out, diag = FALSE)] <- NA

# correct the p-values using hocberg's correction for multiple comparisons
p.out[ (p.out == p.out | is.na(p.out)) ] <- ifelse(p.adjust(p = p.out, method = "hochberg") < 0.05, "*", "") 

# convert to a tibble
p.out <- as_tibble(p.out)

# prepare the correlation data for plotting
jena.corr <- 
  jena.corr %>%
  pivot_longer(
    cols = -term,
    names_to = "colname",
    values_to = "corr"
  ) %>%
  mutate(rowname = fct_inorder(term),
         colname = fct_inorder(colname))

# prepare the p-value data for plotting
p.out$term <- row.names(p.corr$P) 

p.out <- 
  p.out %>%
  pivot_longer(
    cols = -term,
    names_to = "colname",
    values_to = "corr"
  ) %>%
  mutate(rowname = fct_inorder(term),
         colname = fct_inorder(colname))

ggplot(jena.corr, aes(rowname, fct_rev(colname),
                 fill = corr)) +
  geom_tile() +
  geom_text(aes(
    label = format(round(corr, 2), nsmall = 2),
    color = abs(corr) < .75
  )) +
  geom_text(data = p.out, 
            mapping = aes(label = corr)) +
  coord_fixed(expand = FALSE) +
  scale_color_manual(values = c("white", "black"),
                     guide = "none") +
  scale_fill_distiller(
    palette = "PuOr", na.value = "white",
    direction = 1, limits = c(-1, 1)
  ) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(.85, .8))

