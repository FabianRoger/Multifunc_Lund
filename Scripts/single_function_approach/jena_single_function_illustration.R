
# Project: Review of multifunctionality in ecology, conservation and ecosystem service science

# Title: Illustrate the single function approach using the Jena data

# load relevant libraries
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(DHARMa)
library(corrr)
library(tidyr)
library(forcats)
library(viridis)
library(piecewiseSEM)
library(patchwork)

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
# plot a bivariate plot for each

# write a function to do this

# write a function to make a decent scatterplot with a fitted model

biv_plotter <- function(data.in, model.input, 
                        response, predictor.plot, predictor.mod, ylab, xlab,
                        r2.method = "nagelkerke", mean.pred, sd.pred,
                        x.breaks = c(1, 2, 4, 8, 12, 16) ) {
  
  if( class(model.input) == "lm") {
    
    x <- summary(model.input)
    r <- paste("r2 = ", round(x$r.squared, 1) , sep = "")
    
    f <- x$fstatistic
    fit.test <- paste("F(", f[2], ",", f[3], ")", " = ", round(f[1], 1), sep = "")
    
    p <- pf(f[1], df1 = f[2], df2 = f[3], lower.tail = F)
    
    if (p < 0.001) {
      
      p.val <- ("P < 0.001")
      
    } else {
      
      p.val <- paste("P = ", round(p, 3))
      
    }
    
  } else if ( "glm" %in% class(model.input) ) {
    
    r <- piecewiseSEM::rsquared(model.input, method = r2.method) 
    r <- paste("r2 = ", round(r$R.squared, 1), sep = "")
    
    y <- drop1(model.input, test = c("Chisq"))
    
    fit.test <- paste("chisq = ", round(y$LRT[2], 1))
    
    p <- y$`Pr(>Chi)`[2]
    
    if (p < 0.001) {
      
      p.val <- ("P < 0.001")
      
    } else {
      
      p.val <- paste("P = ", round(y$`Pr(>Chi)`[2], 3))
      
    }
    
  } else {
    
    stop("error, function only works for lm and glm objects")
    
  }
  
  x.range <- range(data.in[[predictor.mod]])
  new.dat <- data.frame(seq(x.range[1], x.range[2], 0.01))
  names(new.dat) <- predictor.mod
  
  pred.mod <- predict(model.input, newdata = new.dat, se = TRUE, type = "response")
  
  pred.dat <- 
    data.frame((new.dat[[predictor.mod]]*sd.pred) + mean.pred,
               pred.mod$fit, 
               pred.mod$fit - pred.mod$se.fit,
               pred.mod$fit + pred.mod$se.fit )
  
  names(pred.dat) <- c(predictor.plot,
                       paste(response, "pred", sep = "."),
                       paste(response, "se.low", sep = "."),
                       paste(response, "se.upp", sep = ".") )
  
  pred.names <- names(pred.dat)
  
  # non-gaussian glm's
  
  ggplot() +
    geom_jitter(data = data.in,
                mapping = aes_string(x = predictor.plot, y = response ), width = 0.1) +
    geom_line(data = pred.dat,
              mapping = aes_string(pred.names[1], pred.names[2] )) +
    geom_ribbon(data = pred.dat,
                mapping = aes_string(x = pred.names[1], ymax = pred.names[4],
                                     ymin = pred.names[3] ), fill = "grey", alpha = 0.5 ) +
    scale_x_continuous(breaks = x.breaks ) +
    annotate(geom = 'text', label = r, x = Inf, y = -Inf, 
             hjust = 1.25, vjust = -5, size = 3.5) +
    annotate(geom = 'text', label = fit.test, x = Inf, y = -Inf, 
             hjust = 1.1, vjust = -3.3, size = 3.5) +
    annotate(geom = 'text', label = p.val, x = Inf, y = -Inf, 
             hjust = 1.15, vjust = -1.6, size = 3.5) +
    ylab(ylab) +
    xlab(xlab) +
    theme_meta()
  
}


# biomass
hist(sqrt(jena.dat$biomass) )
lm.bm <- lm((biomass) ~ sowndiv_scale, data = jena.dat)
# plot(lm.bm) # assumptions look good

summary(lm.bm)

# soilC
hist(jena.dat$soilorgC)
lm.sc <- lm(soilorgC ~ sowndiv_scale, data = jena.dat)
# plot(lm.sc) # assumptions look good

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
lm.po2 <- MASS::glm.nb(round(poll, 0) ~ sowndiv_scale, data = jena.dat, link = log)

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
drop1(lm.po2, test = "Chisq")

# rootBN
hist(jena.dat$rootBM)

jena.dat <- 
  jena.dat %>%
  mutate(sqrt_rootBM = sqrt(rootBM))

lm.rb <- lm( sqrt_rootBM ~ sowndiv_scale, data = jena.dat)
# plot(lm.rb) # assumptions look good

summary(lm.rb)

# plantCN
hist(jena.dat$plantCN)
lm.cn <- lm(plantCN ~ sowndiv_scale, data = jena.dat)
# plot(lm.cn) # assumptions are okay (normality a bit off)

summary(lm.cn)

# generate bivariate plots for these five models and functions
mod.list <- list(lm.bm, lm.sc, lm.po2, lm.rb, lm.cn)
resp.list <- c("biomass", "soilorgC", "poll", "sqrt_rootBM", "plantCN")
ylab.list <- c("Biomass, g m-2", "Soil carbon, %", "Pollinator abun.",
               "sqrt(Root biomass, g m-2 )", "Plant C/N")

mod.plots <- vector("list", length = length(mod.list))
for (i in 1:length(mod.list)) {
  
  mod.plots[[i]] <- 
    biv_plotter(data.in = jena.dat, model.input = mod.list[[i]], 
                response = resp.list[[i]], predictor.plot = "sowndiv", 
                predictor.mod = "sowndiv_scale", ylab= ylab.list[[i]], xlab = "Sown species richness",
                r2.method = "nagelkerke", mean.pred = mean(jena.dat$sowndiv), sd.pred = sd(jena.dat$sowndiv),
                x.breaks = c(1, 2, 4, 8, 12, 16))
  
}
  
mod.plots[[5]]

# generate a correlation plot between these functions
jena.corr.raw <- jena.dat[, func.names]
names(jena.corr.raw) <- c("Biomass", "Soil Org. C", "Pollination", "Root biomass", "C/N")

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
# # https://www.cedricscherer.com/2019/08/05/a-ggplot2-tutorial-for-beautiful-plotting-in-r/#charts
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

jena.corr$p.val <- p.out$corr

jena.corr.plot <- 
  ggplot(jena.corr, aes(rowname, fct_rev(colname),
                 fill = corr)) +
  geom_tile() +
  geom_text(aes(
    label = format(round(corr, 2), nsmall = 2),
    color = abs(corr) == 1
  ), size = 4) +
  geom_text(aes(
    label = (p.val) 
    ), nudge_y = 0.1,  nudge_x = 0.3, size = 6) +
  coord_fixed(expand = FALSE) +
  scale_color_manual(values = c("black", "white"),
                     guide = "none") +
  scale_fill_viridis(
    option = "D",
    na.value = "white",
    direction = -1, alpha = 0.8, begin = 0.1,
    name = "rho") +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(.85, .8),
        axis.text = element_text(colour = "black"))

# merge these plots in a logical way that doesn't get too ridiculous space-wise
mod.plots[[1]] + mod.plots[[2]] + mod.plots[[3]] +
  mod.plots[[4]] + mod.plots[[5]] + jena.corr.plot +
  plot_layout(ncol = 2, nrow = 3, widths = c(1, 1))

### END
