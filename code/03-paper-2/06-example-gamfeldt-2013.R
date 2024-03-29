#'
#' @title Example of the workflow
#' 
#' @description Reanalyse data from Gamfeldt et al. (2013) to show how we might
#' operationalise some of the steps recommended
#' 

# load the data
mf_df <- readr::read_delim("data/gamfeldt_2013_forest_data.txt", delim = "\t")
head(mf_df)

# check the data
# View(mf_df)

# subset out the relevant columns
mf_df <-
  mf_df |>
  dplyr::select(year, trakt, plot, S, production, Cstock, bilberries)

# check the initial relationships among the data
# pairs(mf_df[,-c(1, 2, 3) ])

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
  dplyr::sample_n(size = 1) |>
  dplyr::ungroup()

# remove plots with zero tree species
mf_df <-
  mf_df |>
  dplyr::filter(S > 0)
  
# remove one datapoint with extremely high Cstock
mf_df <- 
  mf_df |>
  dplyr::filter(Cstock < 40000)

# convert Cstock to 100 g per m2 and log transform
mf_df <-
  mf_df |>
  dplyr::mutate(Cstock = Cstock/1000) |>
  dplyr::mutate(Cstock = log(Cstock))

# logit transform bilberries
mf_df <-
  mf_df |>
  dplyr::mutate(bilberries = log( (bilberries+0.001)/(1 - (bilberries+0.001) ) ))

# sqrt transform production
mf_df <- 
  mf_df |>
  dplyr::mutate(production = sqrt(production))

# check the raw relationships again
# pairs(mf_df[,-c(1, 2, 3) ])

# how many samples do we have?
nrow(mf_df)

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
print(ci1)

# try different conditional independence testing
bnlearn::ci.test("S", "bilberries", "production", data = mf_df, test = "cor")
with(mf_df, {
  GeneralisedCovarianceMeasure::gcm.test(X = S, Y = bilberries, Z = production, alpha = 0.05, regr.method = "xgboost")
})

# SolC _||_ TrSR | TrPr
lm1 <- lm(Cstock ~ production, data = mf_df)
lm2 <- lm(Cstock ~ production + S, data = mf_df)

# get the likelihood ratio test
ci2 <- lmtest::lrtest(lm2, lm1)
p2 <- ci2$`Pr(>Chisq)`[2]
print(ci2)

# try different conditional independence testing
bnlearn::ci.test("S", "Cstock", "production", data = mf_df, test = "cor")
with(mf_df, {
  GeneralisedCovarianceMeasure::gcm.test(X = S, Y = Cstock, Z = production, alpha = 0.05, regr.method = "xgboost")
  })
  
# conduct p-value correction
p_vals <- p.adjust(p = c(p1, p2), method = "holm")
print(p_vals)

# set-up plotting colours
cols <- wesanderson::wes_palette(name = "Darjeeling2", n = 5)
print(cols)

# select the relevant colours
cols <- cols[2:5]

# set-up axis labels
S <- "Tree species richness"
SC <- expression("ln[ Soil carbon"~(kg~m^{2})~"]")
SC_res <- expression("ln[ Soil carbon"~(kg~m^{2})~"] *")
B <- "logit[ Bilberry cover ]"
B_res <- "logit[ Bilberry cover ] *"
P <- expression(sqrt("Tree production"~(kg~m^{2}~year^{-1})))
P_res <- expression(sqrt("Tree production"~(kg~m^{2}~year^{-1}))~"*")

# get the coefficients for each path
plot(dag1)

# adjustment set required for the effect of TreeSR on TreeProd
dagitty::adjustmentSets(x = dag1, exposure = "TreeSR", outcome = "TreeProd", effect = "direct")

# fit the model
lm1 <- lm(production ~ S, data = mf_df)
summary(lm1)

# get predictive distribution
pred1 <- data.frame(S = seq(1, 6, 0.1))
pred1 <- dplyr::bind_cols(pred1, predict(lm1, pred1, interval = "confidence"))

q1 <- 
  ggplot() +
  geom_jitter(data = mf_df, mapping = aes(x = S, y = production), width = 0.1,
              colour = cols[1]) +
  geom_line(data = pred1, mapping = aes(x = S, y = fit),
            colour = cols[1]) +
  geom_ribbon(data = pred1, mapping = aes(x = S, ymin = lwr, ymax = upr), 
              alpha = 0.1, fill = cols[1]) +
  ylab(P) +
  xlab(S) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.4,
           label = lm_eqn(lm1), parse = TRUE) +
  theme_meta()

# adjustment set required for the effect of TreeProd on SoilC
dagitty::adjustmentSets(x = dag1, exposure = "TreeProd", outcome = "SoilC", effect = "direct")

# fit the model
lm2 <- lm(Cstock ~ production, data = mf_df)
summary(lm2)

# get predictive distribution
pred2 <- data.frame(production = seq(0.01, 1.3, 0.02))
pred2 <- dplyr::bind_cols(pred2, predict(lm2, pred2, interval = "confidence"))

q2 <- 
  ggplot() +
  geom_point(data = mf_df, mapping = aes(x = production, y = Cstock),
             colour = cols[2]) +
  geom_line(data = pred2, mapping = aes(x = production, y = fit),
            colour = cols[2]) +
  geom_ribbon(data = pred2, mapping = aes(x = production, ymin = lwr, ymax = upr), 
              alpha = 0.1, fill = cols[2]) +
  ylab(SC) +
  xlab(P) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.4,
           label = lm_eqn(lm2), parse = TRUE) +
  theme_meta()

# adjustment set required for the effect of TreeProd on Bilberries
dagitty::adjustmentSets(x = dag1, exposure = "TreeProd", outcome = "Bilberries", effect = "direct")

# fit the model
lm3 <- lm(bilberries ~ production + Cstock, data = mf_df)
summary(lm3)

# fit models to get the residual plots
lm3a <- lm(bilberries ~ Cstock, data = mf_df)
summary(lm3a)
lm3b <- lm(production ~ Cstock, data = mf_df)
summary(lm3b)

# fit a model of the residuals on production on Cstock
lm3c_df <- data.frame(production_res = residuals(lm3b),
                      bilberries_res = residuals(lm3a))
lm3c <- lm(bilberries_res ~ production_res, data = lm3c_df)
summary(lm3c)

# get predictive distribution
pred3 <- data.frame(production_res = seq(-0.6, 0.71, 0.02))
pred3 <- dplyr::bind_cols(pred3, predict(lm3c, pred3, interval = "confidence"))

q3 <-
  ggplot() +
  geom_point(data = lm3c_df, mapping = aes(x = production_res, y = bilberries_res),
             colour = cols[3]) +
  geom_line(data = pred3, mapping = aes(x = production_res, y = fit),
            colour = cols[3]) +
  geom_ribbon(data = pred3, mapping = aes(x = production_res, ymin = lwr, ymax = upr), 
              alpha = 0.2, 
              fill = cols[3]) +
  ylab(B_res) +
  xlab(P_res) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.4,
           label = lm_eqn(lm3c), parse = TRUE) +
  theme_meta()

# adjustment set required for the effect of SoilC on Bilberries
dagitty::adjustmentSets(x = dag1, exposure = "SoilC", outcome = "Bilberries", effect = "direct")

# fit the model
lm4 <- lm(bilberries ~ production + Cstock, data = mf_df)
summary(lm4)

# fit models to get the residual plots
lm4a <- lm(bilberries ~ production, data = mf_df)
summary(lm4a)
lm4b <- lm(Cstock ~ production, data = mf_df)
summary(lm4b)

# fit a model of the residuals on production on Cstock
lm4c_df <- data.frame(Cstock_res = residuals(lm4b),
                      bilberries_res = residuals(lm4a))
lm4c <- lm(bilberries_res ~ Cstock_res, data = lm4c_df)
summary(lm4c)

# get predictive distribution
pred4 <- data.frame(Cstock_res = seq(-1.81, 2.1, 0.02))
pred4 <- dplyr::bind_cols(pred4, predict(lm4c, pred4, interval = "confidence"))

q4 <- 
  ggplot() +
  geom_point(data = lm4c_df, mapping = aes(x = Cstock_res, y = bilberries_res),
             colour = cols[4]) +
  geom_line(data = pred4, mapping = aes(x = Cstock_res, y = fit),
            colour = cols[4]) +
  geom_ribbon(data = pred4, mapping = aes(x = Cstock_res, ymin = lwr, ymax = upr), 
              alpha = 0.1,
              fill = cols[4]) +
  ylab(B_res) +
  xlab(SC_res) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.4,
           label = lm_eqn(lm4c), parse = TRUE) +
  theme_meta()

# output these plots
q_list <- list(q1, q2, q3, q4)

# combine into a single plot
q14 <- 
  cowplot::plot_grid(plotlist = q_list, align = "hv", ncol = 2, nrow = 2)
plot(q14)

ggsave("figures-paper-2/fig_4.pdf", q14,
       units = "cm", width = 14, height = 14)

# test the effect of tree species richness on average multifunctionality

# standardise the functions
funcs <- dplyr::select(mf_df, production, Cstock, bilberries)
funcs <- apply(funcs, 2, function(x) (x-min(x))/(max(x)-min(x)) )

# calculate average multifunctionality
mf_df$mf_ave <- apply(funcs, 1, mean)

# calculate ENFQ1
funcs <- as.data.frame(funcs)
mf_df$mf_ENFQ1 <- multifunc::getMF_eff(data = funcs, vars = names(funcs), q = 1, standardized = FALSE)
  
# calculate ENFQ2
mf_df$mf_ENFQ2 <- multifunc::getMF_eff(data = funcs, vars = names(funcs), q = 2, standardized = FALSE)

# is there an effect of biodiversity on ENFQ1 - multifunctionality?
lm5 <- lm(mf_ENFQ1 ~ S, data = mf_df)
summary(lm5)

# get predictive distribution
pred5 <- data.frame(S = seq(1, 6, 0.1))
pred5 <- dplyr::bind_cols(pred5, predict(lm5, pred5, interval = "confidence"))

# plot the results
r1 <- 
  ggplot() +
  geom_jitter(data = mf_df, mapping = aes(x = S, y = mf_ENFQ1), width = 0.1) +
  geom_line(data = pred5, mapping = aes(x = S, y = fit)) +
  geom_ribbon(data = pred5, mapping = aes(x = S, ymin = lwr, ymax = upr), 
              alpha = 0.1) +
  ylab("ENF-Q1 multifunctionality") +
  xlab(S) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.3, vjust = -1.4,
           label = lm_eqn(lm5), parse = TRUE) +
  theme_meta()
plot(r1)

# is there an effect of biodiversity on ENF-Q2 multifunctionality?
lm6 <- lm(mf_ENFQ2 ~ S, data = mf_df)
summary(lm6)

# get predictive distribution
pred6 <- data.frame(S = seq(1, 6, 0.1))
pred6 <- dplyr::bind_cols(pred6, predict(lm6, pred6, interval = "confidence"))

# plot the results
r2 <- 
  ggplot() +
  geom_jitter(data = mf_df, mapping = aes(x = S, y = mf_ENFQ2), width = 0.1) +
  geom_line(data = pred6, mapping = aes(x = S, y = fit)) +
  geom_ribbon(data = pred6, mapping = aes(x = S, ymin = lwr, ymax = upr), 
              alpha = 0.1) +
  ylab("ENF-Q2 multifunctionality") +
  xlab(S) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.4,
           label = lm_eqn(lm6), parse = TRUE) +
  theme_meta()
plot(r2)

# combine the plots
r12 <- cowplot::plot_grid(r1, r2, align = "hv", ncol = 2, nrow = 1,
                          labels = c("a", "b"), label_fontface = "plain", label_size = 11)
plot(r12)

ggsave("figures-paper-2/fig_5.pdf", r12,
       units = "cm", width = 16, height = 8)

