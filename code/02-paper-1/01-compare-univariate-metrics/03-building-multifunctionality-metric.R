#'
#' @title Building a multifunctionality metric
#' 
#' @description Outline six steps that we believe are important to think about
#' when definining a univariate metric of ecosystem multifunctionality
#' 

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(multifunc)
library(vegan)

# customised plotting theme
theme_meta1 <- 
  function(base_size = 12, base_family = "") {
    theme(panel.background = element_rect(fill = "white"), 
          panel.border = element_blank(),
          axis.line.x = element_line(color="black", size = 0.4),
          axis.line.y = element_line(color="black", size = 0.4),
          panel.grid.major =  element_blank(),
          panel.grid.minor =  element_blank(),
          axis.ticks.length = unit(-0.16, "cm"),
          axis.title.x = element_text(colour ="black", size = 9, face = "plain", margin=margin(3.5,0,0,0,"pt")),
          axis.title.y = element_text(colour = "black", size = 9, face = "plain", margin=margin(0,3.5,0,0,"pt")),
          axis.text.x = element_text(colour = "black", size=8, face = "plain",  margin=margin(3.5,0,0,0,"pt")),
          axis.text.y = element_text(colour ="black", size=8, face = "plain", margin=margin(0,3.5,0,0,"pt")),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.text = element_text(colour = "black", size=10, face = "plain"),
          legend.title = element_text(colour = "black", size=10, face = "plain"),
          legend.key = element_rect(fill = NA))
  }


# make a data.frame of hypothetical raw functions values
df <- dplyr::tibble(com = c(1, 2),
                    F1 = c(1.4, 6.7),
                    F2 = c(12, 5))
head(df)

# transform the two functions

# F1

# set the threshold
thresh <- 0.50
df$F1_t <- with(df, ifelse( F1 < (max(F1)*thresh), 0, 1 ))

# get additional points required for plotting the threshold curve
df_p1 <- dplyr::bind_rows(dplyr::select(df, F1, F1_t),
                          dplyr::tibble(F1 = rep(thresh*max(df$F1), 2),
                                        F1_t = c(0, 1) ))

# plot the the standardisation function for F1
p1 <- 
  ggplot(data = df_p1,
         mapping = aes(x = F1, y = F1_t)) +
  geom_line(colour = "#92D050", linewidth = 1) +
  ylab("F1-stand") +
  xlab("F1") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(limits = c(1.3, 6.8), 
                     breaks = c(1.5, 3.5, 6)) +
  theme_meta1()
plot(p1)

# output these two figures
ggsave(filename = "figures-paper-1/fig_3i.svg", p1,
       units = "cm", width = 3.5, height = 3.5)

# F2
df$F2_t <- with(df, F2/max(F2))
head(df)

# plot the standardisation function for F2
p2 <- 
  ggplot(data = df,
       mapping = aes(x = F2, y = F2_t)) +
  geom_line(colour = "#FFD966", linewidth = 1) +
  ylab("F2-stand") +
  xlab("F2") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(limits = c(4.9, 12.1), breaks = c(5, 8, 12)) +
  theme_meta1()
plot(p2)

# output these two figures
ggsave(filename = "figures-paper-1/fig_3ii.svg", p2,
       units = "cm", width = 3.5, height = 3.5)

# step 2: apply weightings
df$F1_t_w <- df$F1_t*0.7*2
df$F2_t_w <- df$F2_t*0.3*2

# step 3: calculate the effective number of functions

# get a sequence of q-parameters
q_pars <- seq(0, 2, 0.05)

# extract the relevant columns
df_in <- df[,c(6,7)]/rowSums(df[,c(6,7)])
names(df_in) <- c("pF1", "pF2")

q_df <- vector("list", length = length(q_pars))
for(i in 1:length(q_pars)) {
  
  ENF <- vegan::renyi(x = df_in, scales = q_pars[i], hill = TRUE)/ncol(df_in)
  x <- tibble(q = q_pars[i], com = c(1, 2), ENF = ENF)
  q_df[[i]] <- x
  
}

# bind into a data.frame
q_df <- bind_rows(q_df)

# plot the effective number of functions for com1
p1 <- 
  ggplot(data = dplyr::filter(q_df, com == 1) ,
         mapping = aes(x = q, y = ENF)) +
  geom_line(colour = "red", linewidth = 1) +
  ylab(expression("ENF  "^"q"*"N") ) +
  xlab("q-parameter") +
  scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_meta1()
plot(p1)

# output these two figures
ggsave(filename = "figures-paper-1/fig_3iii.svg", p1,
       units = "cm", width = 3.5, height = 3.5)

# plot the effective number of functions for com2
p2 <- 
  ggplot(data = dplyr::filter(q_df, com == 2) ,
         mapping = aes(x = q, y = ENF)) +
  geom_line(colour = "red", linewidth = 1) +
  ylab(expression("ENF  "^"q"*"N") ) +
  xlab("q-parameter") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  theme_meta1()
plot(p2)

# output these two figures
ggsave(filename = "figures-paper-1/fig_3iv.svg", p2,
       units = "cm", width = 3.5, height = 3.5)

# step 4: account for correlation

# with two communities, the correlation calculations don't really work...

# we just illustrate it

# add some negative error
q_df$ENF_cor <- q_df$ENF - 0.15

# plot the effective number of functions for com1
p3 <- 
  ggplot() +
  geom_line(data = dplyr::filter(q_df, com == 1) ,
            mapping = aes(x = q, y = ENF), 
            colour = "red", linewidth = 0.5, alpha = 0.5) +
  geom_line(data = dplyr::filter(q_df, com == 1) ,
            mapping = aes(x = q, y = ENF_cor), linewidth = 1, colour = "blue") +
  ylab(expression("ENF  "^"q"*"N") ) +
  xlab("q-parameter") +
  scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_meta1()
plot(p3)

# output these two figures
ggsave(filename = "figures-paper-1/fig_3v.svg", p3,
       units = "cm", width = 3.5, height = 3.5)

# plot the effective number of functions for com2
p4 <- 
  ggplot() +
  geom_line(data = dplyr::filter(q_df, com == 2) ,
            mapping = aes(x = q, y = ENF), 
            colour = "red", linewidth = 0.5, alpha = 0.5) +
  geom_line(data = dplyr::filter(q_df, com == 2) ,
            mapping = aes(x = q, y = ENF_cor), linewidth = 1, colour = "blue") +
  ylab(expression("ENF  "^"q"*"N") ) +
  xlab("q-parameter") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_x_continuous(limits = c(0, 2), breaks = c(0, 1, 2)) +
  theme_meta1()
plot(p4)

# output these two figures
ggsave(filename = "figures-paper-1/fig_3vi.svg", p4,
       units = "cm", width = 3.5, height = 3.5)

# calculate the effective number of functions with different qs
A <- rowSums(df[,c(6, 7)])/2

q_df |>
  dplyr::filter(q %in% c(0, 1, 2)) |>
  dplyr::mutate(A = rep(A, 3)) |>
  dplyr::mutate(MF = A*ENF)



