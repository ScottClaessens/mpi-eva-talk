library(ape)       # 5.8-1
library(brms)      # 2.22.0
library(btw)       # 2.0
library(coevolve)  # 0.1.0.9008
library(ggplot2)   # 3.5.2
library(ggtree)    # 3.14.0
library(tidybayes) # 3.0.7
library(tidyverse) # 2.0.0

# get tree for title page
set.seed(1)
title_tree <-
  ggtree(ape::rcoal(60)) +
  ggtree::geom_tree(colour = "white", size = 0.5) +
  theme(
    panel.background = element_rect(fill = "black", colour = "black"),
    plot.background = element_rect(fill = "black", colour = "black")
  )
ggsave(
  file = "assets/title_tree.png",
  plot = title_tree,
  width = 2,
  height = 3
)

# simulate data
# x -> y but not vice versa
selection_matrix <- 
  matrix(
    data = c(0.8, 0.8, 0, 0.8),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c("x", "y"),
      c("x", "y")
    )
  )
sim <-
  coev_simulate_coevolution(
    n = 100,
    variables = c("x", "y"),
    selection_matrix = selection_matrix,
    drift = c(x = 1, y = 1),
    prob_split = 0.05
  )
# standardise data
sim$data$x <- as.numeric(scale(sim$data$x))
sim$data$y <- as.numeric(scale(sim$data$y))

# plot simulated data
plot_sim_data <-
  ggplot(
    data = sim$data,
    mapping = aes(
      x = x,
      y = y
    )
  ) +
  geom_point() +
  labs(
    x = "Agriculture (std)",
    y = "Political complexity (std)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 14)
  )
ggsave(
  filename = "assets/simulated_data.png",
  plot = plot_sim_data,
  height = 4,
  width = 3.5
)

# plot simulated tree
ggsave(
  file = "assets/simulated_tree.png",
  plot = ggtree(sim$tree),
  width = 3,
  height = 4
)

# phylogenetic GLMM 1: x -> y
m1 <-
  brm(
    formula = y ~ x + (1 | gr(species, cov = A)),
    data = sim$data,
    data2 = list(A = ape::vcv.phylo(sim$tree, corr = TRUE)),
    prior = c(
      prior(normal(0, 0.5), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(2), class = sd),
      prior(exponential(2), class = sigma)
    ),
    cores = 4,
    chains = 4,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    file = "m1.rds",
    seed = 2113
  )
post <- posterior_samples(m1)
plot_glmm1 <-
  tibble(slope = post$b_x) |>
  ggplot(aes(x = slope)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  stat_slabinterval() +
  xlab("Effect of agriculture on political complexity") +
  ylab(NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12)
  )
ggsave(
  plot = plot_glmm1,
  file = "assets/glmm1.png",
  height = 3,
  width = 4
)

# phylogenetic GLMM 2: y -> x
m2 <-
  brm(
    formula = x ~ y + (1 | gr(species, cov = A)),
    data = sim$data,
    data2 = list(A = ape::vcv.phylo(sim$tree, corr = TRUE)),
    prior = c(
      prior(normal(0, 0.5), class = Intercept),
      prior(normal(0, 0.5), class = b),
      prior(exponential(2), class = sd),
      prior(exponential(2), class = sigma)
    ),
    cores = 4,
    chains = 4,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    file = "m2.rds",
    seed = 2113
  )
post <- posterior_samples(m2)
plot_glmm2 <-
  tibble(slope = post$b_y) |>
  ggplot(aes(x = slope)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  stat_slabinterval() +
  xlab("Effect of political complexity on agriculture") +
  ylab(NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 12)
  )
ggsave(
  plot = plot_glmm2,
  file = "assets/glmm2.png",
  height = 3,
  width = 4
)

# mean split continuous data (character for bayes traits)
d <- data.frame(id = sim$data$species)
d$x <- as.character(as.numeric(sim$data$x > mean(sim$data$x)))
d$y <- as.character(as.numeric(sim$data$y > mean(sim$data$y)))

# discrete model in bayes traits - independent
bt_ind <- 
  bayestraits(
    data = d,
    tree = sim$tree,
    commands = c(
      "2",               # independent model
      "2",               # MCMC
      "PriorAll exp 10", # set priors
      "Stones 100 1000", # stepping stone sampler
      "it 5000000"       # number of iterations
    ),
    silent = FALSE
  )
bt_ind$Stones$logMarLH # -140.4435

# discrete model in bayes traits - dependent
bt_dep <-
  bayestraits(
    data = d,
    tree = sim$tree,
    commands = c(
      "3",               # dependent model
      "2",               # MCMC
      "PriorAll exp 10", # set priors
      "Stones 100 1000", # stepping stone sampler
      "it 5000000"       # number of iterations
    ),
    silent = FALSE
  )
bt_dep$Stones$logMarLH # -120.9033

# rates
bt_dep$Log$results |>
  as_tibble() |>
  summarise(across(q12:q43, median))

# log bayes factor
2 * (bt_dep$Stones$logMarLH - bt_ind$Stones$logMarLH) # 39.08038

# coevolve
d <- 
  sim$data |>
  transmute(
    id = species,
    agr = x,
    pol = y
  )
tree <- sim$tree
m3 <- coev_fit(
  data = d,
  variables = list(
    agr = "normal",
    pol = "normal"
  ),
  id = "id",
  tree = tree,
  estimate_correlated_drift = FALSE,
  chains = 4,
  parallel_chains = 4,
  seed = 1
)
save_coevfit(m3, "m3.rds")
m3 <- readRDS("m3.rds")
summary(m3)

# plot delta theta
plot_delta_theta <- coev_plot_delta_theta(m3)
ggsave(
  filename = "assets/delta_theta.png",
  plot = plot_delta_theta,
  height = 3.5,
  width = 5
)

# plot selection gradient
plot_selection_matrix <- 
  coev_plot_selection_gradient(
    object = m3,
    var1 = "agr",
    var2 = "pol"
  )
ggsave(
  filename = "assets/selection_gradient.png",
  plot = plot_selection_matrix,
  height = 3.5,
  width = 5
)

# plot predictive series
plot_pred_series <-
  coev_plot_pred_series(
    object = m3,
    eta_anc = list(agr = 0, pol = 0),
    intervention_values = list(agr = 2, pol = NA)
  )
ggsave(
  filename = "assets/pred_series.png",
  plot = plot_pred_series,
  height = 4,
  width = 5.5
)
