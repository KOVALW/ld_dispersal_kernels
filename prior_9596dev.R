require(tidyverse)
require(rstan)
require(loo)

options(mc.cores = parallel::detectCores()) #run cores in parallel
rstan_options(auto_write = TRUE) #auto-write compiled code to hard drive

loo_extract <- function(stan_obj) {
  log_lik <- extract_log_lik(stan_obj,
                             merge_chains = F)
  r_eff <- relative_eff(exp(log_lik))
  loo_obj <- loo(log_lik, r_eff = r_eff)
  return(loo_obj)
}

ret_pc <- function(pcafit, index, namelist) {
  pcadf <- extract(pcafit, pars = namelist) %>%
    as.data.frame()
  
  indexed_lhoods <- order(pcadf$lp__, decreasing = T)
  ordered <- pcadf[indexed_lhoods, index]
  posteriors <- prcomp(ordered, scale = T)
  
  return(posteriors)
} 

qlogtransf <- function(q) {
  log(q / (1 - q))
}

invqtransf <- function(q_T) {
  exp(q_T) / (1 + exp(q_T))
}

#PCA over wind variables---------------------------------------------
wind_data <- read.csv("./data_shared/climate/wind/wind_AH_1997.csv") %>%
  mutate(U = as.numeric(U),
         V = as.numeric(V)) %>%
  drop_na() %>%
  filter(WS > 0) %>%
  mutate(release = ifelse(day > 145, 2, 1))


wind_stan <- list(
  N = nrow(wind_data),
  u = wind_data$U,
  v = wind_data$V#,
  #indices = fitdata$release_no
) 

wind_init <- function() {
  list(
    beta = runif(1,-5,5), 
    gamma = runif(1,-5,5),
    psi = runif(1,-pi*1/6, pi*2/3)
  )
}

wind_fit <- stan(file = "./code/final_models/windv0.stan", data = wind_stan, init = wind_init, iter = 3000)

wind_draws <- extract(wind_fit, pars = c("logbeta", "loggamma", "psi", "lp__")) %>% as.data.frame()

indexed_lhoods <- order(wind_draws$lp__, decreasing = T)
wind_draws_ordered <- wind_draws[indexed_lhoods, 1:3]
wind_pc_init <- prcomp(wind_draws_ordered,scale=TRUE)

mag_diff_mean <- mean(sqrt(wind_draws$logbeta^2 + wind_draws$loggamma^2))
#mag_diff_sd <- sqrt(mean((wind_draws$logbeta - wind_draws$loggamma)^2) - mean(wind_draws$logbeta - wind_draws$loggamma)^2)
mag_diff_var <- mean(wind_draws$logbeta^2 + wind_draws$loggamma^2) - mean(sqrt(wind_draws$logbeta^2 + wind_draws$loggamma^2))^2


scale_diff <- mag_diff_var / mag_diff_mean
shape_diff <- mag_diff_mean / scale_diff
# used_par <- c()
# pc_par <- wind_pc_init$x[1,]
# 
# for (j in 1:3) {
#   curr = wind_pc_init$rotation[j,];
#   used_par[j] = wind_pc_init$center[j] + (pc_par %*% curr) * wind_pc_init$scale[j]
# }
# 
# 
# scaling_factor <- 1.1

#running pca fit routine-----------------------------------------
windR1_pca_stan <- list(
  N = wind_data %>% filter(release == 1) %>% nrow(),
  u = wind_data$U[which(wind_data$release == 1)],
  v = wind_data$V[which(wind_data$release == 1)],
  P = 3,
  rot = wind_pc_init$rotation,
  sca = wind_pc_init$scale,
  ctr = wind_pc_init$center,
  mag_penalty = c(shape_diff,
                  1/scale_diff)
  #  sds = wind_pc_init$sdev * scaling_factor
)

windR1_pcafit <- stan(file = "./code/final_models/wind_pca.stan", 
                      data = windR1_pca_stan, iter = 3000)

windR2_pca_stan <- list(
  N = wind_data %>% filter(release == 2) %>% nrow(),
  u = wind_data$U[which(wind_data$release == 2)],
  v = wind_data$V[which(wind_data$release == 2)],
  P = 3,
  rot = wind_pc_init$rotation,
  sca = wind_pc_init$scale,
  ctr = wind_pc_init$center,
  mag_penalty = c(shape_diff,
                  1/scale_diff)
  #  sds = wind_pc_init$sdev * scaling_factor
)

windR2_pcafit <- stan(file = "./code/final_models/wind_pca.stan", 
                      data = windR2_pca_stan, iter = 3000)

windall_tmt_pca_stan <- list(
  N = wind_data %>% nrow(),
  u = wind_data$U,
  v = wind_data$V,
  release = wind_data$release,
  P = 3,
  winds = max(wind_data$release),
  wind_rots = wind_pc_init$rotation,
  wind_sca = wind_pc_init$scale,
  wind_ctr = wind_pc_init$center,
  mag_penalty = c(shape_diff,
                  1/scale_diff)
  #  sds = wind_pc_init$sdev * scaling_factor
)

windall_tmt_pcafit <- stan(file = "./code/final_models/wind_pca_releases.stan", 
                           data = windall_tmt_pca_stan, iter = 3000)
# 
windall_pca_stan <- list(
  N = nrow(wind_data),
  u = wind_data$U,
  v = wind_data$V,
  P = 3,
  rot = wind_pc_init$rotation,
  sca = wind_pc_init$scale,
  ctr = wind_pc_init$center,
  mag_penalty = c(shape_diff,
                  1/scale_diff)
  #sds = wind_pc_init$sdev * scaling_factor
)

windall_pcafit <- stan(file = "./code/final_models/wind_pca.stan", 
                       data = windall_pca_stan, iter = 3000)

if (compare_loos){
  pca_loo <- loo_extract(windall_pcafit)
  corr_loo <- loo_extract(wind_fit)
  tmt_loo <- loo_extract(windall_tmt_pcafit)
  windpca_comp <- loo_compare(pca_loo, corr_loo, tmt_loo)
}


aggwind_final_posteriors <- ret_pc(windall_pcafit, 1:3, c("logbeta", "loggamma", "psi", "lp__"))
windR1_final_posteriors <- ret_pc(windR1_pcafit, 1:3, c("logbeta", "loggamma", "psi", "lp__"))
windR2_final_posteriors <- ret_pc(windR2_pcafit, 1:3, c("logbeta", "loggamma", "psi", "lp__"))

wr1 <- matrix(windR1_final_posteriors$rotation, nrow = 3, ncol = 3)
wr2 <- matrix(windR2_final_posteriors$rotation, nrow = 3, ncol = 3)
w2a <- array(NA, c(2,3,3))
w2a[1,,] <- wr1
w2a[2,,] <- wr2

#97data-------------------
#Loading 1997 data---------------------
post_data_correct <- read_csv("./data_shared/literature/Hunter_2001_97_corrected.csv") %>%
  mutate(dir = ifelse(dir %in% c("Cent", "Center","Cen","CENTER"), "Center", dir))
post_data_checked <- read_csv("./data_shared/literature/Hunter_2001_97_checkedsites.csv")
reldens <- read_csv("./data_shared/literature/Hunter_2001_reldens.csv")
angle_dictionary <- data.frame(dir = c("Center", "E", "ENE", 
                                       "ESE", "N", "NNE", "NE", 
                                       "NW", "S", "SE", 
                                       "SSE", "SSW", "SW", "W"),
                               theta = c(0, 0, pi / 8, -pi/8,
                                         pi / 2, 3 * pi / 8, pi / 4,
                                         3 * pi / 4, -pi / 2, 
                                         -pi / 4, -3 * pi / 8, 
                                         -5 * pi / 8, -3 * pi / 4, pi))
site_names <- unique(post_data_correct$site)
all_checked <- data.frame()

for (i in 1:length(site_names)) {
  
  all_checked <- bind_rows(all_checked, data.frame(site = site_names[i],
                                                   r = post_data_checked$r,
                                                   dir = post_data_checked$dir))
}
zero_adj_data <- full_join(all_checked, post_data_correct) %>%
  mutate(leaves = ifelse(is.na(larvae), 1, leaves)) %>%
  mutate(larvae = ifelse(is.na(larvae), 0, larvae)) %>%
  full_join(angle_dictionary) %>%
  mutate(r = ifelse(is.na(r) & dir == "Center", 0, r)) %>%
  mutate(x_r = r * cos(theta),
         y_r = r * sin(theta)) 

final_97data <- zero_adj_data %>%
  mutate(density_tmt_id = substr(site,2,2),
         release_id = substr(site,1,1)) %>%
  mutate(indices = ifelse(release_id == "A", 
                          as.numeric(density_tmt_id) - 2, 
                          as.numeric(density_tmt_id) + 1)) %>%
  left_join(reldens) %>%
  group_by(release_id) %>%
  mutate(density_label = ifelse(est_released == min(est_released), "Low Density", ifelse(est_released == max(est_released), "High Density", "Medium Density"))) %>%
  ungroup() %>%
  mutate(release_label = ifelse(release_id == "A", "Early Release", "Late Release")) %>%
  mutate(density_label = factor(density_label, levels = c("Low Density", "Medium Density", "High Density")))

index_release <- final_97data %>% 
  group_by(indices, release_id) %>% 
  summarise(obs = n()) %>% 
  ungroup() %>% 
  arrange(indices) %>% 
  mutate(t = ifelse(release_id == "A", 1, 2)) %>% 
  select(t)

nonr0_97data <- final_97data %>%
  filter(r > 0 | tree_spp == "BO")

shifts_97 <- nonr0_97data %>%
  group_by(indices) %>%
  summarise(xshift = weighted.mean(r * cos(theta), larvae),
            yshift = weighted.mean(r * sin(theta), larvae)) %>%
  ungroup()

rexpect <- final_97data %>% 
  group_by(indices) %>% 
  summarise(expect_x = sum(x_r * larvae / leaves) / sum(larvae / leaves),
            expect_y = sum(y_r * larvae / leaves) / sum(larvae / leaves)) %>%
  ungroup() %>%
  mutate(r = sqrt(expect_x^2 + expect_y^2))

rexpect_scale <-(sd(rexpect$r) * 0.5)^2 / mean(rexpect$r)
rexpect_shape <- mean(rexpect$r) / rexpect_scale


shifts_97 <- nonr0_97data %>%
  group_by(indices) %>%
  summarise(xshift = weighted.mean(r * cos(theta), larvae),
            yshift = weighted.mean(r * sin(theta), larvae)) %>%
  ungroup()

sd_scale <- 1.5

#9596----------------------------------------------------

prior_data <- read_csv("data_shared/literature/Hunter_2001_95-96.csv")

shifts_9596 <- prior_data %>%
  group_by(release_no) %>%
  summarise(xshift = weighted.mean(r * cos(theta), observations),
            yshift = weighted.mean(r * sin(theta), observations)) %>%
  ungroup()

abc_prior9596_stan <- list(
  N = nrow(prior_data),
  n = prior_data$observations,
  r = prior_data$r,
  theta = prior_data$theta,
  leaf = prior_data$leaves,
  indices = prior_data$release_no,
  K = length(unique(prior_data$release_no)),#,
  rel_density = c(5.0, 5.7, 4.6,4.1,5.1,3.8)/min(c(5.0, 5.7, 4.6,4.1,5.1,3.8)),
  drift = array(matrix(c(shifts_9596$yshift, shifts_9596$xshift),
                       nrow = 6,
                       ncol = 2,
                       byrow = F), dim = c(6,2)),
  standard_N = min(c(5.0, 5.7, 4.6,4.1,5.1,3.8))*1e5,
  tmt = c(1,2,3,1,2,3)
)

drift_norm_fit <- stan(file = "./code/final_models/abc_Taylor1978.stan",
                       data = c(abc_prior9596_stan,
                                c_order = 2))
drift_exp_fit <- stan(file = "./code/final_models/abc_Taylor1978.stan",
                       data = c(abc_prior9596_stan,
                                c_order = 1))

multi_exp_fit <- stan(file = "./code/final_models/multia_Taylor1978.stan",
                      data = c(abc_prior9596_stan,
                               c_order = 1,
                               M = 3))
multi_norm_fit <- stan(file = "./code/final_models/multia_Taylor1978.stan",
                      data = c(abc_prior9596_stan,
                               c_order = 2,
                               M = 3))

source("./code/R/standata_maker.R")

model_vector <- c()

for (tbl in table_vector) {
  model_name <- paste0(substr(tbl,1,nchar(tbl) - 5), "fit")
  model_vector <- c(model_vector, model_name)
  if(substr(tbl,4,5) == "w0") {
    assign(model_name, stan(file = "./code/final_models/pca_Taylor1978.stan",
                            data = get(tbl)))
  } else {
    assign(model_name, stan(file = "./code/final_models/vp_Taylor1978.stan",
                            data = get(tbl)))
  }
  print(paste0("successful fit of ", model_name))
}

loo_vector<- c()

for (mdl in model_vector) {
  assign(paste0(substr(mdl,1, nchar(mdl) - 3), "loo"),
         loo_extract(get(mdl)))
  loo_vector <- c(loo_vector, paste0(substr(mdl,1, nchar(mdl) - 3), "loo"))
  print(paste0(substr(mdl,1, nchar(mdl) - 3), "loo"))
}

loo_compare(a1_w0_dnorm_loo,
            a1_w1_dnorm_loo,
            a1_w2_dnorm_loo,
            a2_w0_dnorm_loo,
            a2_w1_dnorm_loo,
            a2_w2_dnorm_loo,
            a1_w0_dexp_loo,
            a1_w1_dexp_loo,
            a1_w2_dexp_loo,
            a2_w0_dexp_loo,
            a2_w1_dexp_loo,
            a2_w2_dexp_loo)


# c1_loo <- loo_extract(drift_exp_fit)
# c2_loo <- loo_extract(drift_norm_fit)
# c1a3_loo <- loo_extract(multi_exp_fit)
# c2a3_loo <- loo_extract(multi_norm_fit)
# 
# loo_compare(c1_loo, c2_loo,c1a3_loo,c2a3_loo)
# 
# final_97data %>%
#   full_join(shifts_97) %>%
#   mutate(obs_r = sqrt((x_r - xshift)^2 + (y_r - yshift)^2)) %>%
#   mutate(expect_norm = est_released * 0.07282321 * exp(-0.002874959 * obs_r^2) / 3.8e5,
#          expect_dexp = est_released * 0.09886729 * exp(-0.07277637 * obs_r) / 3.8e5) %>%
#   ggplot(aes(x = obs_r)) +
#   geom_point(aes(y = larvae / leaves, col = indices)) +
#   facet_wrap(indices~.) +
#   geom_line(aes(y = expect_dexp, col = indices, group = indices))
# 
# final_97data %>%
#   full_join(shifts_97) %>%
#   mutate(obs_r = sqrt((x_r - xshift)^2 + (y_r - yshift)^2)) %>%
#   mutate(expect_norm = est_released * 0.07282321 * exp(-0.002874959 * obs_r^2) / 3.8e5,
#          expect_dexp = est_released * 0.09886729 * exp(-0.07277637 * obs_r) / 3.8e5) %>%
#   mutate(residual_norm = (expect_norm - (larvae/leaves)),
#          residual_dexp = (expect_dexp - (larvae/leaves))) %>%
#   ggplot(aes(x = theta, y = residual_norm)) +
#   geom_hline(yintercept = 0) +
#   geom_point() +
#   geom_smooth() +
#   coord_polar()
# 
# 
