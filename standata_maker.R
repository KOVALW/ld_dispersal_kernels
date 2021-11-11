nonr0_97data <- final_97data %>%
  filter(r > 0 | tree_spp == "BO") %>%
  mutate(r = ifelse(r == 0, 1e-3, r))

densities97 <- final_97data %>%
  group_by(indices, est_released) %>%
  summarise(obs =  n()) %>%
  ungroup() %>%
  mutate(rel_released = est_released / 3.8e5)

table_vector <- c()

for (a_vals in 1:2) {
  for (wind in 0:2) {
    for (distro in c("norm", "exp")) {
      table_name <-
        paste0("a", a_vals, "_w", wind, "_d", distro, "_input")
      
      if (a_vals == 1) {
        a_name <- paste0("drift_", distro, "_fit")
        a_pars <- c("log_b_coef", "log_a_intercept",  "lp__")
        a_indices <- rep(1, max(final_97data$indices))
      } else {
        a_name <- paste0("multi_", distro, "_fit")
        a_pars <-
          c("log_b_coef",
            "log_a_intercept[1]",
            "log_a_intercept[3]",
            "lp__")
        a_indices <- index_release$t
      }
      table_vector <- c(table_vector, table_name)
      
      lrs <- extract(get(a_name), pars = "log_resid_sigma")  %>%
        as.data.frame() %>%
        summarise(mu = mean(log_resid_sigma),
                  sdev = sd(log_resid_sigma))
      
      a_priors <- ret_pc(get(a_name), 1:(a_vals + 1), a_pars)
      
      
      
      if (wind == 0) {
        assign(
          table_name,
          list(
            N = nrow(nonr0_97data),
            r = nonr0_97data$r,
            n = nonr0_97data$larvae,
            leaf = nonr0_97data$leaves,
            theta = nonr0_97data$theta,
            indices = nonr0_97data$indices,
            rel_density = densities97$rel_released,
            c_order = ifelse(distro == "exp", 1, 2),
            K = max(final_97data$indices),
            tmt = rep(1, max(final_97data$indices)),
            qs = a_vals + 1,
            q_sds = c(a_priors$sdev * sd_scale, lrs$mu, lrs$sdev),
            q_sca = a_priors$scale,
            q_ctr = a_priors$center,
            q_rots = a_priors$rotation,
            drift = t(array(
              matrix(
                c(shifts_97$xshift,
                  shifts_97$yshift),
                nrow = 2,
                ncol = 6,
                byrow = T
              ),
              c(2, 6)
            ))
          )
        )
      } else if (wind == 1) {
        assign(
          table_name,
          list(
            N = nrow(nonr0_97data),
            r = nonr0_97data$r,
            n = nonr0_97data$larvae,
            leaf = nonr0_97data$leaves,
            theta = nonr0_97data$theta,
            indices = nonr0_97data$indices,
            K = max(final_97data$indices),
            qpoint = a_indices,
            windpoint = rep(1, max(final_97data$indices)),
            rel_density = densities97$rel_released,
            c_order = ifelse(distro == "exp", 1, 2),
            winds = 1,
            qs = a_vals + 1,
            wind_sds = array(aggwind_final_posteriors$sdev, c(1, 3)) * sd_scale,
            wind_sca = array(aggwind_final_posteriors$scale, c(1, 3)),
            wind_ctr = array(aggwind_final_posteriors$center, c(1, 3)),
            wind_rots = array(
              matrix(
                aggwind_final_posteriors$rotation,
                byrow = T,
                nrow = 3,
                ncol = 3
              ),
              c(1, 3, 3)
            ),
            q_sds = c(a_priors$sdev * sd_scale, lrs$mu, lrs$sdev),
            q_sca = a_priors$scale,
            q_ctr = a_priors$center,
            q_rots = a_priors$rotation
          )
        )
        
      } else {
        assign(
          table_name,
          list(
            N = nrow(nonr0_97data),
            r = nonr0_97data$r,
            n = nonr0_97data$larvae,
            leaf = nonr0_97data$leaves,
            theta = nonr0_97data$theta,
            indices = nonr0_97data$indices,
            K = max(final_97data$indices),
            qpoint = a_indices,
            windpoint = index_release$t,
            rel_density = densities97$rel_released,
            c_order = ifelse(distro == "exp", 1, 2),
            winds = 2,
            qs = a_vals + 1,
            wind_sds = t(array(
              c(
                windR1_final_posteriors$sdev,
                windR2_final_posteriors$sdev
              ),
              c(3, 2)
            )) * sd_scale,
            wind_sca = t(array(
              c(
                windR1_final_posteriors$scale,
                windR2_final_posteriors$scale
              ),
              c(3, 2)
            )),
            wind_ctr = t(array(
              c(
                windR1_final_posteriors$center,
                windR2_final_posteriors$center
              ),
              c(3, 2)
            )),
            wind_rots = w2a,
            q_sds = c(a_priors$sdev * sd_scale, lrs$mu, lrs$sdev),
            q_sca = a_priors$scale,
            q_ctr = a_priors$center,
            q_rots = a_priors$rotation
          )
        )
      }
    }
  }
}
