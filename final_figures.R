require(gridExtra)
require(ggfortify)
require(scico)
require(grid)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 


circle_fun <- function(center = c(0,0),diameter = 2, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

distort <- function(x, y, beta, gamma){
  term1 = (2*gamma^2 + 1)*((x/beta)^2 + 2*gamma*x/beta)
  term2 = (gamma^2 + 1)*(y^2 + 2*gamma^2)
  term3 = 2*gamma*(x/beta + gamma)
  rootterm = (gamma^2 + 1)*((x/beta + gamma)^2 + y^2 + 1)
  return(sqrt(term1+term2-term3*sqrt(rootterm)))
}

distort_distance<- function(x, y, beta, gamma, psi, dx = 0, dy = 0, dr = 0){
  if (dr) {
    x = x - dr * cos(psi)
    y = y - dr * sin(psi)
  } else if (dx & dy) {
    x = x - dx
    y = y - dy
  }
  x_adj = x*cos(psi) + y*sin(psi);
  y_adj = x*sin(psi) - y*cos(psi);
  return(distort(x_adj, y_adj, beta, gamma));
}


#Shape pca (wind and larval beta,gamma,psi posteriors)------------------------------------

posteriors_shape <- bind_rows(
  extract(windall_pcafit, pars = c("beta", "gamma", "psi", "lp__")) %>%
    as.data.frame() %>%
    mutate(release = "all_pca_wind"),
  extract(wind_fit, pars = c("beta", "gamma", "psi", "lp__")) %>%
    as.data.frame() %>%
    mutate(release = "all_unadj_wind"),
  extract(windR1_pcafit, pars = c("beta", "gamma", "psi", "lp__")) %>%
    as.data.frame() %>%
    mutate(release = "early_wind"),
  extract(windR2_pcafit, pars = c("beta", "gamma", "psi", "lp__")) %>%
    as.data.frame() %>%
    mutate(release = "late_wind"),
  extract(a1_w2_dnorm_fit, pars = c("beta", "gamma", "psi", "lp__")) %>%
    as.data.frame() %>%
    mutate(obs= row_number()) %>%
    gather(-lp__, -obs, key = "parm", value = "value") %>%
    mutate(par_name = substr(parm,1,nchar(parm)-2),
           par_id = substr(parm, nchar(parm), nchar(parm))) %>%
    select(-parm) %>%
    spread(par_name, value) %>%
    mutate(release = ifelse(par_id == 1, "early_larvae", "late_larvae")) %>%
    select(beta, gamma, psi, lp__, release),
  extract(a2_w1_dnorm_fit, pars = c("beta", "gamma", "psi", "lp__")) %>%
    as.data.frame() %>%
    mutate(release = "all_larvae")
)


indexed_lhoods <- order(posteriors_shape$lp__, decreasing = T)
shape_draws_ordered <- posteriors_shape[indexed_lhoods, 1:3]
shape_final_posteriors <- prcomp(shape_draws_ordered,scale=TRUE)

all_shape_plot <- autoplot(shape_final_posteriors, data = posteriors_shape[indexed_lhoods,],
         colour = "release", loadings = TRUE, loadings.colour = 'blue', 
         size = 2, shape = 15, loadings.label = TRUE, loadings.label.colour = 'black', 
         loadings.label.size = 5) +
  scale_color_manual(values = c("#963900", "#ff9959","#f25d00", "#350565","#6e2ab2","#00614b","#009c79"),
                     labels = c("Larvae Aggregate","Wind Aggregate", "Wind Unadj.", "Larvae Early", "Wind Early","Larvae Late", "Wind Late")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(col = "Model Set")

aggregate_shape_plot <- autoplot(shape_final_posteriors, data = posteriors_shape[indexed_lhoods,],
         colour = "release", loadings = TRUE, loadings.colour = 'blue', 
         size = 2, loadings.label = TRUE, loadings.label.colour = 'black', 
         loadings.label.size = 5) +
  scale_color_manual(values = c("#963900", "#ff9959","#f25d00", NA,NA,NA,NA),
                     labels = c("Larvae Aggregate","Wind Aggregate", "Wind Unadj.", "", "Wind Early","Larvae Late", "Wind Late")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(col = "Model Set")

timing_shape_plot <- autoplot(shape_final_posteriors, data = posteriors_shape[indexed_lhoods,],
         colour = "release", loadings = TRUE, loadings.colour = 'blue', 
         size = 2, loadings.label = TRUE, loadings.label.colour = 'black', 
         loadings.label.size = 5) +
  scale_color_manual(values = c(NA,NA,NA,#c("#963900", "#ff9959","#f25d00", NA,NA,NA,NA),#
                                "#350565","#6e2ab2","#00614b","#009c79"),
                     labels = c("Larvae Aggregate","Wind Aggregate", "Wind Unadj.", "Larvae Early", "Wind Early","Larvae Late", "Wind Late")) +
  theme_bw() +
  #theme(legend.position = c(1,1), legend.justification = c(1.5,3)) +
  theme(legend.position = "none") +
  labs(col = "Model Set")

legend <- g_legend(all_shape_plot) 

grid.arrange(
  grid.arrange(
    aggregate_shape_plot,
    timing_shape_plot,
    nrow = 1,
    ncol = 2
  ),
  grid.arrange(
    ggplot() + theme_void(),
    legend,
    ggplot() + theme_void(),
    nrow = 1,
    ncol = 3,
    widths = c(1, 8, 1)
  ),
  nrow = 2,
  ncol = 1,
  heights = c(10, 1)
)

#line fits for 2d kernel a1w2_dnorm----------------------------------------
best_posteriors <- extract(a1_w2_dnorm_fit, pars = c(paste0("a_intercept"),
                                                     "b_coef",
                                                     paste0("beta[",1:2,"]"),
                                                     paste0("gamma[",1:2,"]"),
                                                     paste0("psi[",1:2,"]"),
                                                     "lp__")) %>%
  as.data.frame()

best_posteriors[sample.int(nrow(best_posteriors),50),] -> set_of_parms
best_posteriors[which(best_posteriors$lp__ == max(best_posteriors$lp__)),] -> best_parms

mapping_df <- data.frame(theta = rep(sort(rep(seq(-pi,pi,length.out = 201),241)),6),
                         r = rep(rep(seq(0,60,length.out = 241),201),6),
                         indices = sort(rep(1:6,241*201))) %>%
  mutate(x = r * cos(theta),
         y = r * sin(theta))

vars_df <- data.frame(indices = 1:6,
                      a = as.numeric(c(rep(best_parms$a_intercept, 3),rep(best_parms$a_intercept, 3))),
                      b = as.numeric(best_parms$b_coef),
                      beta = as.numeric(c(rep(best_parms$beta.1., 3),rep(best_parms$beta.2., 3))),
                      gamma = as.numeric(c(rep(best_parms$gamma.1., 3),rep(best_parms$gamma.2., 3))),
                      psi = as.numeric(c(rep(best_parms$psi.1., 3),rep(best_parms$psi.2., 3))))

plane_plot_df <- final_97data %>%
  mutate(plane = ifelse(dir %in% c("S","N"), "SN",
                        ifelse(dir %in% c("SW","NE"), "SWNE",
                               ifelse(dir %in% c("SE", "NW"), "SENW",
                                      ifelse(dir %in% c("W","E"), "WE",
                                             dir))))) %>%
  bind_rows(final_97data %>% 
              filter(dir == "Center") %>%
              mutate(plane = "SN"),
            final_97data %>% 
              filter(dir == "Center") %>%
              mutate(plane = "SWNE"),
            final_97data %>% 
              filter(dir == "Center") %>%
              mutate(plane = "SENW"),
            final_97data %>% 
              filter(dir == "Center") %>%
              mutate(plane = "WE")) %>%
  mutate(plane_r = ifelse(dir %in% c("S", "SW","SE","W"), -r, r)) %>%
  filter(plane != "Center")

all_angles <- plane_plot_df %>%
  group_by(plane, dir, theta, indices, site, est_released) %>%
  summarise(obs = n()) %>%
  ungroup() %>%
  full_join(vars_df)

all_mapping_fits <- data.frame()

for (val_id in 1:nrow(set_of_parms)) {
    curr_parms <- set_of_parms[val_id,]
    
    vars_df <- data.frame(indices = 1:6,
                          a = as.numeric(c(rep(curr_parms$a_intercept, 3),rep(curr_parms$a_intercept, 3))),
                          b = as.numeric(curr_parms$b_coef),
                          beta = as.numeric(c(rep(curr_parms$beta.1., 3),rep(curr_parms$beta.2., 3))),
                          gamma = as.numeric(c(rep(curr_parms$gamma.1., 3),rep(curr_parms$gamma.2., 3))),
                          psi = as.numeric(c(rep(curr_parms$psi.1., 3),rep(curr_parms$psi.2., 3))))
    
    curr_df <- all_angles %>%
      select(-a, -b, -beta, -gamma, -psi) %>%
      full_join(vars_df)
      
    for (r in seq(0.5, 60, by = 0.5)) {
      all_mapping_fits <- bind_rows(all_mapping_fits,
      curr_df %>%
        mutate(dist = distort_distance(r * cos(theta), r * sin(theta), beta, gamma, psi)) %>%
        mutate(kde2 = a * est_released * exp(-b*dist^2)/(3.8e5 * 2*pi*beta*sqrt(gamma^2+1))) %>%
        mutate(plane_r = ifelse(dir %in% c("S", "SW","SE","W"), -r, r),
               par_set = val_id)
      )
    }
    print(val_id)
}

mapping_1d <- data.frame()

for(i in 1:nrow(all_angles)) {
  
  mapping_1d <- bind_rows(mapping_1d, 
                          bind_cols(all_angles[i,], 
                                    data.frame(r = seq(0,60, length.out = 200))) %>%
                            mutate(dist = distort_distance(r * cos(theta), r * sin(theta), beta, gamma, psi)) %>%
                            mutate(kde2 = a * est_released * exp(-b*dist^2)/(3.8e5 * 2*pi*beta*sqrt(gamma^2+1))) %>%
                            mutate(plane_r = ifelse(dir %in% c("S", "SW","SE","W"), -r, r))
  )
}

mapping_df <- full_join(mapping_df, vars_df) %>%
  mutate(dist = distort_distance(x, y, beta, gamma, psi)) %>%
  mutate(kde2 = a * exp(-b*dist^2)/(2*pi*beta*sqrt(gamma^2+1))) %>%
  mutate(density_tmt_id = (indices-1)%%3) %>%
  mutate(density_tmt = ifelse(density_tmt_id == 0, "Low Density", ifelse(density_tmt_id == 1, "Medium Density", "High Density")),
         release_label = ifelse(ceiling(indices/3) == 1, "Early Release", "Late Release")) %>%
  mutate(density_tmt = factor(density_tmt, levels = c("Low Density", "Medium Density", "High Density")))

mapping_1d <- mapping_1d %>%
  mutate(density_tmt_id = (indices-1)%%3) %>%
  mutate(density_tmt = ifelse(density_tmt_id == 0, "Low Density", ifelse(density_tmt_id == 1, "Medium Density", "High Density")),
         release_tmt = ifelse(ceiling(indices/3) == 1, "Early", "Late"),
         release_label = ifelse(ceiling(indices/3) == 1, "Early Release", "Late Release")) %>%
  mutate(density_tmt = factor(density_tmt, levels = c("Low Density", "Medium Density", "High Density")))

all_mapping_fits <- all_mapping_fits %>%
  mutate(density_tmt_id = (indices-1)%%3) %>%
  mutate(density_tmt = ifelse(density_tmt_id == 0, "Low Density", ifelse(density_tmt_id == 1, "Medium Density", "High Density")),
         release_tmt = ifelse(ceiling(indices/3) == 1, "Early", "Late"),
         release_label = ifelse(ceiling(indices/3) == 1, "Early Release", "Late Release")) %>%
  mutate(density_tmt = factor(density_tmt, levels = c("Low Density", "Medium Density", "High Density")))

plane_plot_df <- plane_plot_df %>%
  mutate(density_tmt_id = (indices-1)%%3) %>%
  mutate(release_tmt = ifelse(release_id == "A", "Early", "Late"))

distrofits_plot <- plane_plot_df %>%
  filter(!dir %in% c("Center", "SSE", "SSW", "ESE")) %>%
  full_join(vars_df) %>%
  mutate(plane.f = factor(plane, levels = c("SENW","SN", "NNE", "SWNE","ENE","WE"))) %>%
  #filter(release_label == "Late Release") %>%
  ggplot(aes(x = plane_r)) +
  geom_ribbon(aes(ymin = mink, ymax = maxk, fill = plane.f), alpha = 0.4, col = NA, data = all_mapping_fits %>%
              filter(plane %in% c("SENW","SN", "NNE", "SWNE","ENE","WE"))  %>%
                group_by(plane_r, plane, density_tmt_id, release_tmt) %>% 
                summarise(maxk = max(kde2),
                          mink = min(kde2)) %>% 
                ungroup() %>%
              mutate(plane.f = factor(plane, levels = c("SENW","SN", "NNE", "SWNE","ENE","WE")))
              ) +
  geom_line(aes(y = medkde, col = plane.f), lwd = 1.5, data  = all_mapping_fits %>% 
              filter(plane %in% c("SENW","SN", "NNE", "SWNE","ENE","WE")) %>%
              group_by(plane_r, plane, density_tmt_id, release_tmt) %>% 
              summarise(medkde = median(kde2)) %>% 
              ungroup() %>%
              mutate(plane.f = factor(plane, levels = c("SENW","SN", "NNE", "SWNE","ENE","WE")))) +
  geom_line(aes(y = kde2, col = plane.f), data = mapping_1d %>%
              filter(!dir %in% c("Center", "SSE", "SSW", "ESE")) %>%
              mutate(plane.f = factor(plane, levels = c("SENW","SN", "NNE", "SWNE","ENE","WE")))) +
  geom_point(aes(y = larvae/leaves, fill = plane.f), pch = 21) +
  facet_grid(interaction(density_tmt_id, release_tmt,  sep = "' ")~plane.f) +
  geom_vline(xintercept = 0, lty = 3) +
  scale_color_scico_d(palette = "lapaz", end = 0.8) +
  scale_fill_scico_d(palette = "lapaz", end = 0.8) +
  theme_bw() +
  theme(legend.position = "none") + 
  labs(x = "Distance from center (m)", y = "Density (larvae/leaf)")

mapping_plot <- mapping_df %>%
  filter(indices %in% c(1,4)) %>%
  mutate(theta = ifelse(theta < 0, theta + 2 * pi, theta)) %>%
  ggplot(aes(x = theta, y = r)) + 
  geom_contour_filled(aes(z = kde2),bins = 15) + 
  facet_wrap(release_label~., ncol = 1) + 
  coord_polar() + 
  scale_fill_scico_d(palette = "grayC") +
  theme_minimal() +
  labs(y = "Distance from center (m)", x = "") +
#  theme(legend.position = "none",
#        axis.title.x=element_blank(),
#        axis.text.x=element_blank(),
#        axis.ticks.x=element_blank()) +
  geom_point(data = plane_plot_df %>% 
               filter(indices %in% c(1,4)) %>%
               mutate(theta = ifelse(theta < (pi/2), 2 * pi + (theta + 3 * pi / 2) * -1, 2*pi + (theta - pi / 2) * -1)) %>%
               mutate(plane.f = factor(plane, levels = c("SENW","SN", "NNE", "SWNE","ENE","WE"))) %>%
               filter(!is.na(plane.f)), size = 3.5) +
  # geom_point(data = final_97data %>% 
  #              filter(indices %in% c(1,4)) %>%
  #              mutate(theta = ifelse(theta < 0, theta + 2 * pi, theta)), pch = 4, col = "white)
  geom_point(aes(color = plane.f), data = plane_plot_df %>% 
               filter(indices %in% c(1,4)) %>%
               mutate(theta = ifelse(theta < (pi/2), 2 * pi + (theta + 3 * pi / 2) * -1, 2*pi + (theta - pi / 2) * -1)) %>%
               mutate(plane.f = factor(plane, levels = c("SENW","SN", "NNE", "SWNE","ENE","WE"))), size = 3) +
  geom_contour(aes(z = kde2),bins = 5, col = "white") + 
  scale_color_scico_d(palette = "lapaz", end = 0.8) +
  theme(legend.position = "none")

grid.arrange(mapping_plot, distrofits_plot, nrow = 1, ncol = 2, widths = c(0.4,0.6))

#1 m distortion------------------------------------------
ctr_fxn <- function(centers, release_name, nexp = T) {
  df <- centers %>%
    as.data.frame() %>%
    rownames_to_column(var = "par_name") %>%
    rename(value = ".") %>%
    spread(par_name, value) %>%
    mutate(release = release_name) 
  
  if (nexp) {
    df <- df %>%
      mutate(beta = exp(logbeta),
             gamma = exp(loggamma))
  }
  
  df %>% select(beta, gamma, psi, release)
}

a1_w2_dnorm_pca <- ret_pc(a1_w2_dnorm_fit, 1:6, c("beta", "gamma", "psi", "lp__"))
a1_w1_dnorm_pca <- ret_pc(a1_w1_dnorm_fit, 1:3, c("beta", "gamma", "psi", "lp__"))

ctrs_all <- bind_rows(
a1_w2_dnorm_pca$center %>%
  as.data.frame() %>%
  rownames_to_column(var = "parm") %>%
  rename(value = ".") %>%
  mutate(par_name = substr(parm,1,nchar(parm)-2),
         par_id = substr(parm, nchar(parm), nchar(parm))) %>%
  select(-parm) %>%
  spread(par_name, value) %>%
  mutate(release = ifelse(par_id == 1, "early_larvae", "late_larvae")) %>%
  select(beta, gamma, psi, release),
ctr_fxn(a1_w1_dnorm_pca$center, "agg_larvae", F),
ctr_fxn(aggwind_final_posteriors$center, "agg_wind"),
ctr_fxn(windR1_final_posteriors$center, "early_wind"),
ctr_fxn(windR2_final_posteriors$center, "late_wind")
)

contouring_data <- data.frame()

for (ctr in 1:nrow(ctrs_all)){
  curr_ctr <- ctrs_all[ctr,]
  
  contouring_data <- bind_rows(contouring_data,
  data.frame(x = sort(rep(seq(-3,6,length.out = 500),500)),
             y = rep(seq(-3,6, length.out = 500), 500),
             beta = curr_ctr$beta, 
             gamma = curr_ctr$gamma,
             psi = curr_ctr$psi,
             release = curr_ctr$release) %>%
    mutate(dist = distort_distance(x,y,beta,gamma,psi))
  )
}

contouring_data %>%
  filter(abs(dist-1)<1e-2) %>%
  arrange(angle = atan2(y,x)) %>%
  ggplot(aes(x = x, y = y, col = release)) +
  theme_bw() +
  #scale_color_manual(values = c("#ff865d","grey", "#193893", "#aa2e03")) +
  geom_point(aes(x = U, y = V), data = bind_rows(wind_data %>% mutate(release = ifelse(release == 1, "early_wind", "late_wind"))),alpha = 0.5) +
  geom_path(lwd = 1.5) +
  scale_color_manual(values = c("#963900", "#ff9959", "#350565","#6e2ab2","#00614b","#009c79"),
                     labels = c("Larvae Aggregate","Wind Aggregate", "Larvae Early", "Wind Early","Larvae Late", "Wind Late"))

dual_contour_plot <- contouring_data %>%
  mutate(facet_release = ifelse(release %in% c("agg_larvae", "agg_wind"), "Aggregate", "Temporal")) %>%
  filter(abs(dist-1)<1e-2) %>%
  arrange(angle = atan2(y,x)) %>%
  ggplot(aes(x = x, y = y, col = release)) +
  theme_bw() +
  facet_wrap(facet_release~.) +
  theme(legend.position = "none") +
  #scale_color_manual(values = c("#ff865d","grey", "#193893", "#aa2e03")) +
  geom_point(aes(x = U, y = V), data = bind_rows(wind_data %>% mutate(release = ifelse(release == 1, "early_wind", "late_wind"),
                                                                      facet_release = "Temporal"),
                                                 wind_data %>% mutate(release = "agg_wind",
                                                                      facet_release = "Aggregate")),alpha = 0.25) +
  geom_path(lwd = 1.5) +
  scale_color_manual(values = c("#963900", "#ff9959", "#350565","#6e2ab2","#00614b","#009c79"),
                     labels = c("Larvae Aggregate","Wind Aggregate", "Larvae Early", "Wind Early","Larvae Late", "Wind Late")) +
  geom_path(col = 1, data = circle_fun(),lty =2,lwd = 1.5,alpha = 0.75) +
  labs(x = "Distance East (m)", y = "Distance North (m)")

grid.arrange(
  dual_contour_plot,
  grid.arrange(
    aggregate_shape_plot,
    timing_shape_plot,
    nrow = 1,
    ncol = 2
  ),
  grid.arrange(
    ggplot() + theme_void(),
    legend,
    ggplot() + theme_void(),
    nrow = 1,
    ncol = 3,
    widths = c(1, 8, 1)
  ),
  nrow = 3,
  ncol = 1,
  heights = c(5,5, 1)
)

#grid distortion----------------------------------------------------------
mapping_xy <- data.frame(x = rep(sort(rep(seq(-50,50,length.out = 201),241)),6),
                         y = rep(rep(seq(-60,60,length.out = 241),201),6),
                         indices = sort(rep(1:6,241*201)))
mapping_xy <- full_join(mapping_xy, vars_df) %>%
  mutate(dist = distort_distance(x, y, beta, gamma, psi)) %>%
  mutate(kde2 = a * exp(-b*dist^2)/(2*pi*beta*sqrt(gamma^2+1))) %>%
  mutate(density_tmt_id = (indices-1)%%3) %>%
  mutate(density_tmt = ifelse(density_tmt_id == 0, "Low Density", ifelse(density_tmt_id == 1, "Medium Density", "High Density")),
         release_label = ifelse(ceiling(indices/3) == 1, "Early Release", "Late Release")) %>%
  mutate(density_tmt = factor(density_tmt, levels = c("Low Density", "Medium Density", "High Density")))


grid.arrange(
  mapping_xy %>%
    filter((x %% 5) == 0 & (y %% 5)  == 0 & release_label == "Early Release" & density_tmt == "Low Density") %>%
    ggplot(aes(x = x, y = y)) +
    geom_path(aes(group = y, col = y)) +
    geom_path(aes(group = x, col = x)) +
    theme_minimal() +
    geom_point(aes(x = x_r, y = y_r), data = final_97data %>%
                 filter(density_label == "Low Density", release_label == "Early Release"),
               size = 2,alpha = 0.75) +
    geom_vline(xintercept = 0, lty = 3) +
    geom_hline(yintercept = 0, lty = 3)  +
    theme(legend.position = "none") +
    labs(x= "Easting (m)", y = "Northing (m)"),
  
  mapping_xy %>%
    mutate(theta = atan2(y,x)) %>% 
    filter((x %% 5) == 0 & (y %% 5)  == 0 & release_label == "Early Release" & density_tmt == "Low Density") %>%
    ggplot(aes(x = dist * cos(theta), y = dist * sin(theta))) +
    geom_path(aes(group = y, col = y)) +
    geom_path(aes(group = x, col = x)) +
    theme_minimal() +
    geom_point(aes(x = true_dist * cos(theta), y = true_dist * sin(theta)), data = final_97data %>%
                 filter(density_label == "Low Density", release_label == "Early Release") %>%
                 mutate(true_dist = distort_distance(x_r, y_r, best_parms$beta.1.,
                                                     best_parms$gamma.1.,
                                                     best_parms$psi.1.)),
               size = 2,alpha = 0.75) +
    geom_vline(xintercept = 0, lty = 3) +
    geom_hline(yintercept = 0, lty = 3)  +
    theme(legend.position = "none") +
    labs(x = "Distorted distance (u(x,y))", y = "Distorted distance (v(x,y))"),
  
  nrow = 2, ncol = 1
)
