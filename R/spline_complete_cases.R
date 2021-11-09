

library(table.glue)
library(tidyverse)
library(splines)
library(geepack)
library(mice)

source("R/functions.R")

data_nhanes <- table.glue::nhanes |>
 as_tibble() |>
 select(age, bp_sys_mmhg, bp_dia_mmhg, meds_bp) |>
 filter(meds_bp == 'yes', age >= 20) |>
 drop_na() |>
 mutate(
  bp_control = if_else(
   condition = bp_sys_mmhg < 140 & bp_dia_mmhg < 90,
   true = 1,
   false = 0
  ),
  id = seq(n())
 )

fit <- geeglm(bp_control ~ ns(age, df = 4),
              data = data_nhanes,
              family = 'poisson',
              id = id)

basis <- ns(data_nhanes$age, df = 4)

spline_preds <- get_spline_preds(fit, basis,
                                 pattern = '^ns\\(',
                                 x_min = 20,
                                 x_max = 80,
                                 x_ref = 50)

data_segment <- bin_segments(x = data_nhanes$age,
                             y = data_nhanes$bp_control,
                             x_min = 20,
                             x_max = 80,
                             by_y = TRUE,
                             bin_length = 2/3,
                             bin_count = 15,
                             bin_yintercept = 0.45) |>
 mutate(event_status = factor(event_status,
                              levels = c(1, 0),
                              labels = c("Yes", "No")))

# you might need to tune this based on the range of your figure.
nudge_y <- (as.numeric(data_segment$event_status=='Yes')-1/2)*0.01

fig <- ggplot(spline_preds) +
 aes(x = x,
     y = exp(pred),
     ymin = exp(ci_lwr),
     ymax = exp(ci_upr)) +
 labs(x = 'age, years',
      y = 'Prevalence ratio',
      color = 'Blood pressure control') +
 geom_line() +
 geom_ribbon(alpha = 0.2) +
 scale_y_log10(limits = c(0.35, 2.5)) +
 geom_segment(data = data_segment,
              inherit.aes = FALSE,
              size = 12,
              mapping = aes(x = x,
                            y = y,
                            color = event_status,
                            xend = xend,
                            yend = yend)) +
 geom_text(data = data_segment,
           inherit.aes = FALSE,
           show.legend = FALSE,
           nudge_y = nudge_y,
           mapping = aes(x = x,
                         y = yend,
                         vjust = as.numeric(event_status=='No'),
                         color = event_status,
                         label = count)) +
 theme_bw() +
 geom_hline(yintercept = 1, linetype = 2, color = 'grey') +
 theme(panel.grid = element_blank(),
       legend.position = c(.25, 0.80)) +
 scale_color_manual(values = c("grey", "black"))

