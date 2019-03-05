
stan_mode <- function(fit, par, d=3) {
  x <- round(rstan::extract(fit, par)[[1]], d)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' change time interval of the data (from 15min interval). Not very general, use carefully
change_interval <- function(.df, unit="hour") {
  if (is.numeric(unit)) { # use for upsampling
    .df <- tibble(timestamp=seq.POSIXt(.df$timestamp[1], .df$timestamp[nrow(.df)], unit)) %>%
      left_join(.df, by="timestamp") %>% mutate_if(~all(class(.)=="numeric"), funs(zoo::na.approx(.)))
  } else if (unit=="15min") {
    .df
  } else if (unit=="30min") {
    .df %>% mutate_if(is.numeric, funs(RcppRoll::roll_mean(c(.[1], ., .[n()]), weights = c(0.25, 0.5, 0.25)))) %>%
      filter(format(timestamp, "%M") %in% c("00", "30")) %>% return()
  } else if (unit=="hour") {
    .df %>% mutate_if(is.numeric, funs(RcppRoll::roll_mean(c(.[1], ., .[n()]), 3))) %>%
      filter(format(timestamp, "%M") == "00") %>% return()
  } else {
    .df %>%
      #mutate(ts=lubridate::round_date(timestamp, unit)) %>%
      #group_by(ts) %>% summarise_all(mean) 
      mutate(timestamp=lubridate::round_date(timestamp, unit)) %>%
      group_by(timestamp) %>% summarise_all(mean) 
  }
}

