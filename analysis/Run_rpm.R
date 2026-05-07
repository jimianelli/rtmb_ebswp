source(file.path("R", "config.R"))
fit <- nlminb(obj$par, obj$fn, obj$gr); 
#cest <- wham::check_estimability(obj);cest$BadParams[cest$BadParams$Param_check=="Bad",]
#unique(names(parms))
#sdr <- sdreport(obj)
#rtmb_df <- tibble( name = names(sdr$value), sd   = sdr$sd )
