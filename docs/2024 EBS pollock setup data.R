# =============================================================================
# 2024 EBS pollock assessment data for Rceattle (CEATTLE)
# =============================================================================
# Single-sex, single-species: fishery + AVO acoustic index, BTS bottom-trawl and
# ATS acoustic-trawl surveys, ATS age-1 index, and the 1965-76 Japanese CPUE.
# Aligns Rceattle with the ADMB reference ./ADMB/m23_rceattle_full/ and writes
# Data/2024_EBS_pollock_canonical_pm.xlsx (run once; the corrected workflow reads it).
#
# ADMB bridge -- the "pm" edits this config matches (flagged "MODIFIED (m23_...)"
# in ADMB/*/pm.tpl):
#   Structural: log_avg_F off + plain log_F_devs (one free F/yr); BTS sel-dev
#     vectors plain with year 1 pinned; wt_like excluded; initial-age geometric
#     series = equilibrium + init devs (Rceattle initMode = "NonEquilibrium").
#   Likelihood: rec_like full normal (sigr = 1); Ricker rec penalty + steepness off;
#     ATS biomass index and AVO exclude age-1; log_q_avo bounded [-15,0]; terminal
#     ATS age-1 obs dropped from q and fit.
#   Data/penalty: BTS sel random-walk penalty over the survey period only; the
#     PM stabilized multinomial retains nominal observed-proportion weights and
#     applies MN_const only inside logarithms.
#
# Rceattle-side encoding (in the body below; configuration, not ADMB code edits):
#   fishery terminal-year length comp not fit; AVO obs in million-tonnes (absolute-SD
#   normal); CPUE fit once as a survey fleet mirroring the fishery selectivity;
#   BTS/ATS comp sample sizes truncated to integer, 2020 ATS = 1; ATS index
#   Log_sd = sqrt(log(CV^2+1)); AMAK avgsel penalty on (Sel_avgsel_pen = 10).
#
# Rebuild the reference:
#   cd ADMB/m23_rceattle_full && export PATH=/usr/local/bin:$PATH && admb pm && ./pm -nox -iprint 150
# =============================================================================

.libPaths(c(file.path(getwd(), ".r-lib-pm"), .libPaths()))
library(Rceattle)
library(dplyr)
library(readxl)

n_selages_fsh <- 12; bts_styr <- 1982; ats_styr <- 1994

# -----------------------------------------------------------------------------
# Data ----
# -----------------------------------------------------------------------------
mydata <- Rceattle::read_data(file = "Data/2024_EBS_pollock.xlsx")
styr <- mydata$styr
endyr <- mydata$endyr
nages <- mydata$nages
yrs  <- styr:endyr
nyr <- length(yrs)

keep_age <- c("Species_name", "Species", "Sex", "Year", paste0("Age", 1:nages))
mydata$NByageFixed <- mydata$NByageFixed[, intersect(keep_age, colnames(mydata$NByageFixed))]
mydata$spawn_month <- 3 # ADMB yrfrac 0.25

est <- mydata
est$estDynamics <- 0
fcn <- est$fleet_control$Fleet_name

# -- observation errors --------------------------------------------------------
# NOTE: the xlsx index Log_sd is ALREADY a CV / log-sd (0.05-0.56; ADMB sdnr ~1).
# Do NOT divide it by Observation. catch Log_sd = 0.05 (ADMB ctrl_flag(1)=200 =>
# sigma = 1/sqrt(2*200) = 0.05). ATS age-1 index sigma = age1_sigma_ats = 1.
est$catch_data$Log_sd <- 0.05
est$index_data$Log_sd[est$index_data$Fleet_name %in% c("BTS_1", "ATS_1")] <- 1
est$fleet_control$Fleet_type[fcn %in% c("BTS_1", "ATS_1")] <- "Survey"
est$age_error[1:nages, 3:(nages + 2)] <- diag(nages)       # ageing error off
est$sigma_rec_prior <- 1                                   # full-normal rec penalty (ADMB L1)

# -- selectivity forms (AMAK "pm"): Fishery = Ianelli non-parametric,
#    BTS = logistic + free age-1, ATS/AVO = non-parametric ascending-constrained.
#    Penalty weights come from ADMB ctrl_flags / selvar24.dat.
est$fleet_control$Selectivity[fcn == "Fishery"]               <- "NonParametricPM"
est$fleet_control$Time_varying_sel[fcn == "Fishery"]          <- "RandomWalk"
est$fleet_control$N_sel_bins[fcn == "Fishery"]                <- n_selages_fsh
est$fleet_control$Sel_curve_pen1[fcn == "Fishery"]            <- 12.5    # ctrl_flag(13)
est$fleet_control$Sel_curve_pen2[fcn == "Fishery"]            <- 1/60    # ctrl_flag(11)/nch
est$fleet_control$Sel_curve_pen3                              <- 0
est$fleet_control$Sel_curve_pen3[fcn == "Fishery"]            <- 1       # ctrl_flag(10)/group
est$fleet_control$Sel_norm_bin1[fcn == "Fishery"]             <- NA
est$fleet_control$Time_varying_sel_sd_prior[fcn == "Fishery"] <- 0.5     # selvar24.dat

est$fleet_control$Selectivity[fcn == "BTS"]                <- "LogisticPM"
est$fleet_control$Time_varying_sel[fcn == "BTS"]           <- "RandomWalk"
est$fleet_control$Sel_curve_pen1[fcn == "BTS"]             <- 2          # ctrl_flag(26)
est$fleet_control$Sel_curve_pen2[fcn == "BTS"]             <- 0
est$fleet_control$Sel_curve_pen3[fcn == "BTS"]             <- 8          # age-1-dev RW weight
est$fleet_control$Sel_norm_bin1[fcn == "BTS"]              <- 3          # penalty age-range lo
est$fleet_control$Sel_norm_bin2[fcn == "BTS"]              <- 14         # penalty age-range hi
est$fleet_control$Sel_start_year[fcn == "BTS"]             <- bts_styr
est$fleet_control$Bin_first_selected[fcn == "BTS"]         <- 1
est$fleet_control$Time_varying_sel_sd_prior[fcn == "BTS"]  <- 1

for (fl in c("ATS", "AVO")) {
  est$fleet_control$Selectivity[fcn == fl]               <- "NonParametricPM"
  est$fleet_control$Time_varying_sel[fcn == fl]          <- "RandomWalk"
  est$fleet_control$N_sel_bins[fcn == fl]                <- 8
  est$fleet_control$Sel_curve_pen1[fcn == fl]            <- -1           # penalise INCREASING
  est$fleet_control$Sel_curve_pen2[fcn == fl]            <- 1
  est$fleet_control$Sel_curve_pen3[fcn == fl]            <- 0
  est$fleet_control$Sel_norm_bin1[fcn == fl]             <- NA
  est$fleet_control$Bin_first_selected[fcn == fl]        <- 2            # exclude age-1 (ADMB L4/L5)
  est$fleet_control$Sel_pen_first_bin[fcn == fl]         <- 2            # mina_ats
  est$fleet_control$Sel_start_year[fcn == fl]            <- ats_styr
  est$fleet_control$Time_varying_sel_sd_prior[fcn == fl] <- 0.138        # selvar24.dat
}

# AMAK "avgsel" base-level penalty: fff += 10*square(log(mean(exp(base coffs))))
# (pm.tpl:5535) on the type-9 fleets, accumulated once per shared block (lead fleet).
est$fleet_control$Sel_avgsel_pen[fcn %in% c("Fishery", "ATS", "AVO")] <- 10

# -- survey timing + catchability ---------------------------------------------
est$index_data <- est$index_data %>%
  mutate(Month = case_when(Fleet_name %in% c("BTS", "BTS_1", "ATS", "ATS_1") ~ 6, TRUE ~ 0))
est$comp_data <- est$comp_data %>%
  mutate(Month = case_when(Fleet_name == "BTS" ~ 6, Fleet_name == "ATS" ~ 6, TRUE ~ Month))
est$fleet_control$Catchability <- as.character(est$fleet_control$Catchability)
est$fleet_control$Catchability[fcn == "ATS"]                 <- "Estimated"
est$fleet_control$Catchability[fcn %in% c("BTS_1", "ATS_1")] <- "Analytical"        # geometric-mean
est$fleet_control$Index_loglike[fcn == "BTS"] <- "MVN"                              # DoCovBTS
est$fleet_control$Catchability[fcn == "BTS"]  <- "AnalyticalArith"
# BTS survey biomass variance-covariance matrix (42x42, VAST-derived, 1982-2023).
# It is embedded in the written xlsx (index_cov round-trips), so the fit reads it
# from the xlsx; the source matrix ships with the model in Data/.
est$index_cov <- list(BTS = as.matrix(read.table("Data/BTS_survey_covariance_2024.dat")))

# -- ATS biomass index: the xlsx Log_sd is a CV (std/obs), but ADMB's lognormal
#    variance is lvarb_ats = log(CV^2 + 1) (the exact CV -> log-scale-SD conversion,
#    pm.tpl:1689-1691), and Rceattle's lognormal likelihood uses Log_sd directly as
#    the log-scale SD. Convert CV -> sqrt(log(CV^2 + 1)) so the ATS biomass variance
#    matches ADMB exactly (the +0.01 inside-log offset is negligible at this scale).
ats_rows <- est$index_data$Fleet_name == "ATS"
est$index_data$Log_sd[ats_rows] <- sqrt(log(est$index_data$Log_sd[ats_rows]^2 + 1))

# -- AVO acoustic index: ADMB avo_like is a natural-scale normal with an ABSOLUTE
#    observation SD (ob_avo_std, pm_24.dat), not a lognormal CV. Fit it with
#    Index_loglike = "Normal" (residual obs - q*pred ~ N(0, ob_avo_std^2)) and
#    supply ob_avo_std directly in Log_sd (provided, not estimated).
est$fleet_control$Index_loglike[fcn == "AVO"]     <- "Normal"
# All index SDs are provided (not estimated). Set the WHOLE column as a string alias:
# a per-fleet string assignment would coerce the numeric column to character and leave
# "0" strings on the untouched fleets, which strict (dev-line) validators reject.
est$fleet_control$Estimate_index_sd <- "Fixed"
ob_avo_std <- setNames(
  c(0.407974331, 0.79543824, 0.292865177, 0.390095688, 0.579193251, 0.447677778,
    0.371938445, 0.390115995, 0.58024587, 0.406257388, 0.379092753, 0.317389245,
    0.254960502, 0.63539506, 0.529928784, 0.454780316, 0.335349192, 0.250814465),
  as.character(c(2006:2019, 2021:2024)))
avo_rows <- est$index_data$Fleet_name == "AVO"
stopifnot(sum(avo_rows) == length(ob_avo_std))
est$index_data$Log_sd[avo_rows] <- ob_avo_std[as.character(abs(est$index_data$Year[avo_rows]))]
# AVO obs are thousand-tonnes in the base xlsx (~1741) but ADMB's obs_avo is million-
# tonnes (~1.74); with a natural-scale normal + absolute sigma, rescale to million-
# tonnes so obs, sigma and prediction (q*wt_avo*N*sel_ats) are on the same scale.
est$index_data$Observation[avo_rows] <- est$index_data$Observation[avo_rows] / 1000

# -- composition likelihood: canonical AMAK/PM stabilized multinomial. The
# observed proportions retain nominal multinomial weights, while 0.001 enters
# only inside the logarithms to stabilize empty bins.
est$fleet_control$Comp_loglike <- "MultinomialPM"

# -- ADMB reads survey comp sample sizes as integer vectors (init_ivector sam_bts/
#    sam_ats), truncating the fractional (McAllister-Ianelli) weights; the fishery is a
#    float and left as-is. Truncate BTS/ATS to match so the multinomial weights agree.
for (fl in c("BTS", "ATS"))
  est$comp_data$Sample_size[est$comp_data$Fleet_name == fl] <-
    trunc(est$comp_data$Sample_size[est$comp_data$Fleet_name == fl])
# The 2020 ATS age comp (COVID-year survey) is stored with sample size 0 in the
# xlsx (effectively excluded), but ADMB's data file fits it with sample size 1
# (sam_ats(2020) = 1). Restore it so the ATS multinomial matches ADMB exactly.
est$comp_data$Sample_size[est$comp_data$Fleet_name == "ATS" &
                          est$comp_data$Year == 2020] <- 1

# -- BTS age-1: ADMB keeps age-1 IN the BTS comps and has NO BTS age-1 index; the
#    xlsx relocated it into a separate BTS_1 index (verified identical to the raw
#    comp age-1 count). Restore it to the comps and drop the redundant BTS_1 index.
b1 <- est$index_data[est$index_data$Fleet_name == "BTS_1", c("Year", "Observation")]
for (r in which(est$comp_data$Fleet_name == "BTS")) {
  o <- b1$Observation[abs(b1$Year) == abs(est$comp_data$Year[r])]
  if (length(o) == 1) est$comp_data[r, "Comp_1"] <- o
}
est$fleet_control$Fleet_type[fcn == "BTS_1"] <- "Off"

# -- ATS age-1 (ATS_1): drop the terminal 2024 obs (ADMB ignore_last_ats_age1, last-
#    year CV 1.81) via the negative-year convention (predicted, not fitted, out of the
#    analytical q). Flip -2020 -> 2020 so ATS/ATS_1 fit 2020 (as the ATS comps do).
est$index_data$Year[est$index_data$Fleet_name == "ATS_1" & est$index_data$Year == 2024]  <- -2024
est$index_data$Year[est$index_data$Fleet_name %in% c("ATS", "ATS_1") &
                      est$index_data$Year == -2020] <- 2020

# -- Drop the base xlsx's CPUE copy (attached to the fishery fleet). ADMB fits the CPUE
#    once (cpue_like); re-adding it below as its own survey fleet would double-count the
#    only early-period index and over-constrain the initial age structure.
fishery_code <- est$fleet_control$Fleet_code[fcn == "Fishery"]
est$index_data <- est$index_data[est$index_data$Fleet_code != fishery_code, ]

# -- Japanese fishery CPUE index (1965-1976): the only abundance index before the BTS
#    (1982), so it pins the early numbers-at-age / initial age structure. Added as a
#    survey fleet mirroring the fishery selectivity (pred = wt_fsh*natage*sel_fsh*q_cpue)
#    with an estimated q; cpue_like is a natural-scale normal, absolute SD in Log_sd.
# NOTE: the fleet_control column names below (Q_index, Q_prior, Q_sd_prior,
# Time_varying_q_sd_prior, Index_sd_prior) are the deprecated spellings. They still
# work on the release package but are auto-upgraded to canonical names on the dev
# line (dev-data-workflow / perf-profiling), so this block errors there. TODO: migrate
# to the canonical names (Q_index -> Catchability_index, Q_prior -> Catchability_init,
# Q_sd_prior -> Catchability_prior_sd, Time_varying_q_sd_prior -> Time_varying_q_sd,
# Index_sd_prior -> Index_sd) so data.R re-runs on the dev line too.
cpue_row <- est$fleet_control[fcn == "Fishery", ]          # inherit fishery selectivity
cpue_row$Fleet_name <- "CPUE"
cpue_row$Fleet_code <- max(est$fleet_control$Fleet_code) + 1L
cpue_row$Fleet_type <- "Survey"
cpue_row$Q_index    <- max(est$fleet_control$Q_index, na.rm = TRUE) + 1L
cpue_row$Catchability <- "Estimated"                       # estimated q (log_q_cpue)
cpue_row$Index_loglike     <- "Normal"                     # natural-scale normal, absolute SD
cpue_row$Estimate_index_sd <- "Fixed"
avo_r <- which(fcn == "AVO")
for (col in c("Q_prior", "Q_sd_prior", "Time_varying_q", "Time_varying_q_sd_prior",
              "Estimate_index_sd", "Index_sd_prior"))
  cpue_row[[col]] <- est$fleet_control[avo_r, col]
for (col in c("Estimate_catch_sd", "Catch_sd_prior", "proj_F_prop")) cpue_row[[col]] <- NA
est$fleet_control <- rbind(est$fleet_control, cpue_row)

cpue_obs <- c(2816.437428, 3473.580475, 3802.169891, 5257.304601, 6712.468418,
              5679.809828, 5257.331283, 5726.743484, 4787.923949, 4740.992588,
              4271.574460, 4318.523058)
cpue_sd  <- c(563.2874856, 694.716095, 760.4339781, 1051.46092, 1342.493684,
              1135.961966, 1051.466257, 1145.348697, 957.5847898, 948.1985176,
              854.3148919, 863.7046116)
idx0 <- est$index_data[est$index_data$Fleet_name == "AVO", ][1, ]
cpue_idx <- idx0[rep(1, length(cpue_obs)), ]
cpue_idx$Fleet_name  <- "CPUE"
cpue_idx$Fleet_code  <- cpue_row$Fleet_code
cpue_idx$Species     <- 1
cpue_idx$Year        <- 1965:1976
cpue_idx$Month       <- 0
cpue_idx$Observation <- cpue_obs
cpue_idx$Log_sd      <- cpue_sd                            # absolute observation SD (natural-scale normal)
est$index_data <- rbind(est$index_data, cpue_idx)

# -- Fishery length comp: NOT fit. ADMB use_endyr_len = 0 excludes the terminal-year
#    length comp from the objective (pm.rep len_like = 25.795 reported but not summed),
#    so we omit it too -- the terminal fish are already in the fishery AGE comp.

write_data(est, file = "Data/2024_EBS_pollock_canonical_pm.xlsx")
