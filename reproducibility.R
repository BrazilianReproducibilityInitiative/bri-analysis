### ### ### ### ### ### ### ### ### ### ### ###

# reproducibility.r
# This script contains the functions to run the meta-analyses of the replications.

### ### ### ### ### ### ### ### ### ### ### ###

# Paths to necessary files

# Data from experiments
data_paths = c(
  "./replication-results/MTT",
  "./replication-results/PCR",
  "./replication-results/EPM/primary_outcomes"
)

# Metadata about experiments
all_data_MTT = fread(file = "./gathered-data/MTT_data.csv", dec = ",", na.strings = c("","NA"), encoding = "Latin-1") |> as_tibble() |> filter(!is.na(Group))

all_data_PCR = fread(file = "./gathered-data/PCR_data.csv", dec = ",", na.strings = c("","NA", "Undetermined", "undetermined"), encoding = "Latin-1") |> as_tibble() |> filter(!is.na(Group))

all_data_EPM = read_csv2("./gathered-data/EPM_data.csv", show_col_types = F)

meta_data_MTT = read_csv2("./gathered-data/MTT_meta.csv", show_col_types = F)

meta_data_PCR = read_csv2("./gathered-data/PCR_meta.csv", show_col_types = F)

meta_data_EPM = read_csv2("./gathered-data/EPM_meta.csv", show_col_types = F)

# Subjective reproduciblity data
subjective_repro_data_path = "./other-data/Subjective assessment of reproducibility.xlsx"
subjective_repro_data = read_excel(subjective_repro_data_path)

# Original experiments data
original_experiment_data_path = "./other-data/Original Experiments Statistical Summaries.tsv"
original_experiment_data = read_tsv(original_experiment_data_path, show_col_types = F)

# Dictionary mapping experiments to data files
file_dict = rbind(
  read_tsv("./replication-results/PCR_results_dict.tsv", show_col_types = F),
  read_tsv("./replication-results/PCR_ALT_results_dict.tsv", show_col_types = F),
  read_tsv("./replication-results/EPM_results_dict.tsv", show_col_types = F),
  read_tsv("./replication-results/MTT_results_dict.tsv", show_col_types = F)
) |>
  filter(
    !is.na(Filename),
    # Excluding secondary outcomes from EPM experiments
    !str_detect(Filename, "secondary")
  )


# Main function: given an inclusion set, it runs all the meta-analyses
run_all_meta_analyses = function (inclusion_set_column, save_results_to, params, simulated = F) {
  # Creates path to save
  path_suffix = paste0(
    " ", params$ma_dist,
    ifelse(params$exclude_outliers, " no outliers", "")
  )
  
  dir.create(paste0("./output/", save_results_to))
  output_path = paste0("./output/", save_results_to, "/", inclusion_set_column |> str_remove("analysis_"), path_suffix)
  dir.create(output_path)
  
  # Builds the list of meta-analyses to run from inclusion set
  data_list = inclusion_sets
  data_list$toinclude = data_list[[inclusion_set_column]]
  data_list = data_list |> filter(toinclude %in% c("INCLUDE", "MA ONLY"))
  ma_only_list = data_list |> filter(toinclude == "MA ONLY") |> select(EXP, LAB)
  data_list = data_list |> select(LAB, EXP, UNIT) |>
    mutate(UNIT = ifelse(UNIT == "BRI", NA, UNIT)) |>
    left_join(file_dict, by = c("LAB","EXP", "UNIT")) |>
    mutate(ALT = F)
  
  # Adds alternative PCR measures
  pcr_alt = data_list |> filter(str_detect(EXP, "PCR")) |>
    mutate(
      Filename = Filename |> str_replace("PCR/Result", "PCR_ALT/Result") |> str_replace(" PCR", " ALTPCR"),
      EXP = EXP |> str_replace("PCR", "ALTPCR")
    )
  ma_only_list_alt = ma_only_list |> filter(str_detect(EXP, "PCR")) |>
    mutate(
      EXP = EXP |> str_replace("PCR", "ALTPCR")
    )
  
  if (!simulated) {
    data_list = rbind(data_list, pcr_alt)
    ma_only_list = rbind(ma_only_list, ma_only_list_alt)
  }
  
  # data_list = data_list |> filter(str_detect(EXP, "MTT56"))
  exps = unique(data_list$EXP)
  
  # Run each analysis 
  analysis_results = map(exps, perform_analysis, data_list, output_path, params, ma_only_list, simulated = simulated)
  
  # Extracts forest plots and data
  forest_plots = map(analysis_results, "plot")
  rep_summaries = map_dfr(analysis_results, "data")
  rep_summaries_individual = map_dfr(analysis_results, "individual_data")
  
  if (simulated) {
    if (nrow(rep_summaries_individual) < 100) return (NA)
    
    sim_alerts = unlist(map(analysis_results, "sim_failed"))
    
    # Breaks if there were any alerts from the data generation for the simulations
    if (any(sim_alerts == T)) {
      return (NA)
    }
    
    return (list(by_exp = rep_summaries, by_rep = rep_summaries_individual))
  }
  
  # Writes output of the analyses, by experiment and by replication
  write_tsv(rep_summaries, paste0(output_path, "/Replication Assessment by Experiment.tsv"))
  write_tsv(rep_summaries_individual, paste0(output_path, "/Replication Assessment by Replication.tsv"))
}

# Function to perform the meta-analysis of each dataset and plot the effects
perform_analysis = function (EXP_code, df, output_path, params, ma_only_list, simulated = F) {
  print(paste0(EXP_code, " ", params$ma_dist))
  
  # Flags for PCR types and exponents
  is_PCR = str_detect(EXP_code, "PCR") & !str_detect(EXP_code, "ALTPCR")
  is_ALTPCR = str_detect(EXP_code, "ALTPCR")
  
  if (is_PCR | is_ALTPCR) {
    PCR_exponent = ifelse((meta_data_PCR |> filter(EXP == EXP_code) |> pull(`Primary Outcome`))[1] == "Optical Density (a.u)", 1, -1)
  } else {
    PCR_exponent = 0
  }
  
  # Read experiment data
  data_filenames = df |> filter(EXP == EXP_code & !is.na(Filename)) |> pull(Filename)
  
  # Paired experiment data
  uses_lab_unit = any(data_filenames |> str_detect("LABUNIT"))
  if (uses_lab_unit) {
    paired_labs = rbind(
      meta_data_EPM |> filter(EXP == EXP_code) |> select(LAB, Paired),
      meta_data_PCR |> filter(EXP == EXP_code) |> select(LAB, `Paired Alternative`) |> rename(Paired = `Paired Alternative`),
      meta_data_MTT |> filter(EXP == EXP_code) |> select(LAB, Paired)
    ) |> filter(Paired %in% c("yes","Yes")) |> pull(LAB)
  } else {
    paired_labs = rbind(
      meta_data_EPM |> filter(EXP == EXP_code) |> select(LAB, Paired),
      meta_data_PCR |> filter(EXP == EXP_code) |> select(LAB, Paired),
      meta_data_MTT |> filter(EXP == EXP_code) |> select(LAB, Paired)
    ) |> filter(Paired %in% c("yes","Yes")) |> pull(LAB)
  }
  
  if (is_ALTPCR) {
    data_filenames_PCR_ref = df |> filter(EXP == EXP_code |> str_remove("ALT") & !is.na(Filename)) |> pull(Filename)
  }
  
  orig_data = original_experiment_data |> filter(EXP == str_remove(EXP_code, "ALT"))
  
  # Get original effect size
  original_n = orig_data$`Control Sample Size for DFs` + orig_data$`Treated Sample Size for DFs`
  orig_summary = escalc(measure = "ROM", m1i = orig_data$`Treated Mean`, m2i = orig_data$`Control Mean`, sd1i = orig_data$`Treated SD`, sd2i = orig_data$`Control SD`, n1i = orig_data$`Treated Sample Size for DFs`, n2i = orig_data$`Control Sample Size for DFs`)
  
  ## Calculate CIs with degrees of freedom from the original
  if (params$ma_dist == "z") {
    orig_crit = qnorm(0.025, lower.tail = FALSE)
  } else if (params$ma_dist == "t" | params$ma_dist == "knha") {
    orig_crit = qt(0.025, df = orig_data$degrees_of_freedom, lower.tail = FALSE)
  }
  
  original_es = orig_summary$yi
  original_sem = sqrt(orig_summary$vi)
  original_ci_lower = original_es - orig_crit * original_sem
  original_ci_upper = original_es + orig_crit * original_sem
  
  # Analyze replications and extract summaries
  all_rep_es_data = make_rep_es_analysis(data_filenames, simulated, EXP_code, is_PCR, original_es, paired_labs)
  
  if (is_ALTPCR) {
    all_rep_es_data_PCR_ref = make_rep_es_analysis(data_filenames_PCR_ref, simulated, EXP_code, T, original_es, paired_labs)
  }
  
  summaries = all_rep_es_data$summaries
  rep_es = all_rep_es_data$rep_es
  rep_data = all_rep_es_data$rep_data
  rep_enough_ns = all_rep_es_data$rep_enough_ns
  
  if (simulated) {
    if (nrow(summaries) != nrow(rep_es)) {
      return (list(data = tibble(fema_pvalue = NA, indiv_pvalue = NA), plot = ggplot(), individual_data = tibble(), sim_failed = T))
    }
  }
  
  if (nrow(summaries) == 0)
    return (list(data = tibble(fema_pvalue = NA, indiv_pvalue = NA), plot = ggplot(), individual_data = tibble()))
  
  dir.create(paste0(output_path, "/escalc"))
  fn = paste0(output_path, "/escalc/",
              "escalc - ", EXP_code, ".tsv")
  
  # Get coefficient of variation for original and experiments with the summaries using the _Perc columns
  
  rep_es_with_perc = make_rep_es_analysis(data_filenames, simulated, EXP_code, is_PCR, original_es, paired_labs, use_perc = T)
  
  original_cv = orig_data$`Pooled SD` / ((orig_data$`Treated Mean` + orig_data$`Control Mean`) / 2)
  
  rep_cv = rep_es_with_perc$rep_es |>
    mutate(
      pooled_sd = sqrt((sd_group_1 ^ 2 + sd_group_2 ^ 2) / 2),
      cv = pooled_sd / ((mean_group_1 + mean_group_2) / 2)
    ) |>
    pull(cv)
  
  mean_rep_cv = mean(rep_cv, na.rm = T)
  
  if (!simulated) {
    escalc_summary = summary(summaries) |>
      select(yi, ci.lb, ci.ub)
    class(escalc_summary) = "data.frame"
    escalc_summary$LAB = rep_es$LAB
    write_tsv(escalc_summary, file = fn, quote = "all")
  }
  
  summaries_label = all_rep_es_data$summaries
  
  if (is_PCR) {
    es_labels = summary(summaries_label) |>
      mutate(
        s = paste0(
          # "2 ^ Difference between means [95% CI]: ",
          "Log transformed ratio of means [95% CI]: ",
          round(exp_only_PCR(yi, T, PCR_exponent), 2),
          " [", round(exp_only_PCR(ci.lb, T, PCR_exponent), 2), ", ",
          round(exp_only_PCR(ci.ub, T, PCR_exponent), 2), "]"
        )
      ) |> pull(s)
  } else {
    es_labels = summary(summaries_label) |>
      mutate(
        s = paste0(
          "Log transformed ratio of means [95% CI]: ",
          round(yi, 2),
          " [", round(ci.lb, 2), ", ",
          round(ci.ub, 2), "]"
        )
      ) |> pull(s)
  }
  
  # Make graphs for individual experiments
  if (!simulated) {
    print("Making individual graphs ...")
    walk2(data_filenames[str_detect(data_filenames, pattern = paste0(rep_es$LAB, collapse = "|"))], es_labels, plot_individual_results, output_path)
  }
  
  # Write effect size summary to file for all replications
  if (!simulated) {
    dir.create(paste0(output_path, "/summaries"))
    fn = paste0(output_path, "/summaries/",
                "Summary per EXP - ", EXP_code, ".tsv")
    
    write_tsv(file = fn, rep_es, quote = "all")
  }
  
  n_reps = rep_data |> pull(LAB) |> unique() |> length()
  if (n_reps == 0) {
    return (0)
  }
  
  print(paste(rep_data |> pull(LAB) |> unique(), sep = ", "))
  
  
  zval_FEMA = NA
  
  # If there is more than one replication, run meta-analysis
  if (n_reps > 1) {
    # Perform random and fixed effects meta-analysis
    print ("Running meta-analyses and individual tests ...")
    
    ma_welch_df = rep_es$welch_df[1]
    # write_tsv(rep_es, file = paste0("_", EXP_code, "WELCH.tsv"))
    
    if (!simulated) {
      
      if (params$ma_dist == "t") {
        RE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "REML", test = "t", dfs = ma_welch_df, control = list(stepadj = 0.25, maxiter = 500000))
      } else if (params$ma_dist == "z") {
        RE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "REML", test = "z", control = list(stepadj = 0.25, maxiter = 500000))
      } else if (params$ma_dist == "knha") {
        RE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "REML", test = "knha", control = list(stepadj = 0.25, maxiter = 500000))
      }
      
      pred = predict(RE_MA)
      
      rema_es = RE_MA$beta[1,1]
      rema_pi_lower = pred$pi.lb
      rema_pi_upper = pred$pi.ub
      rema_pvalue = RE_MA$pval
      
      rema_I2 = RE_MA$I2
      rema_tau2 = RE_MA$tau2
      rema_Qp = RE_MA$QEp
      
    } else {
      
      rema_es = NA
      rema_pi_lower = NA
      rema_pi_upper = NA
      rema_pvalue = NA
      
      rema_I2 = NA
      rema_tau2 = NA
      rema_Qp = NA
      
    }
    
    t_FE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "FE", test = "t", dfs = ma_welch_df)
    zval_FEMA = t_FE_MA$zval
    
    if (params$ma_dist == "t") {
      FE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "FE", test = "t", dfs = ma_welch_df)
    } else if (params$ma_dist == "z") {
      FE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "FE", test = "z")
    } else if (params$ma_dist == "knha") {
      FE_MA = rma(yi = summaries$yi, vi = summaries$vi, method = "FE", test = "knha")
    }
    
    pred = predict(FE_MA)
    
    fema_es = FE_MA$beta[1,1]
    fema_ci_lower = pred$ci.lb
    fema_ci_upper = pred$ci.ub
    fema_pvalue = FE_MA$pval
  }
  
  # If only one replication is available, no meta-analyses
  if (n_reps == 1) {
    print ("Not enough replications for meta-analysis, running single test ...")
    
    # The rest is defined as NA, will be used to fill the results
    rema_es = NA
    rema_pi_lower = NA
    rema_pi_upper = NA
    rema_pvalue = NA
    
    rema_I2 = NA
    rema_tau2 = NA
    rema_Qp = NA
    
    fema_es = NA
    fema_ci_lower = NA
    fema_ci_upper = NA
    fema_pvalue = NA
    
  }
  
  # Store results and bounds for each experiment individually
  summaries = summary(summaries)
  
  summaries = summaries |> mutate(
    dfi = rep_es$dfi
  )
  
  # Calculate individual experiment intervals
  if (params$ma_dist == "z") {
    summaries = summaries |> mutate(
      crit = qnorm(0.025, lower.tail = FALSE),
      pval = 2 * pnorm(abs(zi), lower.tail = FALSE)
    )
  } else if (params$ma_dist == "t" | params$ma_dist == "knha") {
    summaries = summaries |> mutate(
      crit = qt(0.025, df = dfi, lower.tail = FALSE),
      pval = 2 * pt(abs(zi), df = dfi, lower.tail = FALSE)
    )
  }
  
  summaries = summaries |> mutate(
    ci.lb = yi - crit * sei,
    ci.ub = yi + crit * sei
  )
  
  # Transform PCR estimates
  summaries = summaries |>
    mutate(
      yi = exp_only_PCR(yi, is_PCR, PCR_exponent),
      ci.lb = exp_only_PCR(ci.lb, is_PCR, PCR_exponent),
      ci.ub = exp_only_PCR(ci.ub, is_PCR, PCR_exponent)
    )
  class(summaries) = "data.frame"
  
  summaries$LAB = rep_es$LAB
  
  indiv_estimate = summaries$yi
  indiv_ci_lower = summaries$ci.lb
  indiv_ci_upper = summaries$ci.ub
  indiv_pvalue = summaries$pval
  
  if (!simulated) {
    
    # Make a forest plot including the original effect size
    print ("Making forest plot ...")
    
    original = tibble(
      yi = original_es,
      lb = original_ci_lower,
      ub = original_ci_upper,
      LAB = "Original Effect",
      weight = original_n
    )
    
    psummaries = summaries |>
      select(LAB, yi, ci.lb, ci.ub) |>
      rename(lb = ci.lb, ub = ci.ub) |>
      left_join(rep_enough_ns, by = "LAB") |>
      mutate(weight = n) |>
      select(-n)
    
    summaries = summaries |>
      rename(lb = ci.lb, ub = ci.ub) |>
      left_join(rep_enough_ns, by = "LAB") |>
      mutate(weight = n) |>
      select(-n)
    
    if (exists("RE_MA")) {
      
      ma_weight = (summaries |> pull(weight) |> sum())
      
      ma_pred = predict(RE_MA)
      
      ma_summary_re_ci = tibble(
        yi = ma_pred$pred,
        lb = ma_pred$ci.lb,
        ub = ma_pred$ci.ub,
        LAB = "RE Model CI",
        weight = ma_weight
      ) |> mutate(
        yi = exp_only_PCR(yi, is_PCR, PCR_exponent),
        lb = exp_only_PCR(lb, is_PCR, PCR_exponent),
        ub = exp_only_PCR(ub, is_PCR, PCR_exponent)
      )
      
      ma_summary_re_pi = tibble(
        yi = ma_pred$pred,
        lb = ma_pred$pi.lb,
        ub = ma_pred$pi.ub,
        LAB = "RE Model PI",
        weight = ma_weight
      ) |> mutate(
        yi = exp_only_PCR(yi, is_PCR, PCR_exponent),
        lb = exp_only_PCR(lb, is_PCR, PCR_exponent),
        ub = exp_only_PCR(ub, is_PCR, PCR_exponent)
      )
      
      ma_pred = predict(FE_MA)
      
      ma_summary_fe_ci = tibble(
        yi = ma_pred$pred,
        lb = ma_pred$ci.lb,
        ub = ma_pred$ci.ub,
        LAB = "FE Model CI",
        weight = ma_weight
      ) |> mutate(
        yi = exp_only_PCR(yi, is_PCR, PCR_exponent),
        lb = exp_only_PCR(lb, is_PCR, PCR_exponent),
        ub = exp_only_PCR(ub, is_PCR, PCR_exponent)
      )
      
      plot_data = rbind(original, ma_summary_re_pi, ma_summary_re_ci, ma_summary_fe_ci, psummaries)
    } else {
      plot_data = rbind(original, psummaries)
    }
    
    plot_data = switch_interval_bounds(plot_data)
    
    # dir.create(paste0(output_path, "/plot_data"))
    # fn = paste0(output_path, "/plot_data/",
    #             "plot data - ", ".tsv")
    # write_tsv(plot_data, file = fn, quote = "all")
    
    plot_label = paste0("I2: ", round(rema_I2,2), ", tau2: ", round(rema_tau2, 2), ", Q test p-value: ", round(rema_Qp, 2))
    
    p = forest_plot(plot_data, EXP_code, plot_label, output_path, scale_to_use = "ROM")
  }
  
  # Exponentiate the PCR estimates
  
  if (!simulated) {
    if (exists("RE_MA")) {
      rema_es = exp_only_PCR(rema_es, is_PCR, PCR_exponent)
      rema_pi_lower = exp_only_PCR(rema_pi_lower, is_PCR, PCR_exponent)
      rema_pi_upper = exp_only_PCR(rema_pi_upper, is_PCR, PCR_exponent)
      
      fema_es = exp_only_PCR(fema_es, is_PCR, PCR_exponent)
      fema_ci_lower = exp_only_PCR(fema_ci_lower, is_PCR, PCR_exponent)
      fema_ci_upper = exp_only_PCR(fema_ci_upper, is_PCR, PCR_exponent)
    }
  } else {
    if (exists("FE_MA")) {
      fema_es = exp_only_PCR(fema_es, is_PCR, PCR_exponent)
      fema_ci_lower = exp_only_PCR(fema_ci_lower, is_PCR, PCR_exponent)
      fema_ci_upper = exp_only_PCR(fema_ci_upper, is_PCR, PCR_exponent)
    }
  }
  
  if (n_reps == 1) {
    replication_es = indiv_estimate
  } else {
    if (!simulated) {
      replication_es = rema_es
    } else {
      replication_es = fema_es
    }
  }
  
  if (!simulated) {
    # Perform separate t-tests for each replication, this is only used for the voting criterion
    
    t_tests = summaries |>
      mutate(
        rep_es = indiv_estimate,
        signif = (lb > 0 & ub > 0) | (lb < 0 & ub < 0),
        signif_same_sense = signif & (rep_es / original_es > 0),
        t_statistic = yi / (vi ^ 0.5),
        t_statistic = sign(rep_es / original_es) * abs(t_statistic)
      )
    
    # If there is any replication with n > 1 that is not on the meta analysis only list, then run t test voting on those
    excl = ma_only_list |> filter(EXP == EXP_code, LAB %in% unique(rep_data$LAB)) |> pull(LAB)
    used_tests = t_tests |>
      filter(!(LAB %in% excl))
    n_signif_tests = used_tests |>
      pull(signif_same_sense) |> sum()
    total_tests = nrow(used_tests)
    
    t_tests[t_tests$LAB %in% excl, c("pvalue", "t_statistic")] = NA
    
    # Count the subjective replication votes (independent of the number of replications)
    subjective_evals = subjective_repro_data |>
      filter(`Código do Experimento` == (str_remove(EXP_code, "ALT"))) |>
      rename(LAB = `Código do Laboratório`) |>
      mutate(SUBJ_SUCCESS = `Na opinião do laboratório, os resultados do experimento original foram replicados com sucesso?` == "Sim") |>
      select(LAB, SUBJ_SUCCESS)
    
    t_tests = t_tests |> left_join(subjective_evals, by = "LAB")
    used_tests = used_tests |> left_join(subjective_evals, by = "LAB")
    
    n_subjective_successful_eval = nrow(
      used_tests |>
        filter(SUBJ_SUCCESS)
    )
    
  }
  
  # Switch the bounds of intervals to have upper bounds be the larger ones (for PCRs this might be switched because of the negative exponent)
  if (exists("RE_MA")) {
    if (rema_pi_lower > rema_pi_upper) { swap(rema_pi_lower, rema_pi_upper) }
    if (fema_ci_lower > fema_ci_upper) { swap(fema_ci_lower, fema_ci_upper) }
  } else if (simulated) {
    if (exists("FE_MA")) {
      if (fema_ci_lower > fema_ci_upper) { swap(fema_ci_lower, fema_ci_upper) }
    }
  }
  
  if (all(indiv_ci_lower > indiv_ci_upper)) { swap(indiv_ci_lower, indiv_ci_upper) }
  
  # Adjust the sign of the effect sizes so that the ES ratio makes sense (corrected signs (cs_))
  cs_original_es = abs(original_es) # Original is set to be always positive
  cs_replication_es = sign(original_es / replication_es) * abs(replication_es) # Replication is set to be positive if it has the same sense as the original
  cs_replication_individual_es = sign(original_es / indiv_estimate) * abs(indiv_estimate) # Replication is set to be positive if it has the same sense as the original, now using the individual replications
  
  # Calculate absolute differences between effect sizes
  # Replications
  if (length(indiv_estimate) > 1) {
    if (length(indiv_estimate) == 2) {
      rep_abs_diffs = c(
        indiv_estimate[1] - indiv_estimate[2]
      )
      rep_orig_abs_diffs = c(
        indiv_estimate[1] - original_es,
        indiv_estimate[2] - original_es
      )
    } else if (length(indiv_estimate) == 3) {
      rep_abs_diffs = c(
        indiv_estimate[1] - indiv_estimate[2],
        indiv_estimate[2] - indiv_estimate[3],
        indiv_estimate[3] - indiv_estimate[1]
      )
      rep_orig_abs_diffs = c(
        indiv_estimate[1] - original_es,
        indiv_estimate[2] - original_es,
        indiv_estimate[3] - original_es
      )
    }
    mean_abs_diff_reps = mean(abs(rep_abs_diffs))
    mean_abs_diff_reps_orig = mean(abs(rep_orig_abs_diffs))
    
  } else {
    mean_abs_diff_reps = NA
    mean_abs_diff_reps_orig = NA
  }
  
  # Build summary data frame to return, containing everything that will be needed to evaluate replication success
  if (!simulated) {
    rep_summary = tibble(
      EXP = EXP_code,
      original_es = original_es,
      original_ci_lower = original_ci_lower,
      original_ci_upper = original_ci_upper,
      n_completed_reps = n_reps,
      replication_es = replication_es,
      rema_es = rema_es,
      rema_pi_lower = rema_pi_lower,
      rema_pi_upper = rema_pi_upper,
      rema_pvalue = rema_pvalue,
      rema_I2 = rema_I2,
      rema_tau2 = rema_tau2,
      rema_Qp = rema_Qp,
      fema_es = fema_es,
      fema_ci_lower = fema_ci_lower,
      fema_ci_upper = fema_ci_upper,
      fema_pvalue = fema_pvalue,
      t_score_FEMA = abs(zval_FEMA) * sign(fema_es / original_es),
      t_score_mean = mean(t_tests$t_statistic, na.rm = T),
      t_score = ifelse(is.na(fema_es), t_score_mean, t_score_FEMA),
      corrected_sign_original_es = cs_original_es,
      corrected_sign_replication_es = cs_replication_es,
      log_es_ratio = cs_original_es - cs_replication_es,
      es_ratio = exp(log_es_ratio),
      original_cv = original_cv,
      replication_cv = mean_rep_cv,
      mean_abs_diff_reps = mean_abs_diff_reps,
      mean_abs_diff_reps_orig = mean_abs_diff_reps_orig,
      significant_reps = n_signif_tests,
      n_voting_reps = total_tests,
      subjective_success_reps = n_subjective_successful_eval
    )
    
    indiv_summary = tibble(
      LAB = rep_es |> pull(LAB),
      EXP = EXP_code,
      original_es = original_es,
      original_ci_lower = original_ci_lower,
      original_ci_upper = original_ci_upper,
      replication_es = indiv_estimate,
      ci_lower = indiv_ci_lower,
      ci_upper = indiv_ci_upper,
      log_es_ratio_individual = cs_original_es - cs_replication_individual_es,
      es_ratio_individual = exp(log_es_ratio_individual),
      pvalue = indiv_pvalue,
      t_score_individual = t_tests$t_statistic,
      t_pvalue = t_tests$pvalue,
      original_cv = original_cv,
      replication_cv = rep_cv,
      subjective_success = t_tests$SUBJ_SUCCESS,
      ma_only = LAB %in% excl
    )
  } else {
    indiv_summary = tibble(
      LAB = rep_es |> pull(LAB),
      EXP = EXP_code,
      original_es = original_es,
      replication_es = indiv_estimate,
      pvalue = indiv_pvalue
    )
    
    if (length(indiv_pvalue) > 1) {
      indiv_pvalue = NA
    }
    
    rep_summary = tibble(
      EXP = EXP_code,
      original_es = original_es,
      replication_es = replication_es,
      fema_pvalue = fema_pvalue,
      pvalue = indiv_pvalue
    )
    
    p = NA
  }
  
  list(data = rep_summary, plot = p, individual_data = indiv_summary)
}


# Helper function to make summaries for escalc
summarise_es_by_lab = function(df, paired_labs, use_perc = F) {
  
  if (use_perc) {
    df = df |> group_by(LAB) |>
      group_modify(function(.x, .y) {
        if (.y$LAB %in% paired_labs) {
          df = .x |> filter(!is.na(Group1_Perc) & !is.na(Group2_Perc))
        } else {
          df = .x
        }
        
        df
      }) |>
      summarise(
        mean_group_1 = mean(Group1_Perc, na.rm = T),
        mean_group_2 = mean(Group2_Perc, na.rm = T),
        sd_group_1 = sd(Group1_Perc, na.rm = T),
        sd_group_2 = sd(Group2_Perc, na.rm = T),
        n_group_1 = sum(!is.na(Group1_Perc)),
        n_group_2 = sum(!is.na(Group2_Perc))
      )
  } else {
    df = df |> group_by(LAB) |>
      group_modify(function(.x, .y) {
        if (.y$LAB %in% paired_labs) {
          df = .x |> filter(!is.na(Group1) & !is.na(Group2))
        } else {
          df = .x
        }
        
        df
      }) |>
      summarise(
        mean_group_1 = mean(Group1, na.rm = T),
        mean_group_2 = mean(Group2, na.rm = T),
        sd_group_1 = sd(Group1, na.rm = T),
        sd_group_2 = sd(Group2, na.rm = T),
        n_group_1 = sum(!is.na(Group1)),
        n_group_2 = sum(!is.na(Group2))
      )
  }
  
  df
}

# Function to run escalc while handling PCR vs MTT/EPM as well as paired and non-paired experiments (all of which have different effect size measures)
make_rep_es_analysis = function (data_fns, simulated, EXP_code, is_PCR, original_es, paired_labs, use_perc = F) {
  
  rep_data = map_dfr(data_fns, read_tsv, show_col_types = F)
  
  rep_enough_ns = rep_data |>
    count(LAB) |>
    filter(n > 1)
  
  # Exclude data from replications with N == 1
  rep_data = rep_data |> filter(LAB %in% rep_enough_ns$LAB)
  
  # Calculate correlations for paired experiments' effect sizes (before simulations - the empirical correlation is always used on escalc)
  if (use_perc) {
    rcors_for_escalc = rep_data |> group_by(LAB) |>
      summarise(
        ri = cor(Group1_Perc, Group2_Perc, use = "pairwise.complete.obs", method = "pearson")
      )
  } else {
    rcors_for_escalc = rep_data |> group_by(LAB) |>
      summarise(
        ri = cor(Group1, Group2, use = "pairwise.complete.obs", method = "pearson")
      )
  }
  
  # If simulating results for power analyses, generate data here based on the real distribution
  if (simulated) {
    regen = T
    while (regen) {
      rep_data = generate_sim_data_from_replication(rep_data, original_es, EXP_code, paired_labs)
      
      if (nrow(rep_data) == 0) {
        regen = T
        next
      }
      
      # Summarize effect sizes
      rep_es = rep_data |> summarise_es_by_lab(paired_labs, use_perc)
      
      # If any mean is NaN, regenerate the data until you get a valid dataset
      regen = any(
        c(
          is.nan(c(rep_es$mean_group_1, rep_es$mean_group_2)),
          c(rep_es$mean_group_1, rep_es$mean_group_2) == 0
        )
      )
    }
  } else {
    # Summarize effect sizes
    rep_es = rep_data |> summarise_es_by_lab(paired_labs, use_perc)
  }
  
  
  # Different rma measure for paired and non-paired experiments
  rep_es_unpaired = rep_es |> filter(!LAB %in% paired_labs)
  rep_es_paired = rep_es |> filter(LAB %in% paired_labs)
  
  run_escalc_on_df = function (a_rep_es, m, rcor = 0) {
    if (is.na(rcor)) rcor = 0
    escalc(
      measure = m,
      m1i = a_rep_es$mean_group_2,
      m2i = a_rep_es$mean_group_1,
      sd1i = a_rep_es$sd_group_2,
      sd2i = a_rep_es$sd_group_1,
      n1i = a_rep_es$n_group_2,
      n2i = a_rep_es$n_group_1,
      ni = a_rep_es$n_group_1,
      ri = rcor
    )
  }
  
  # Calculate effect size summaries for meta analysis
  s1 = tibble()
  s2 = tibble()
  
  if (is_PCR) {
    if (nrow(rep_es_paired) > 0) {
      s1 = rep_es_paired |> group_by(LAB) |>
        group_map(function (.x, .y) {
          group_cor = rcors_for_escalc |> filter(LAB == .y$LAB[1])
          run_escalc_on_df(.x, "MC", group_cor$ri)  
        })
    }
    if (nrow(rep_es_unpaired) > 0) {
      s2 = rep_es_unpaired |> group_by(LAB) |>
        group_map(function (.x, .y) {
          run_escalc_on_df(.x, "MD")
        })
    }
  } else {
    if (nrow(rep_es_paired) > 0) {
      s1 = rep_es_paired |> group_by(LAB) |>
        group_map(function (.x, .y) {
          group_cor = rcors_for_escalc |> filter(LAB == .y$LAB[1])
          run_escalc_on_df(.x, "ROMC", group_cor$ri)  
        })
    }
    if (nrow(rep_es_unpaired) > 0) {
      s2 = rep_es_unpaired |> group_by(LAB) |>
        group_map(function (.x, .y) {
          run_escalc_on_df(.x, "ROM")
        })
    }
  }
  
  summaries = do.call(rbind, c(s1,s2))
  summaries = summaries |> filter(!is.na(yi) & !is.na(vi))
  
  # if (nrow(summaries) != nrow(rep_es)) browser()
  
  # Calculate Welch's degrees of freedom
  rep_es = rep_es |>
    mutate(
      vi = summaries$vi,
      wi = 1 / vi,
      ai = wi / sum(wi),
      var_total = sum(ai ^ 2 * vi),
      
      welch_num_i = ((sd_group_1 ^ 2) / n_group_1 + (sd_group_2 ^ 2) / n_group_2) ^ 2,
      welch_den_i = ((sd_group_1 ^ 4) / ((n_group_1 - 1) * n_group_1 ^ 2)) + ((sd_group_2 ^ 4) / ((n_group_2 - 1) * n_group_2 ^ 2)),
      dfi = welch_num_i / welch_den_i,
      
      welch_num = var_total ^ 2,
      welch_den = sum((ai ^ 2 * vi) ^ 2 / dfi),
      welch_df = welch_num / welch_den
    ) |>
    select(-welch_num,-welch_den)
  
  return (
    list(
      summaries = summaries,
      rep_es = rep_es,
      rep_data = rep_data,
      rep_enough_ns = rep_enough_ns
    )
  )
}

# Generates data for experiments based on the original effect sizes, for the post hoc power simulations
# Written to generate data in the same format, so that it can be embedded in the rest of the scripts replacing the actual data
generate_sim_data_from_replication = function (empirical_data, original_es, exp_code, paired_labs) {
  is_PCR = str_detect(exp_code, "PCR")
  
  empirical_es = empirical_data |> summarise_es_by_lab(paired_labs)
  
  # Generates simulated datasets based on the summaries of the empirical results
  all_labs_simdata = pmap_dfr(empirical_es, function (LAB, mean_group_1, mean_group_2, sd_group_1, sd_group_2, n_group_1, n_group_2) {
    
    max_n = max(c(n_group_1, n_group_2), na.rm = T)
    
    if (any(is.nan(c(mean_group_1, mean_group_2, sd_group_1, sd_group_2, n_group_1, n_group_2)))) return (tibble())
    
    if (is_PCR) {
      sim_data = tibble(
        LAB = LAB,
        EXP = exp_code,
        Replicate = 1:max_n,
        Group1 = c(
          rnorm(n = n_group_1, mean = mean_group_1, sd = sd_group_1),
          rep(NA, max_n - n_group_1)
        ),
        Group2 = c(
          rnorm(n = n_group_2, mean = mean_group_1 - log2(exp(original_es)), sd = sd_group_2),
          rep(NA, max_n - n_group_2)
        )
      ) |> mutate(
        # Calculates percentages
        Group1_Perc = Group1 - mean(Group1, na.rm = T),
        Group2_Perc = Group2 - mean(Group1, na.rm = T)
      )
    } else {
      sim_data = tibble(
        LAB = LAB,
        EXP = exp_code,
        Replicate = 1:max_n,
        Group1 = c(
          rnorm(n = n_group_1, mean = mean_group_1, sd = sd_group_1),
          rep(NA, max_n - n_group_1)
        ),
        Group2 = c(
          rnorm(n = n_group_2, mean = mean_group_1 * exp(original_es), sd = sd_group_2),
          rep(NA, max_n - n_group_2)
        )
      ) |> mutate(
        # Truncates data that is negative
        Group1 = ifelse(Group1 < 0, runif (1, min=0.0001, max=0.001), Group1),
        Group2 = ifelse(Group2 < 0, runif (1, min=0.0001, max=0.001), Group2),
        # Calculates percentages
        Group1_Perc = Group1 / mean(Group1, na.rm = T),
        Group2_Perc = Group2 / mean(Group1, na.rm = T)
      )
    }
    sim_data
  })
  
  all_labs_simdata
}

# Function to transform PCR effect sizes
exp_only_PCR = function (x, b, s) {
  if (b) {
    return (log(2 ^ (s * x)))
  } else {
    return (x)
  }
}

# Function to invert CIs
switch_interval_bounds = function (df) {
  df = df |>
    mutate(
      lb0 = lb,
      ub0 = ub,
      lb = ifelse(ub0 > lb0, lb0, ub0),
      ub = ifelse(ub0 < lb0, lb0, ub0)
    ) |>
    select(-lb0, ub0)
  df
}

# Helper function to swap the values of two variables
swap = function(x,y) {
  eval(
    parse(
      text = paste(
        "swap_var = ", substitute(x), ";",
        substitute(x), " = ", substitute(y), ";",
        substitute(y), " = swap_var; rm(swap_var)")
    ),
    env=parent.frame()
  )
}
