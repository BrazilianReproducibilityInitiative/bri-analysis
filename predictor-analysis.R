### ### ### ### ### ### ### ### ### ### ### ###

# predictor-analysis.r
# This script contains the functions to run the analyses of predictors of reproducibility success.

### ### ### ### ### ### ### ### ### ### ### ###

# Function to run the experiment level predictor analysis
run_predictor_analysis_exp_level <- function(replication_datapath, make_plots) {
  # Gets replication outcomes
  REP_DATA <- read_tsv(paste0(replication_datapath, "/Replication Success by Experiment.tsv"))

  REP_DATA <- REP_DATA |>
    select(
      EXP,
      REP_Orig_in_REMA_PI,
      REP_REMA_in_Orig_CI,
      REP_FEMA_SSS_Orig,
      REP_Voting_WithTies,
      REP_Subjective_WithTies,
      SignalErrorAll,
      t_score,
      log_es_ratio,
      OriginalCV
    )

  colnames(REP_DATA) <- c(
    "EXP",
    "Original in replication PI",
    "Replication in original CI",
    "Same sense significance",
    "Voting (with ties)",
    "Subjective (with ties)",
    "Signal error",
    "t Value",
    "Log ES Ratio",
    "Original CV"
  )

  # Merge with predictor data
  dir.create(paste0(replication_datapath, "/predictors/"))
  dir.create(paste0(replication_datapath, "/predictors/by experiment"))
  FDATA <- REP_DATA |> left_join(PRED_DATA_EXP_LEVEL, by = "EXP")
  write_tsv(FDATA, paste0(replication_datapath, "/predictors/by experiment/Table for Predictor Correlations.tsv"))

  # List predictor and outcome variables
  predictors <- c(unique(colnames(PRED_DATA_EXP_LEVEL |> select(-EXP, -DOI))))
  predictors <- append(predictors, "Original CV", after = which(predictors == "p-value"))

  categorical_rep_outcome <- c(
    "Original in replication PI",
    "Replication in original CI",
    "Same sense significance",
    "Voting (with ties)",
    "Subjective (with ties)"
  )
  continuous_rep_outcome <- c(
    "t Value",
    "Log ES Ratio"
  )

  # Correlation table and plot between predictors
  print("Building correlation table, predictors x predictors ...")
  vars1 <- predictors
  vars2 <- vars1

  p_pred <- plot_cortable(FDATA, vars1, vars2, "Correlation between predictors", paste0(replication_datapath, "/predictors/by experiment/predictor x predictor"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between replication criteria
  print("Building correlation table, replication outcomes x replication outcomes ...")
  vars1 <- c(categorical_rep_outcome, continuous_rep_outcome)
  vars2 <- vars1

  p_rep <- plot_cortable(FDATA, vars1, vars2, "Correlation between replication criteria", paste0(replication_datapath, "/predictors/by experiment/rep outcome x rep outcome"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between predictors and replication outcomes
  print("Building correlation table, predictors x replication outcomes ...")
  vars1 <- predictors
  vars2 <- c(categorical_rep_outcome, continuous_rep_outcome)

  p_pred_rep <- plot_cortable_cluster(FDATA, vars1, vars2, "Correlation between predictors and replication", paste0(replication_datapath, "/predictors/by experiment/predictor x rep outcome"), show_label = T, make_individual_plots = make_plots)

  return(list(p_pred, p_rep, p_pred_rep))
}

# Function to run the replication level predictor analysis
run_predictor_analysis_rep_level <- function(replication_datapath, make_plots) {
  # Gets replication outcomes
  REP_DATA <- read_tsv(paste0(replication_datapath, "/Replication Success by Replication.tsv"))

  REP_DATA <- REP_DATA |>
    select(
      EXP,
      LAB,
      REP_Individual_Rep_in_Orig_CI,
      REP_Individual_SSS,
      REP_Individual_IndivSubjective,
      SignalErrorAll_Individual,
      t_score_individual,
      log_es_ratio_individual
    )

  colnames(REP_DATA) <- c(
    "EXP",
    "LAB",
    "Replication in original CI",
    "Same sense significance",
    "Subjective",
    "Signal error",
    "t Value",
    "Log ES Ratio"
  )

  # Merge with predictor data
  dir.create(paste0(replication_datapath, "/predictors/"))
  dir.create(paste0(replication_datapath, "/predictors/by replication/"))
  FDATA <- REP_DATA |> left_join(PRED_DATA_REP_LEVEL, by = c("EXP", "LAB"))
  write_tsv(FDATA, paste0(replication_datapath, "/predictors/by replication/Table for Predictor Correlations.tsv"))

  # List predictor and outcome variables
  predictors <- unique(colnames(PRED_DATA_REP_LEVEL |> select(-EXP, -LAB)))

  categorical_rep_outcome <- c(
    "Replication in original CI",
    "Same sense significance",
    "Subjective"
  )
  continuous_rep_outcome <- c(
    "t Value",
    "Log ES Ratio"
  )

  # Correlation table and plot between predictors
  print("Building correlation table, predictors x predictors ...")
  vars1 <- predictors
  vars2 <- vars1

  p_pred <- plot_cortable(FDATA, vars1, vars2, "Correlation between predictors", paste0(replication_datapath, "/predictors/by replication/predictor x predictor"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between replication criteria
  print("Building correlation table, replication outcomes x replication outcomes ...")
  vars1 <- c(categorical_rep_outcome, continuous_rep_outcome)
  vars2 <- vars1

  p_rep <- plot_cortable(FDATA, vars1, vars2, "Correlation between replication criteria", paste0(replication_datapath, "/predictors/by replication/rep outcome x rep outcome"), show_label = T, make_individual_plots = F)


  # Correlation table and plot between predictors and replication outcomes
  print("Building correlation table, predictors x replication outcomes ...")
  vars1 <- predictors
  vars2 <- c(categorical_rep_outcome, continuous_rep_outcome)

  p_pred_rep <- plot_cortable_cluster(FDATA, vars1, vars2, "Correlation between predictors and replication", paste0(replication_datapath, "/predictors/by replication/predictor x rep outcome"), show_label = T, make_individual_plots = make_plots)

  return(list(p_pred, p_rep, p_pred_rep))
}

# Function to run the predictor analysis in both levels (EXP and REP)
run_predictor_analysis_both_level <- function(replication_datapath, make_plots) {
  dir.create(paste0(replication_datapath, "/predictors/"))
  dir.create(paste0(replication_datapath, "/predictors/both"))

  # EXP - Gets predictor data
  FDATA_EXP <- read_tsv(paste0(replication_datapath, "/predictors/by experiment/Table for Predictor Correlations.tsv")) |>
    select(-c(2:9), -"DOI") |>
    relocate(`Original CV`, .after = `p-value`) |>
    rename(`University Ranking (Experiment)` = `Institution Ranking`)

  # REP - Gets predictor data
  FDATA_REP <- read_tsv(paste0(replication_datapath, "/predictors/by replication/Table for Predictor Correlations.tsv")) |>
    select(-c(3:8)) |>
    rename(`University Ranking (Replication)` = `Institution Ranking`)

  # Join
  FDATA <- full_join(FDATA_EXP, FDATA_REP, by = "EXP")
  write_tsv(FDATA, paste0(replication_datapath, "/predictors/both/Table for Predictor Correlations.tsv"))

  # List predictor and outcome variables
  global_vars <- setdiff(colnames(FDATA_EXP), c("EXP", "Original CV"))
  global_vars <- append(global_vars, "Original CV", after = which(global_vars == "p-value"))
  individual_vars <- setdiff(colnames(FDATA_REP), c("EXP", "LAB"))

  # Correlation table and plot between predictors
  print("Building correlation table with both levels, predictors x predictors ...")
  vars1 <- c(global_vars, individual_vars)
  vars2 <- c(global_vars, individual_vars)

  # Make the plot
  plot_cortable_alternative(FDATA, FDATA_EXP, FDATA_REP, vars1, vars2, "Correlation between predictors", paste0(replication_datapath, "/predictors/both/predictor x predictor"), show_label = T, make_individual_plots = F)
}

# Function to run a correlation between two predictors/outcomes
run_correlation = function (FDATA, c1, c2, fn, make_plot) {
  print(paste(c1, c2, sep = " X "))
  
  col1 = FDATA[[c1]]
  col2 = FDATA[[c2]]
  
  c1_cat = is.character(col1) | is.logical(col1)
  c2_cat = is.character(col2) | is.logical(col2)
  
  if (c1_cat) col1 = as.integer(as.factor(col1)) - 1
  if (c2_cat) col2 = as.integer(as.factor(col2)) - 1
  
  df = tibble(Method = FDATA$EXP |> str_extract("(MTT|EPM|PCR)"), vx = col1, vy = col2) |>
    filter(!is.na(vx), !is.na(vy))
  
  if (length(unique(df$vx)) < 2 | length(unique(df$vy)) < 2 | nrow(df) < 3) {
    return (
      tibble(
        col1 = c1,
        col2 = c2,
        rho = NA,
        p.value = NA,
        pearson_r = NA,
        pearson_p = NA,
        n = nrow(df)
      )
    )
  }
  
  corr = cor.test(col1, col2, method = "spearman", use = "pairwise.complete.obs")
  corr_pearson = cor.test(col1, col2, method = "pearson", use = "pairwise.complete.obs")
  
  if (make_plot) {
    plot_scatterplot_predictors(df, corr, c1, c2, c1_cat, c2_cat, fn)
  }
  
  tibble(
    col1 = c1,
    col2 = c2,
    rho = unname(corr$estimate),
    p.value = corr$p.value,
    pearson_r = unname(corr_pearson$estimate),
    pearson_p = corr_pearson$p.value,
    n = nrow(df)
  )
}

# Function to run a correlation between predictors and outcomes simultaneously (this function is used inside the function run_predictor_analysis_both_level)
run_correlation_alternative = function (FDATA, c1, c2, fn, make_plot, global_vars) {
  print(paste(c1, c2, sep = " X "))
  
  if(c1 %in% global_vars & c2 %in% global_vars){
    FDATA <-  FDATA |>
      group_by(EXP) |>
      summarise(across(all_of(global_vars), ~ mean(.x, na.rm = TRUE))) |>
      ungroup()
    
    col1 = FDATA[[c1]]
    col2 = FDATA[[c2]]
    
    c1_cat = is.character(col1) | is.logical(col1)
    c2_cat = is.character(col2) | is.logical(col2)
    
    if (c1_cat) col1 = as.integer(as.factor(col1)) - 1
    if (c2_cat) col2 = as.integer(as.factor(col2)) - 1
    
    df = tibble(Method = FDATA$EXP |> str_extract("(MTT|EPM|PCR)"), vx = col1, vy = col2) |>
      filter(!is.na(vx), !is.na(vy))

    if (length(unique(df$vx)) < 2 | length(unique(df$vy)) < 2 | nrow(df) < 3) {
      return (
        tibble(
          col1 = c1,
          col2 = c2,
          rho = NA,
          p.value = NA,
          pearson_r = NA,
          pearson_p = NA,
          n = nrow(df)
        )
      )
    }
    
    corr = cor.test(col1, col2, method = "spearman", use = "pairwise.complete.obs")
    corr_pearson = cor.test(col1, col2, method = "pearson", use = "pairwise.complete.obs")
    
    if (make_plot) {
      plot_scatterplot_predictors(df, corr, c1, c2, c1_cat, c2_cat, fn)
    }
    
    tibble(
      col1 = c1,
      col2 = c2,
      rho = unname(corr$estimate),
      p.value = corr$p.value,
      pearson_r = unname(corr_pearson$estimate),
      pearson_p = corr_pearson$p.value,
      n = nrow(df)
    )
  } else{
    col1 = FDATA[[c1]]
    col2 = FDATA[[c2]]
    
    c1_cat = is.character(col1) | is.logical(col1)
    c2_cat = is.character(col2) | is.logical(col2)
    
    if (c1_cat) col1 = as.integer(as.factor(col1)) - 1
    if (c2_cat) col2 = as.integer(as.factor(col2)) - 1
    
    df = tibble(Method = FDATA$EXP |> str_extract("(MTT|EPM|PCR)"), vx = col1, vy = col2) |>
      filter(!is.na(vx), !is.na(vy))
    
    if (length(unique(df$vx)) < 2 | length(unique(df$vy)) < 2 | nrow(df) < 3) {
      return (
        tibble(
          col1 = c1,
          col2 = c2,
          rho = NA,
          p.value = NA,
          pearson_r = NA,
          pearson_p = NA,
          n = nrow(df)
        )
      )
    }
    
    corr = cor.test(col1, col2, method = "spearman", use = "pairwise.complete.obs")
    corr_pearson = cor.test(col1, col2, method = "pearson", use = "pairwise.complete.obs")
    
    if (make_plot) {
      plot_scatterplot_predictors(df, corr, c1, c2, c1_cat, c2_cat, fn)
    }
    
    tibble(
      col1 = c1,
      col2 = c2,
      rho = unname(corr$estimate),
      p.value = corr$p.value,
      pearson_r = unname(corr_pearson$estimate),
      pearson_p = corr_pearson$p.value,
      n = nrow(df)
    )
  }
}
