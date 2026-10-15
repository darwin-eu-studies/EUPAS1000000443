# Copyright 2026 European Medicines Agency
#
# This software is developed by the DARWIN EU Coordination Centre
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# This study code contains unmodified open source dependencies.

info(logger, "Prepare sccs data")
nsaidsData <- cdm$cohort_nsaids |>
  dplyr::inner_join(
    cdm$sccs_cohort |> 
      dplyr::distinct(.data$subject_id),
    by = "subject_id"
  ) |>
  PatientProfiles::addCohortName() |>
  dplyr::select(!"cohort_definition_id") |>
  dplyr::compute() 

cdm <- omopgenerics::insertTable(
  cdm = cdm, 
  name = "set",
  table = omopgenerics::settings(cdm$sccs_cohort), 
  temporary = TRUE
)

inflammatoryDisease <- c(
  "inflammatory_bowel_disease", "behcets_syndrome", "lupus_erythematosus", 
  "rheumatoid_arthritis", "cancer"
)
external = c(
  "edema_of_lower_leg", "hospitalisation", "surgery", "trauma_fracture",
  "wheelchair_reduced_mobility"
)

sccsData <- cdm$sccs_cohort |> 
  addDeathFlag(indexDate = "cohort_end_date")|>
  PatientProfiles::addDateOfBirth() |>
  dplyr::inner_join(
    cdm$set |>
      dplyr::select(
        "cohort_definition_id", "denominator_cohort" = "original_cohort",
        "exposure_cohort", "outcome_cohort"
      ),
    by = "cohort_definition_id"
  ) |>
  dplyr::inner_join(
    nsaidsData |> 
      dplyr::rename(
        exposure_cohort = "cohort_name", exposure_start = "cohort_start_date",
        exposure_end = "cohort_end_date"
      ),
    by = c("subject_id", "exposure_cohort")
  ) |>
  dplyr::filter(
    .data$exposure_start <= .data$cohort_end_date &
      .data$cohort_start_date <= .data$exposure_end 
  ) |>
  dplyr::mutate(exposure_end = dplyr::if_else(
    .data$exposure_end <= .data$cohort_end_date,
    .data$exposure_end,
    .data$cohort_end_date
  )) |>
  left_join(
    cdm$cohort_indication |>
      PatientProfiles::addCohortName() |>
      dplyr::group_by(.data$subject_id) |>
      dplyr::summarise(
        covariate_inflammatory_disease = min(.data$cohort_start_date[.data$cohort_name %in% !!inflammatoryDisease], na.rm = TRUE),
        covariate_external = min(.data$cohort_start_date[.data$cohort_name %in% !!external], na.rm = TRUE)
      ),
    by = "subject_id"
  ) |>
  dplyr::collect()

# prepare dates
sccsData <- sccsData |>
  dplyr::mutate(
    across(starts_with("covariate"), ~ as.numeric(. - date_of_birth), .names = "a{.col}"),
    aevent = as.numeric(outcome_date - date_of_birth),
    astart = as.numeric(cohort_start_date  - date_of_birth),
    aend = as.numeric(cohort_end_date  - date_of_birth),
    adrug_start = as.numeric(exposure_start  - date_of_birth),
    adrug_end = as.numeric(exposure_end  - date_of_birth),
    subject_id = as.numeric(subject_id),
    dob = format(date_of_birth, "%d%m%Y")
  ) |>
  dplyr::select(
    "aevent", "astart", "aend", "adrug_start", "adrug_end", "subject_id", "dob",
    "acovariate_inflammatory_disease", "acovariate_external", "death",
    "primary" = "denominator_cohort", "drug" = "exposure_cohort",
    "outcome" = "outcome_cohort"
  )

# generate cohort names for all nsaids and all HC comb
primary_cohorts <- c("hc_low_no_medium_no_high", "hc_medium_no_high", "hc_high")
combinations <- omopgenerics::settings(cdm$sccs_cohort) |>
  dplyr::select("primary" = "original_cohort", "drug" = "exposure_cohort", "outcome" = "outcome_cohort") |>
  dplyr::mutate(id = dplyr::row_number())

# outcome events related to exposure initiation
info(logger, "Event counts related to exposure date")
sccsAssumption1 <- sccsData |>
  dplyr::filter(.data$outcome == "vte") |>
  dplyr::mutate(days_since_nsaids = .data$aevent - .data$adrug_start) |>
  dplyr::group_by(.data$primary, .data$drug, .data$outcome, .data$days_since_nsaids) |>
  dplyr::summarise(num_times = dplyr::n(), .groups = "drop") |>
  dplyr::rename(variable_level = "days_since_nsaids") |>
  dplyr::mutate(
    cdm_name = omopgenerics::cdmName(cdm),
    variable_name = "Days between event and exposure",
    result_type = "sccs_plot_assumptions"
  ) |>
  omopgenerics::transformToSummarisedResult(
    group = c("primary", "drug", "outcome"),
    estimates = "num_times",
    settings = "result_type"
  )

# exposure events are symmetrical
info(logger, "Exposure events are symmetrical")
sccsAssumption2 <- sccsData |>
  dplyr::filter(.data$outcome == "vte") |>
  dplyr::select("primary", "drug", "outcome", "subject_id", "astart", "aend", "adrug_start") |>
  dplyr::rowwise() |>
  dplyr::mutate(nsaids_followup = list(seq(.data$astart, .data$aend, by = 1) - .data$adrug_start)) |>
  tidyr::unnest("nsaids_followup") |>
  dplyr::group_by(.data$primary, .data$drug, .data$outcome, .data$nsaids_followup) |>
  dplyr::summarise(num_times = dplyr::n(), .groups = "drop") |>
  dplyr::rename(variable_level = "nsaids_followup") |>
  dplyr::mutate(
    cdm_name = omopgenerics::cdmName(cdm),
    variable_name = "Individuals in observation respect exposure start",
    result_type = "sccs_plot_assumptions"
  ) |>
  omopgenerics::transformToSummarisedResult(
    group = c("primary", "drug", "outcome"),
    estimates = "num_times",
    settings = "result_type"
  )

# individual that die followup
info(logger, "Individuals in observation respect outcome (death = 1)")
sccsAssumption3 <- sccsData |>
  dplyr::filter(.data$outcome == "vte") |>
  dplyr::filter(.data$death == 1)
if (nrow(sccsAssumption3) > 0) {
  sccsAssumption3 <- sccsAssumption3 |>
    dplyr::select("primary", "drug", "outcome", "subject_id", "astart", "aend", "aevent") |>
    dplyr::rowwise() |>
    dplyr::mutate(outcome_followup = list(seq(.data$astart, .data$aend, by = 1) - .data$aevent)) |>
    tidyr::unnest("outcome_followup") |>
    dplyr::group_by(.data$primary, .data$drug, .data$outcome, .data$outcome_followup) |>
    dplyr::summarise(num_times = dplyr::n(), .groups = "drop") |>
    dplyr::rename(variable_level = "outcome_followup") |>
    dplyr::mutate(
      cdm_name = omopgenerics::cdmName(cdm),
      variable_name = "Individuals (that die) in observation respect outcome date",
      result_type = "sccs_plot_assumptions"
    ) |>
    omopgenerics::transformToSummarisedResult(
      group = c("primary", "drug", "outcome"),
      estimates = "num_times",
      settings = "result_type"
    )
} else {
  sccsAssumption3 <- omopgenerics::emptySummarisedResult()
}

# individual that die followup
info(logger, "Individuals in observation respect outcome (death = 0)")
sccsAssumption4 <- sccsData |>
  dplyr::filter(.data$outcome == "vte") |>
  dplyr::filter(.data$death == 0)
if (nrow(sccsAssumption4) > 0) {
  sccsAssumption4 <- sccsAssumption4 |>
    dplyr::select("primary", "drug", "outcome", "subject_id", "astart", "aend", "aevent") |>
    dplyr::rowwise() |>
    dplyr::mutate(outcome_followup = list(seq(.data$astart, .data$aend, by = 1) - .data$aevent)) |>
    tidyr::unnest("outcome_followup") |>
    dplyr::group_by(.data$primary, .data$drug, .data$outcome, .data$outcome_followup) |>
    dplyr::summarise(num_times = dplyr::n(), .groups = "drop") |>
    dplyr::rename(variable_level = "outcome_followup") |>
    dplyr::mutate(
      cdm_name = omopgenerics::cdmName(cdm),
      variable_name = "Individuals (that do not die) in observation respect outcome date",
      result_type = "sccs_plot_assumptions"
    ) |>
    omopgenerics::transformToSummarisedResult(
      group = c("primary", "drug", "outcome"),
      estimates = "num_times",
      settings = "result_type"
    )
} else {
  sccsAssumption4 <- omopgenerics::emptySummarisedResult()
}

sccsAssumptions <- omopgenerics::bind(sccsAssumption1, sccsAssumption2, sccsAssumption3, sccsAssumption4)

# sccs parameters
month <- c("0101", "0104", "0108", "0112")
covariates <- colnames(sccsData) |>
  purrr::keep(\(x) stringr::str_detect(x, "acovariate"))
age_groups_years <- c(15, 20, 25, 30, 35, 40, 45)
age_groups_days <- age_groups_years * 365.25

# 111 sccs analyses
sccsResult <- pmap(combinations, function(primary, drug, outcome, id) {
  info(logger, glue::glue("SCCS {id}/{nrow(combinations)}; denominator: {primary}, exposure: {drug}, outcome: {outcome}"))
  
  cohort_data <- sccsData |> 
    dplyr::filter(
      .data$primary == .env$primary, 
      .data$drug == .env$drug, 
      .data$outcome == .env$outcome
    )
  
  # Define the covariates and filter them based on unique values
  valid_covariates <- covariates[sapply(cohort_data[covariates], function(x) length(unique(na.omit(x))) >= 2)]
  
  cohort_data[is.na(cohort_data)] <- 1e20
  
  # Create the dynamic formula for SCCS
  dynamic_formula <- reformulate(
    c("adrug_start", "adrug_end", "astart", "age", "season", valid_covariates), 
    response = "event"
  )
  
  info(logger, "Fit full model")
  # Fit SCCS model, trycatch to avoid error
  sccsFull <- tryCatch({
    sccs_model <- SCCS::standardsccs(
      formula = dynamic_formula,
      indiv = subject_id,
      aevent = aevent,
      astart = astart,
      aend = aend,
      adrug = c(
        list(adrug_start, adrug_end + 1),
        lapply(valid_covariates, function(x) cohort_data[[x]])
      ),
      aedrug = c(
        list(adrug_end, adrug_end + 30),
        rep(list(aend), length(valid_covariates))
      ),
      expogrp = c(
        list(c(-14, 0, 1, 8)),
        list(c(0)), 
        rep(list(0), length(valid_covariates))
      ),
      data = cohort_data,
      agegrp = age_groups_days,
      dob = dob,
      seasongrp = month
    )
    
    # Format SCCS data for further processing
    sccs_data <- SCCS::formatdata(
      indiv = subject_id,
      aevent = aevent,
      astart = astart,
      aend = aend,
      adrug = list(adrug_start, adrug_end + 1),
      aedrug = list(adrug_end, adrug_end + 30),
      expogrp = list(c(-14, 0, 1, 8), c(0)),
      data = cohort_data,
      agegrp = age_groups_days
    )
    
    sccsToSummarisedResult(sccs_model, sccs_data, cdm, primary, drug, outcome, "sccs_full_model")
    
    }, error = function(e) {
      info(logger, paste("Error in full SCCS", as.character(e)))
      omopgenerics::emptySummarisedResult()
    })
  
  info(logger, "Fit age model")
  sccsAge <- tryCatch({
    sccs_model_age <- SCCS::standardsccs(
      formula = event ~ adrug_start + adrug_end + astart + age,
      indiv = subject_id,
      aevent = aevent,
      astart = astart,
      aend = aend,
      adrug = list(adrug_start, adrug_end + 1),
      aedrug = list(adrug_end, adrug_end + 30),
      expogrp = list(c(-14, 0, 1, 8), 0),
      data = cohort_data,
      agegrp = age_groups_days
    )
    
    # Format SCCS data for further processing
    sccs_data_age <- SCCS::formatdata(
      indiv = subject_id,
      aevent = aevent,
      astart = astart,
      aend = aend,
      adrug = list(adrug_start, adrug_end + 1),
      aedrug = list(adrug_end, adrug_end + 30),
      expogrp = list(c(-14, 0, 1, 8), c(0)),
      data = cohort_data,
      agegrp = age_groups_days
    )
    
    sccsToSummarisedResult(sccs_model_age, sccs_data_age, cdm, primary, drug, outcome, "sccs_age_model")

  }, error = function(e) {
    info(logger, paste("Error in age SCCS", as.character(e)))
    omopgenerics::emptySummarisedResult()
  })
  
  omopgenerics::bind(sccsFull, sccsAge)
}) |>
  omopgenerics::bind()
