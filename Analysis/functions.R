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

excludeCohort <- function(cohort1, cohort2, id1, id2, cohortName, nm) {
  id1 <- omopgenerics::validateCohortIdArgument(id1, cohort1)
  id2 <- omopgenerics::validateCohortIdArgument(id2, cohort2)
  pref <- omopgenerics::tmpPrefix()
  x <- cohort1 |>
    dplyr::filter(.data$cohort_definition_id == .env$id1) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  y <- cohort2 |>
    dplyr::filter(.data$cohort_definition_id %in% .env$id2) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  # not repeated individuals
  subjectsWithOverlap <- x |>
    dplyr::inner_join(y, by = "subject_id", suffix = c("_x", "_y")) |>
    dplyr::filter(
      .data$cohort_start_date_x <= .data$cohort_end_date_y &
        .data$cohort_start_date_y <= .data$cohort_end_date_x
    ) |>
    dplyr::distinct(.data$subject_id) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  x0 <- x |>
    dplyr::anti_join(subjectsWithOverlap, by = "subject_id") |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  xy <- x |>
    dplyr::inner_join(subjectsWithOverlap, by = "subject_id") |>
    dplyr::union_all(
      y |>
        dplyr::inner_join(subjectsWithOverlap, by = "subject_id")
    ) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  x1 <- splitOverlaps(xy, pref) |>
    dplyr::inner_join(
      xy |> 
        dplyr::rename(
          exposure_start = "cohort_start_date", exposure_end = "cohort_end_date" 
        ), 
      by = "subject_id"
    ) |>
    dplyr::filter(
      .data$exposure_start <= .data$cohort_end_date &
        .data$exposure_end >= .data$cohort_start_date
    ) |>
    dplyr::group_by(.data$subject_id, .data$id) |>
    dplyr::filter(dplyr::n() == 1L) |>
    dplyr::ungroup() |>
    dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") |>
    dplyr::filter(.data$cohort_definition_id == .env$id1) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  res <- x0 |>
    dplyr::union_all(x1) |>
    dplyr::compute(name = nm) |>
    omopgenerics::newCohortTable(
      cohortSetRef = omopgenerics::settings(cohort1) |>
        dplyr::filter(.data$cohort_definition_id == .env$id1) |>
        dplyr::mutate(cohort_name = .env$cohortName),
      cohortAttritionRef = NULL, 
      cohortCodelistRef = NULL
    )
  cdm <- omopgenerics::cdmReference(cohort1)
  omopgenerics::dropSourceTable(cdm = cdm, name = dplyr::starts_with(pref))
  return(res)
}
splitOverlaps <- function(xx, pref) {
  starts <- xx |>
    dplyr::select("subject_id", start_date = "cohort_start_date") |>
    dplyr::union_all(
      xx |>
        dplyr::select("subject_id", start_date = "cohort_end_date") |>
        dplyr::mutate(start_date = as.Date(clock::add_days(.data$start_date, 1L)))
    ) |>
    dplyr::group_by(.data$subject_id) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$start_date) |>
    dplyr::mutate(id = as.integer(dplyr::row_number())) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  ends <- xx |>
    dplyr::select( "subject_id", end_date = "cohort_end_date") |>
    dplyr::union_all(
      xx |>
        dplyr::select("subject_id", end_date = "cohort_start_date") |>
        dplyr::mutate(end_date = as.Date(clock::add_days(.data$end_date, -1L)))
    ) |>
    dplyr::group_by(.data$subject_id) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$end_date) |>
    dplyr::mutate(id = as.integer(dplyr::row_number()) - 1L) |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
  starts |>
    dplyr::inner_join(ends, by = c("subject_id", "id")) |>
    dplyr::rename(cohort_start_date = "start_date", cohort_end_date = "end_date") |>
    dplyr::compute(name = omopgenerics::uniqueTableName(pref))
}
sccsToSummarisedResult <- function(sccs_model, sccs_data, cdm, primary, drug, outcome, rt) {
  # Coefficients and IRR calculations
  coef_df <- as.data.frame(sccs_model$coefficients) |>
    tibble::rownames_to_column("Variable") |>
    dplyr::mutate(
      irr = exp(coef), 
      irr_lower = exp(coef - 1.96 * `se(coef)`), 
      irr_upper = exp(coef + 1.96 * `se(coef)`)
    ) |>
    dplyr::select("variable_level" = "Variable", "irr", "irr_lower", "irr_upper", p_value = "Pr(>|z|)") |>
    dplyr::mutate(
      primary = .env$primary, 
      drug = .env$drug, 
      outcome = .env$outcome,
      variable_name = "sccs coefficients",
      cdm_name = omopgenerics::cdmName(cdm),
      result_type = .env$rt
    ) |>
    omopgenerics::transformToSummarisedResult(
      group = c("primary", "drug", "outcome"),
      estimates = c("irr", "irr_lower", "irr_upper", "p_value"),
      settings = "result_type"
    )
  # Summarize events and person-time for each cohort
  sccs_summary_data <- sccs_data |>
    group_by(adrug_start, adrug_end) |>
    dplyr::summarise(
      event_count = as.integer(sum(event)), 
      person_time = as.integer(sum(interval)), 
      .groups = "drop"
    ) |>
    dplyr::mutate(
      variable_name = "sccs counts",
      variable_level = dplyr::case_when(
        .data$adrug_start == 1 & .data$adrug_end == 0 ~ "Pre exposure [-14, -1]",
        .data$adrug_start == 2 & .data$adrug_end == 0 ~ "Day of exposure [0, 0]",
        .data$adrug_start == 3 & .data$adrug_end == 0 ~ "Main risk period [1, 7]",
        .data$adrug_start == 4 & .data$adrug_end == 0 ~ "Extended risk [8, end exp]",
        .data$adrug_start == 0 & .data$adrug_end == 1 ~ "Post exposure [1, 30]",
        .data$adrug_start == 0 & .data$adrug_end == 0 ~ "Baseline",
        .data$adrug_start == 1 & .data$adrug_end == 1 ~ "Overlap between post and pre exposure"
      ),
      primary = .env$primary, 
      drug = .env$drug, 
      outcome = .env$outcome,
      cdm_name = omopgenerics::cdmName(cdm),
      result_type = .env$rt
    ) |>
    dplyr::select(!c("adrug_start", "adrug_end")) |>
    omopgenerics::transformToSummarisedResult(
      group = c("primary", "drug", "outcome"),
      estimates = c("event_count", "person_time"),
      settings = "result_type"
    )
  omopgenerics::bind(coef_df, sccs_summary_data)
}
