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

# ------------------------ OBJ 1 ------------------------
outcomes <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "outcome")), type = "csv")
exposures <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "exposure")), type = "csv")
sensitivity <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "sensitivity")), type = "csv")
conditionToExclude <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "exclusion")), type = "csv")
covariate <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "characteristics")), type = "csv")
cormo_cs <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "comorbidities")), type = "csv")
md_cs <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "medications")), type = "csv")
preg <- importCodelist(path = normalizePath(here("..", "Codelists", "Study", "other")), type = "csv")

conceptSet <- exposures[c("any_nsaids", "ibuprofen", "diclofenac", "naproxen", "paracetamol")] |>
  omopgenerics::newCodelist()

info(logger, "Create nsaids cohort")

cdm$cohort_nsaids <- conceptCohort(cdm = cdm,
                                   conceptSet = conceptSet,
                                   name = "cohort_nsaids")

cdm$cohort_nsaids <- cdm$cohort_nsaids |> collapseCohorts(gap = 30)

info(logger, "Create hc cohort")
hc_code <- exposures[c("hc_high", "hc_medium", "hc_low")] |>
  omopgenerics::newCodelist()

cdm$hc <- conceptCohort(cdm = cdm,
                        conceptSet = hc_code,
                        name = "hc")

cdm$hc <- cdm$hc |> collapseCohorts(gap = 30)

info(logger, "Create hc low implants cohort")
nm <- omopgenerics::uniqueTableName()
cdm <- omopgenerics::insertTable(
  cdm = cdm, 
  name = nm, 
  table = dplyr::tibble(drug_concept_id = exposures$hc_low_3y, years_to_add = 3L) |>
    dplyr::union_all(dplyr::tibble(drug_concept_id = exposures$hc_low_5y, years_to_add = 5L)) |>
    dplyr::union_all(dplyr::tibble(drug_concept_id = exposures$hc_low_8y, years_to_add = 8L))
)
cdm$low_implants <- cdm$drug_exposure |>
  dplyr::inner_join(cdm[[nm]], by = "drug_concept_id") |>
  dplyr::mutate("exposure_end" = clock::add_years(.data$drug_exposure_start_date, .data$years_to_add)) |>
  dplyr::select(
    "subject_id" = "person_id", 
    "exposure_start" = "drug_exposure_start_date", 
    "exposure_end"
  ) |>
  dplyr::compute(name = "low_implants") |>
dplyr::inner_join(
    cdm$observation_period |>
      dplyr::select(
        "observation_period_id",
        "subject_id" = "person_id", 
        "obs_start" = "observation_period_start_date",
        "obs_end" = "observation_period_end_date"
      ),
    by = "subject_id"
  ) |>
  dplyr::mutate(
    cohort_start_date = dplyr::if_else(
      .data$obs_start >= .data$exposure_start, .data$obs_start, .data$exposure_start
    ),
    cohort_end_date = dplyr::if_else(
      .data$obs_end <= .data$exposure_end, .data$obs_end, .data$exposure_end
    ),
    cohort_definition_id = 1L
  ) |>
  dplyr::filter(.data$cohort_start_date <= .data$cohort_end_date) |>
  dplyr::select(dplyr::all_of(omopgenerics::cohortColumns("cohort")), "observation_period_id") |>
  dplyr::compute(name = "low_implants") |>
  CohortConstructor:::joinOverlap(name = "low_implants", gap = 30, by = c("cohort_definition_id", "subject_id", "observation_period_id")) |>
  dplyr::select(!"observation_period_id") |>
  dplyr::compute(name = "low_implants") |>
  omopgenerics::newCohortTable(cohortSetRef = dplyr::tibble(
    cohort_definition_id = 1L, cohort_name = "low_implants"
  ))

info(logger, "merge low")
cdm <- omopgenerics::bind(cdm$low_implants, cdm$hc, name = "hc")
cdm$hc2 <- cdm$hc |>
  CohortConstructor::unionCohorts(
    cohortId = c("low_implants", "hc_low"),
    gap = 30,
    cohortName = "low_merged",
    name = "hc2", 
    keepOriginalCohorts = TRUE
  )

info(logger,"remove medium+high use from low and high use from medium")
cdm$hc_high <- cdm$hc2 |>
  subsetCohorts(cohortId = "hc_high", name = "hc_high")
cdm$hc_medium <- cdm$hc2 |>
  excludeCohort(
    cohort2 = cdm$hc2, 
    id1 = "hc_medium", 
    id2 = "hc_high", 
    cohortName = "hc_medium_no_high", 
    nm = "hc_medium"
  )
cdm$hc_low <- cdm$hc2 |>
  excludeCohort(
    cohort2 = cdm$hc2, 
    id1 = "low_merged", 
    id2 = c("hc_medium", "hc_high"), 
    cohortName = "hc_low_no_medium_no_high", 
    nm = "hc_low"
  )

cdm <- bind(cdm$hc_low, cdm$hc_medium, cdm$hc_high, name = "hc")
cdm$hc <- cdm$hc |>
  omopgenerics::newCohortTable(
    cohortSetRef = omopgenerics::settings(cdm$hc) |>
      dplyr::select("cohort_definition_id", "cohort_name"), 
    .softValidation = TRUE
  )
omopgenerics::dropSourceTable(cdm = cdm, name = c(nm, "hc_low", "hc_medium", "hc_high", "hc2", "low_implants"))

info(logger, "Instantiate prior conditions to exclude")
cdm$cohort_hist_con <- conceptCohort(cdm = cdm, 
                                     conceptSet = conditionToExclude, 
                                     name = "cohort_hist_con")

info(logger, "Instantiate indication cohorts")
cdm$cohort_indication <- conceptCohort(cdm = cdm,
                                       name = "cohort_indication",
                                       conceptSet = covariate)
cancerId <- omopgenerics::getCohortId(cdm$cohort_indication, cohortName = "cancer") |>
  unname()

info(logger, "Instantiate vte cohorts")
cdm$outcome <- conceptCohort(cdm = cdm, 
                             conceptSet = outcomes, 
                             name = "outcome")
vteId <- omopgenerics::getCohortId(cdm$outcome, cohortName = "vte") |>
  unname()

info(logger, "Require no prior conditions")
nme <- omopgenerics::uniqueTableName()
subjectsToExclude <- cdm$hc |>
  dplyr::select("subject_id", "cohort_start_date") |>
  dplyr::inner_join(
    cdm$cohort_hist_con |>
      dplyr::select("subject_id", "cond_date" = "cohort_start_date") |>
      dplyr::union_all(
        cdm$cohort_indication |>
          dplyr::filter(.data$cohort_definition_id %in% .env$cancerId) |>
          dplyr::select("subject_id", "cond_date" = "cohort_start_date")
      ) |>
      dplyr::union_all(
        cdm$outcome |>
          dplyr::filter(.data$cohort_definition_id %in% .env$vteId) |>
          dplyr::select("subject_id", "cond_date" = "cohort_start_date")
      ),
    by = "subject_id"
  ) |>
  dplyr::filter(.data$cond_date <= .data$cohort_start_date) |>
  dplyr::distinct(.data$subject_id, .data$cohort_start_date) |>
  dplyr::compute(name = nme)
cdm$hc <- cdm$hc |> 
  dplyr::anti_join(subjectsToExclude, by = c("subject_id", "cohort_start_date")) |>
  dplyr::compute(name = "hc") |>
  omopgenerics::recordCohortAttrition(reason = "No prior historical condition")
omopgenerics::dropSourceTable(cdm = cdm, name = nme)

info(logger, "Create pregnancy cohorts")
cdm$pregnancy <- conceptCohort(cdm = cdm, name = "pregnancy", conceptSet = preg)

info(logger, "Restrict hc to prior obs of 365")
cdm$hc <- cdm$hc |>
  requireObservationBeforeDrug(days = 365)

info(logger, "Restrict hc to  females")
cdm$hc <- cdm$hc |>
  requireSex(sex = "Female")

info(logger, "Restrict hc aged 15-49")
cdm$hc <- cdm$hc |>
  CohortConstructor::trimDemographics(ageRange = c(15, 49))

info(logger, "Restrict hc startDate")
cdm$hc <- cdm$hc |>
  trimToDateRange(dateRange = as.Date(c(startDate, NA)))

info(logger, "Censor hc at pregnancy")
cdm$hc <- cdm$hc |>
  addCohortIntersectDate(targetCohortTable = "pregnancy", nameStyle = "pregnancy_date", name = "hc") |>
  CohortConstructor::exitAtFirstDate(c("cohort_end_date", "pregnancy_date"))

info(logger, "Instantiate comorbidities cohorts")
cdm$cohort_cormo <- conceptCohort(cdm = cdm,
                                  name = "cohort_cormo",
                                  conceptSet = cormo_cs, 
                                  subsetCohort = "hc")

info(logger, "Instantiate medications cohorts")
cdm$cohort_md <- conceptCohort(cdm = cdm,
                               name = "cohort_md",
                               conceptSet = md_cs, 
                               subsetCohort = "hc")

info(logger, "Restrict hc to  no nsaids users at index date")
cdm$hc_paracetamol <- cdm$hc |>
  CohortConstructor::copyCohorts(name = "hc_paracetamol") |>
  CohortConstructor::requireCohortIntersect(
    targetCohortTable = "cohort_nsaids",
    intersections = c(0, 0),
    targetCohortId = "paracetamol",
    window = c(-90, 0)
  )

cdm$hc <- cdm$hc |>
  CohortConstructor::requireCohortIntersect(
    targetCohortTable = "cohort_nsaids",
    intersections = c(0, 0),
    targetCohortId = "any_nsaids",
    window = c(-90, 0)
  )

info(logger, "Restrict nsaids to 365 days prior hist")
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids |>
  requireObservationBeforeDrug(days = 365, name = "cohort_nsaids_in_hc")

info(logger, "Restrict nsaids to females")
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids_in_hc |>
  requireSex(sex = "Female")

info(logger, "Restrict nsaids to age 15-49")
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids_in_hc |>
  requireAge(ageRange = list(c(15, 49)))

info(logger, "Restrict nsaids to start date")
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids_in_hc |>
  requireInDateRange(dateRange = as.Date(c(startDate, NA)))

info(logger, "Update names")
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids_in_hc |>
  CohortConstructor::copyCohorts(name = "cohort_nsaids_in_hc", n = 3)
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids_in_hc |>
  omopgenerics::newCohortTable(
    cohortSetRef = omopgenerics::settings(cdm$cohort_nsaids_in_hc) |>
      dplyr::mutate(cohort_name = dplyr::case_when(
        endsWith(.data$cohort_name, "_1") ~ gsub("_1", "_mid", .data$cohort_name),
        endsWith(.data$cohort_name, "_2") ~ gsub("_2", "_high", .data$cohort_name),
        TRUE ~ paste0(.data$cohort_name, "_low")
      )) |>
      dplyr::select("cohort_definition_id", "cohort_name")
  )

info(logger, "Intersect nsaids with hc")
cdm$cohort_nsaids_in_hc <- cdm$cohort_nsaids_in_hc |>
  CohortConstructor::requireCohortIntersect(
    cohortId = c("diclofenac_low", "ibuprofen_low", "naproxen_low"),
    targetCohortTable = "hc", 
    targetCohortId = "hc_low_no_medium_no_high",
    window = c(0, 0), 
  ) |>
  CohortConstructor::requireCohortIntersect(
    cohortId = c("diclofenac_mid", "ibuprofen_mid", "naproxen_mid"),
    targetCohortTable = "hc", 
    targetCohortId = "hc_medium_no_high",
    window = c(0, 0), 
  ) |>
  CohortConstructor::requireCohortIntersect(
    cohortId = c("diclofenac_high", "ibuprofen_high", "naproxen_high"),
    targetCohortTable = "hc", 
    targetCohortId = "hc_high",
    window = c(0, 0), 
  )

info(logger, "Sensitivity vte definition with no prior trauma/cancer/surgery/hosp")
cdm$outcome_sens1 <- cdm$outcome |> 
  requireCohortIntersect(
    targetCohortTable = "cohort_indication",
    targetCohortId = "cancer",
    intersections = c(0, 0),
    window = c(-180, 0),
    name = "outcome_sens1"
  ) |>
  requireCohortIntersect(
    targetCohortTable = "cohort_indication",
    targetCohortId = "hospitalisation",
    intersections = c(0, 0),
    window = c(-180, 0)
  ) |>
  requireCohortIntersect(
    targetCohortTable = "cohort_indication",
    targetCohortId = "surgery",
    intersections = c(0, 0),
    window = c(-180, 0)
  ) |>
  requireCohortIntersect(
    targetCohortTable = "cohort_indication",
    targetCohortId = "trauma_fracture",
    intersections = c(0, 0),
    window = c(-180, 0)
  ) |>
  omopgenerics::newCohortTable(
    cohortSetRef = omopgenerics::settings(cdm$outcome) |>
      dplyr::mutate(cohort_name = paste0(cohort_name, "_no_indication"))
  )

info(logger, "Sensitivity vte definition no prior death")
cdm$outcome_sens2 <- cdm$outcome |> 
  requireDeathFlag(window = c(0, 90), name = "outcome_sens2", negate = TRUE) |>
  omopgenerics::newCohortTable(
    cohortSetRef = omopgenerics::settings(cdm$outcome) |>
      dplyr::mutate(cohort_name = paste0(cohort_name, "_no_death"))
  )

cdm <- bind(cdm$outcome, cdm$outcome_sens1, cdm$outcome_sens2, name = "outcome")
omopgenerics::dropSourceTable(cdm = cdm, name = dplyr::starts_with("outcome_sens"))

info(logger, "Create sccs cohorts")
nsaidsNoPara <- omopgenerics::getCohortName(cdm$cohort_nsaids) |>
  purrr::keep(\(x) x != "paracetamol")
# create all combinations
combinations <- tidyr::expand_grid(
  exposure_cohort = omopgenerics::getCohortName(cdm$cohort_nsaids),
  outcome_cohort = omopgenerics::getCohortName(cdm$outcome)
) |>
  dplyr::filter(.data$exposure_cohort != "paracetamol" | .data$outcome_cohort == "vte")

# generate cohort names for all nsaids and all HC comb
pref <- omopgenerics::tmpPrefix()
nm1 <- omopgenerics::uniqueTableName(prefix = pref)
cdm[[nm1]] <- cdm$hc |>
  PatientProfiles::addCohortIntersectCount(
    targetCohortTable = "cohort_nsaids",
    targetCohortId = nsaidsNoPara,
    censorDate = "cohort_end_date",
    nameStyle = "counts_{cohort_name}",
    name = nm1
  ) |>
  PatientProfiles::addCohortIntersectDate(
    targetCohortTable = "outcome",
    censorDate = "cohort_end_date",
    nameStyle = "date_{cohort_name}",
    name = nm1
  )
nm2 <- omopgenerics::uniqueTableName(prefix = pref)
cdm[[nm2]] <- cdm$hc_paracetamol |>
  PatientProfiles::addCohortIntersectCount(
    targetCohortTable = "cohort_nsaids",
    targetCohortId = "paracetamol",
    censorDate = "cohort_end_date",
    nameStyle = "counts_{cohort_name}",
    name = nm2
  ) |>
  PatientProfiles::addCohortIntersectDate(
    targetCohortTable = "outcome",
    censorDate = "cohort_end_date",
    nameStyle = "date_{cohort_name}",
    name = nm2
  )

sccsCohorts <- purrr::pmap(combinations, \(exposure_cohort, outcome_cohort) {
  msg <- paste0("Exposure: ", exposure_cohort, "; outcome: ", outcome_cohort)
  info(logger, msg)
  cli::cli_inform(msg)
  newName <- omopgenerics::uniqueTableName(prefix = pref)
  colDate <- paste0("date_", outcome_cohort)
  colCounts <- paste0("counts_", exposure_cohort)
  nm <- ifelse(exposure_cohort == "paracetamol", nm2, nm1)
  cdm[[nm]] |>
    CohortConstructor::copyCohorts(name = newName) |>
    dplyr::filter(.data[[colCounts]] > 0) |>
    dplyr::compute(name = newName) |>
    omopgenerics::recordCohortAttrition("Exposed with: {exposure_cohort} during followup") |>
    dplyr::filter(!is.na(.data[[colDate]])) |>
    dplyr::select(
      "cohort_definition_id", "subject_id", "cohort_start_date", 
      "cohort_end_date", "number_exposure_periods" = dplyr::all_of(colCounts),
      "outcome_date" = dplyr::all_of(colDate)
    ) |>
    dplyr::compute(name = newName) |>
    omopgenerics::recordCohortAttrition("Outcome: {outcome_cohort} during followup") |>
    omopgenerics::newCohortTable(
      cohortSetRef = omopgenerics::settings(cdm$hc) |>
        dplyr::mutate(
          original_cohort = .data$cohort_name,
          cohort_name = paste0(
            .data$cohort_name, "_", .env$exposure_cohort, "_", .env$outcome_cohort
          ),
          exposure_cohort = .env$exposure_cohort,
          outcome_cohort = .env$outcome_cohort
        ) |>
        dplyr::select(
          "cohort_definition_id", "cohort_name", "original_cohort",
          "exposure_cohort", "outcome_cohort"
        )
    )
}) |>
  unname()
sccsCohorts$name <- "sccs_cohort"
cdm <- do.call(omopgenerics::bind, sccsCohorts)
omopgenerics::dropSourceTable(cdm = cdm, name = dplyr::starts_with(pref))

info(logger, "Export counts and attrition of cohorts")
attritions <- list()
counts <- list()
cohorts <- c(
  "sccs_cohort", "hc", "cohort_nsaids", "outcome", "cohort_indication", 
  "cohort_hist_con", "cohort_nsaids_in_hc"
)
for (nm in cohorts) {
  counts[[nm]] <- CohortCharacteristics::summariseCohortCount(cdm[[nm]])
  attritions[[nm]] <- CohortCharacteristics::summariseCohortAttrition(cdm[[nm]])
}

counts <- omopgenerics::bind(counts)
attritions <- omopgenerics::bind(attritions)
