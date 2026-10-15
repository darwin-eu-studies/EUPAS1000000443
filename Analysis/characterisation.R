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

# drug utilisation
info(logger, "DrugUtilisation for ibuprofen")
dus1 <- cdm$cohort_nsaids_in_hc |>
  DrugUtilisation::summariseDrugUtilisation(
    cohortId = c("ibuprofen_low", "ibuprofen_mid", "ibuprofen_high"),
    gapEra = 30, 
    indexDate = "cohort_start_date",
    censorDate = "cohort_end_date",
    estimates = c("q25", "median", "q75", "mean", "min", "max"),
    conceptSet = conceptSet["ibuprofen"], 
    ingredientConceptId = NULL, 
    initialDailyDose = FALSE, 
    cumulativeDose = FALSE
  )

info(logger, "DrugUtilisation for diclofenac")
dus2 <- cdm$cohort_nsaids_in_hc |>
  DrugUtilisation::summariseDrugUtilisation(
    cohortId = c("diclofenac_low", "diclofenac_mid", "diclofenac_high"),
    gapEra = 30, 
    indexDate = "cohort_start_date",
    censorDate = "cohort_end_date",
    estimates = c("q25", "median", "q75", "mean", "min", "max"),
    conceptSet = conceptSet["diclofenac"], 
    ingredientConceptId = NULL, 
    initialDailyDose = FALSE, 
    cumulativeDose = FALSE
  )

info(logger, "DrugUtilisation for naproxen")
dus3 <- cdm$cohort_nsaids_in_hc |>
  DrugUtilisation::summariseDrugUtilisation(
    cohortId = c("naproxen_low", "naproxen_mid", "naproxen_high"),
    gapEra = 30, 
    indexDate = "cohort_start_date",
    censorDate = "cohort_end_date",
    estimates = c("q25", "median", "q75", "mean", "min", "max"),
    conceptSet = conceptSet["naproxen"], 
    ingredientConceptId = NULL, 
    initialDailyDose = FALSE, 
    cumulativeDose = FALSE
  )

# indication
info(logger, "Indications for nsaids")
indication <- cdm$cohort_nsaids_in_hc |>
  DrugUtilisation::summariseIndication(
    indicationCohortName = "cohort_indication",
    unknownIndicationTable = "condition_occurrence",
    indicationWindow = list(c(-7, 0), c(-30, 0)),
    mutuallyExclusive = FALSE
  )

# characterisation
characterisation <- function(cohort, cohortId = NULL, indexDate = "cohort_start_date") {
  msg <- paste0("Characterising cohort: ", omopgenerics::tableName(cohort))
  info(logger, msg)
  cohortId <- omopgenerics::validateCohortIdArgument(cohortId = cohortId, cohort = cohort)
  lsc <- cohort |>
    CohortCharacteristics::summariseLargeScaleCharacteristics(
      cohortId = cohortId,
      indexDate = indexDate,
      window = list(c(-Inf, -366), c(-365, -31), c(-30, -1)),
      eventInWindow = "condition_occurrence",
      episodeInWindow = "drug_exposure"
    )
  if (indexDate != "cohort_start_date") {
    cohort <- cohort |>
      dplyr::select("cohort_definition_id", "subject_id", cohort_start_date = dplyr::all_of(indexDate)) |>
      dplyr::mutate("cohort_end_date" = .data$cohort_start_date)
  }
  characteristics <- cohort |>
    CohortCharacteristics::summariseCharacteristics(
      cohortId = cohortId,
      counts = TRUE, 
      demographics = TRUE,
      cohortIntersectFlag = list(
        "Comorbidity any time prior" = list(
          targetCohortTable = "cohort_cormo", window = c(-Inf, 0)
        ),
        "Indication cohorts any time prior" = list(
          targetCohortTable = "cohort_indication", window = c(-Inf, 0)
        ),
        "Medications prior year" = list(
          targetCohortTable = "cohort_md", window = c(-365, 0)
        )
      )
    )
  omopgenerics::bind(lsc, characteristics)
}
characterisationHc <- characterisation(cdm$hc)
characterisationNsaids <- characterisation(cdm$cohort_nsaids_in_hc)
ids <- omopgenerics::settings(cdm$sccs_cohort) |>
  dplyr::filter(.data$outcome_cohort == "vte") |>
  dplyr::pull("cohort_definition_id")
characterisationSccs <- characterisation(cdm$sccs_cohort, ids, "outcome_date")
