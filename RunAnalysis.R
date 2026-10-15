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

# start log ----
release <- 'v1'
resultsFolder <- here("Results")
log_file <- paste0(here::here(resultsFolder), "/log_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), "_.txt")
logger <- create.logger()
logfile(logger) <- log_file
level(logger) <- "INFO"
source(here("Analysis", "functions.R"))

# create cdm ----
info(logger, "CREATE CDM OBJECT")
cdm <- cdmFromCon(
  con = db,
  cdmSchema = c(schema = cdmSchema),
  writeSchema = writeSchema,
  writePrefix = writePrefix,
  cdmName = dbName
)

# database specific ----
if (omopgenerics::cdmName(cdm) == "SIDIAP") {
  cdm$drug_exposure <- cdm$drug_exposure |>
    dplyr::filter(.data$drug_type_concept_id == 32839)
}

if (omopgenerics::cdmName(cdm) == "NLHR") {
  startDate <- "2018-01-01"
} else {
  startDate <- "2014-01-01"
}

# snapshot ----
info(logger, "EXTRACT SNAPSHOT OF CDM")
snapshot <- OmopSketch::summariseOmopSnapshot(cdm = cdm)

# instantiate cohorts ----
info(logger, "INSTANTIATE COHORTS")
source(here("Analysis", "instantiateCohorts.R"))

# characterisation ----
info(logger, "CHARCATERISE COHORTS OF INTEREST")
source(here("Analysis", "characterisation.R"))

# sccs ----
info(logger, "RUN SCCS")
source(here("Analysis", "SCCS.R"))

# export results ----
info(logger, "EXPORT RESULTS")
omopgenerics::exportSummarisedResult(
  snapshot,
  counts,
  attritions,
  dus1,
  dus2,
  dus3,
  indication,
  characterisationHc,
  characterisationNsaids,
  characterisationSccs,
  sccsAssumptions,
  sccsResult,
  path = resultsFolder,
  fileName = paste0("p3_c3_008_", release, "_{cdm_name}"),
  minCellCount = minCellCount
)
