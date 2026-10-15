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

# 1. Install the environment --------------------
#install.packages("renv") # if not already installed, install renv from CRAN

# this should prompt you to install the various packages required for the study
renv::activate()
renv::restore() # if asked to install or update any package, answer "y" in the console

# packages -----
library(DBI)
library(dplyr)
library(log4r)
library(CDMConnector)
library(here)
library(DrugUtilisation)
library(tidyr)
library(SqlRender)
library(CodelistGenerator)
library(PatientProfiles)
library(CohortConstructor)
library(omopgenerics)
library(readr)
library(purrr)
library(CohortCharacteristics)
library(OmopSketch)
library(clock)
library(ggplot2)
# 2. Database metadata and connection details -----
# name of the database, use acronym in capital letters (eg. "CPRD_GOLD")
dbName <- "..."

# Establishing connection
db <- dbConnect("...")

# parameters to connect to create cdm object ----
# name of the schema where cdm tables are located
cdmSchema <- "..."

# name of a schema in the database where you have writing permission
writeSchema <- "..."

# combination of at least 5 letters + _ (eg. "abcde_") that will lead any table
# written in the write schema
writePrefix <- "..."

minCellCount <- 5L

# 3. Run analysis ----
source(here("RunAnalysis.R"))
