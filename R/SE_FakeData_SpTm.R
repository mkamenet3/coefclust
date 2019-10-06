#'@title
#'South-Eastern United Status - Fake Incidence Counts across Time and Space
#'@description
#'Fake dataset of disease incidece in South-Eastern United States with time. This dataset is used for illustration of the spatio-temporal framework in the vignette.
#'
#'@docType data
#'
#'@usage data(SE_FakeData_SpTm)
#'
#'@format An object of class "data.frame". Includes the following variables:
#'\itemize{
#'    \item 1) State - Seven states of South-Eastern United States (Alabama, Florida, Georgia, Mississippi, North Carolina, South Carolina, and Tennessee); factor
#'    \item 2) FIPS - County-level FIPS codes in South-Eastern United States
#'    \item 3) long - longitude
#'    \item 4) lat - latitude
#'    \item 5) x1-x3 - covariates
#'    \item 6) y1-y3 - responses across three time periods
#'}
#'
#'@keywords datasets
#'
#'
#'@aliases SE_FakeData_SpTm
#'
#'@examples
#'data(SE_FakeData_SpTm)
#'
"SE_FakeData_SpTm"
