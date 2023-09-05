#' Data of 36 participants in a masked orientation discrimination experiment
#'
#' In each trial, participants were shown a sinusoidal grating oriented either horizontally or vertically, followed by a mask after varying stimulus-onset-asynchronies.
#' Participants were instructed to report the orientation and their degree of confidence as accurately as possible
#'
#' @docType data
#'
#' @usage data(MaskOri)
#'
#' @format A data.frame with 34,430 rows and 8 variables:
#' \describe{
#' \item{participant}{integer values as unique participant identifier}
#'   \item{stimulus}{orientation of the grating (90: vertical, 0: horizontal)}
#'   \item{correct}{0-1 column indicating whether the discrimination response was correct (1) or not (0)}
#'   \item{rating}{factor 4-point confidence scale. The four confidence categories were labelled as "not at all", "a little", "nearly sure" and "sure".}
#'   \item{diffCond}{stimulus-onset-asynchrony in ms (i.e. time between stimulus and mask onset)}
#'   \item{gender}{gender of the participant: "w" for female; "m" for male participants}
#'   \item{age}{the age of participants in years}
#'   \item{trialNo}{Enumeration of trials per participant}
#' }
#' @keywords datasets
#'
# @references TO DO
# @source \url{hier evtl OSF link  noch einfügen}

#' @examples
#' data(MaskOri)
#' summary(MaskOri)
"MaskOri"
