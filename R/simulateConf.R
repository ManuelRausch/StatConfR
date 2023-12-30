#' @importFrom plyr mdply
#' @export
simConf <- function(model = "SDT",  paramDf) {
  AllModels <-
    c('WEV', 'SDT','IG','ITGc', 'ITGcm', 'Noisy', 'PDA', 'logN')
  if (! model %in% AllModels){
    stop(paste(paste(setdiff(models, AllModels),collapse = " and "), " not implemented!"))
  }
  SimFun <- switch(model,
                   'WEV' = generateDataWEV,
                   'SDT' = generateDataSDT,
                   'IG' = generateData2Chan,
                   'ITGc' = generateDataIndTruncF,
                   'ITGcm' = generateDataIndTruncML,
                   'Noisy' = generateDataNoisy,
                   'PDA' = generateDataISDT,
                   'logN' = generateDataLognorm)
  SimData <- SimFun(paramDf)
  SimData <- cbind(paramDf, SimData)
  SimData

}



