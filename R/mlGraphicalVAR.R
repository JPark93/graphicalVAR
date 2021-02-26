# Multi-level like graphical VAR:
  
mlGraphicalVAR <- function(
  data,
  vars,
  beepvar,
  dayvar,
  idvar,
  scale = TRUE,
  centerWithin = TRUE,
  gamma = 0.5, # Gamma used in glasso in qgraph
  verbose = TRUE,
  subjectNetworks = TRUE, # Or a vector of which subjects to use!
  btwnNetworks = FALSE, # Adding an option so the between subjects networks aren't estimated for now to save time.
  redo = FALSE,
  lambda_min_kappa_fixed = 0.001,
  lambda_min_beta_fixed = 0.001,
  lambda_min_kappa = 0.05,
  lambda_min_beta = lambda_min_kappa,
  lambda_min_glasso = 0.01,
  mlgvclust = FALSE,
  ... # Args sent to graphicalVAR
){
  if (missing(idvar)) stop("'idvar' must be assigned")
  
  # Prep data:
  # if(mlgvclust == FALSE){
  dataPrepped <- tsData(data,vars=vars,beepvar=beepvar,
                        dayvar=dayvar,idvar=idvar,scale=scale,
                        centerWithin=centerWithin)
  # }else if(mlgvclust == TRUE){dataPrepped = data}
  if (verbose){
    message("Estimating fixed networks")
  }
  
  # Fixed effects:
  ResFixed <- graphicalVAR(dataPrepped, lambda_min_kappa = lambda_min_kappa_fixed, 
                           lambda_min_beta = lambda_min_beta_fixed, gamma=gamma,...)
  
  # Between-subjects:
  if(btwnNetworks == TRUE){
    if (verbose){
      message("Estimating between-subjects network")
    }
    meansData <- dataPrepped$data_means
    meansData <- meansData[,names(meansData) != idvar]
    meansData <- meansData[rowMeans(is.na(meansData))!=1,]
    ResBetween <- qgraph::EBICglasso(cov(meansData),nrow(meansData),gamma,returnAllResults = TRUE,lambda.min.ratio=lambda_min_glasso)
  }else{ResBetween = NULL
  ResBetween$optnet = NULL}
  # Computing model per person:
 
  
  IDs <- unique(dataPrepped$data[[idvar]])
  idResults <- list()
  
  # if(redo == TRUE){x = 1}else if(redo == FALSE){x = 9999}
  # if(x <= 9999){
    if (!identical(subjectNetworks,FALSE)){
      if (isTRUE(subjectNetworks)){
        subjectNetworks <- IDs
      }
      if (verbose){
        message("Estimating subject-specific networks")
        pb <- txtProgressBar(max=length(subjectNetworks),style=3)
      }
      for (i in seq_along(subjectNetworks)){
        capture.output({idResults[[i]] <- try(suppressWarnings(graphicalVAR(dataPrepped$data[dataPrepped$data[[idvar]] == subjectNetworks[i],],
                                                            vars=dataPrepped$vars,
                                                            beepvar=dataPrepped$beepvar,
                                                            dayvar=dataPrepped$dayvar,
                                                            idvar=dataPrepped$idvar,
                                                            scale = scale,
                                                            lambda_min_kappa=lambda_min_kappa,
                                                            lambda_min_beta=lambda_min_beta,
                                                            gamma=gamma,
                                                            centerWithin = centerWithin,...,verbose = FALSE)))})
        if (verbose){
          setTxtProgressBar(pb,i)
        }
        if (is(idResults[[i]], "try-error")){
          idResults[[i]] <- list()
        }
      }   
      if (verbose){
        close(pb)
      }
    } else {
      idResults <- lapply(seq_along(IDs),function(x)list())
    }
  #}else {idResults = NULL}

  # Aggregate results:
  Results <- list(fixedPCC = ResFixed$PCC, 
                  fixedPDC = ResFixed$PDC,
                  fixedbeta = ResFixed$beta,
                  fixedkappa = ResFixed$kappa,
                  fixedResults = ResFixed,
                  betweenNet = ResBetween$optnet,
                  betweenResults = ResBetween,
                  ids = IDs,
                  subjectPCC = lapply(idResults, '[[', 'PCC'),
                  subjectPDC = lapply(idResults, '[[', 'PDC'),
                  subjectbeta = lapply(idResults, '[[', 'beta'),
                  subjectkappa = lapply(idResults, '[[', 'kappa'),
                  subjecResults = idResults)
  class(Results) <- "mlGraphicalVAR"
  return(Results)
}