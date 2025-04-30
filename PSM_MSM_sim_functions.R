fit.models1 <- function(formula = NULL, data , distr = NULL, method = "mle", bhazard = NULL, ...) {
  #call = match.call()
  exArgs <- list(...)
  exArgs$formula <- formula
  exArgs$data = data
  exArgs$bhazard =bhazard
  # exArgs$call = call
  # if (is.null(formula)) {
  #   stop("You need to specify a model 'formula', e.g. 'formula=Surv(time,event)~treat'")
  # }
  method <- tolower(method)
  if (!method %in% c("hmc", "inla", "mle")) {
    stop("Methods available for use are 'mle', 'hmc' or 'inla'")
  }
  method = check_distributions(method, distr)
  if (method == "mle") {
    res <- format_output_fit.models1(lapply(distr, function(x) runMLE1(x, 
                                                                     exArgs)), method, distr, formula, data)
  }
  return(res)
}

check_distributions <- function(method,distr) {
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  distr3 <- manipulate_distributions(distr)$distr3
  
  # If 'method' is either 'inla' or 'hmc but we're trying to run a model that is not available, then
  # falls back to 'mle'
  if(method %in% c("inla","hmc")) {
    if(!all(distr3 %in% availables[[method]])) {
      ####modelsString <- unname(labelTable[availables[[method]]])
      modelsString <- unname(manipulate_distributions(availables[[method]])$labs)
      modelsString[length(modelsString)] = paste0("or ", modelsString[length(modelsString)])
      message(paste0(
        "NB: ",toupper(method)," can only fit ",
        paste(modelsString, collapse = ", "),
        " parametric survival models. Falling back on MLE analysis")
      )
      method <- "mle"
    }
  }
  
  # 'mle' can implement all the possible models, except the PolyWeibull
  # In this case, I choose to *stop* execution, rather than falling back to 'hmc'!
  if (method == "mle") {
    if(!all(distr3 %in% availables[[method]])) {
      stop(paste0("The Poly-Weibull model is only implemented under method='hmc'.
       Please set this option in your call to 'fit.models'"))
    }
  }
  return(method)
}

load_availables <- function() {
  # INLA can only do a limited set of models (for now) so if user has selected
  # one that is not available, then falls back on MLE analysis
  availables=list(
    mle=c("genf" = "gef",
          "genf.orig" = "gof",
          "gengamma" = "gga",
          "gengamma.orig" = "ggo",
          "exp" = "exp",
          "weibull" = "wei",
          "weibullPH" = "wph",
          "lnorm" = "lno",
          "gamma" = "gam",
          "gompertz" = "gom",
          "llogis" = "llo",
          "lognormal" = "lno",
          "rps" = "rps"
    ),
    inla=c("exponential" = "exp",
           "weibull" = "wei",
           "weibullPH" = "wph",
           "lognormal" = "lno",
           "loglogistic" = "llo",
           "rps" = "rps",
           "gompertz" = "gom"      # added Mar 19, 2021
    ),
    hmc=c("Exponential" = "exp",
          "Gamma" = "gam",
          "GenF" = "gef",
          "GenGamma" = "gga",
          "Gompertz" = "gom",
          "PolyWeibull" = "pow",
          "RP" = "rps",
          "WeibullAF" = "wei",
          "WeibullPH" = "wph",
          "logLogistic" = "llo",
          "logNormal" = "lno"
    )
  )
  return(availables)
}

format_output_fit.models1 <- function(output,method,distr,formula,data) {

  # Uses the helper 'manipulated_distributions' to create the vector labs
  labs <- manipulate_distributions(distr)$labs

  # Model output
  models <- lapply(output, function(x) x$model)
  # Model fitting statistics
  model.fitting <- list(
    aic=unlist(lapply(output,function(x) x$aic)),
    bic=unlist(lapply(output,function(x) x$bic)),
    #dic=unlist(lapply(output,function(x) x$dic)),
    aicc=unlist(lapply(output,function(x) x$aicc)),
    bicc=unlist(lapply(output,function(x) x$bicc))
  )
  # Miscellanea
  misc <- list(
    time2run= unlist(lapply(output, function(x) x$time2run)),
    formula=formula,
    data=data,
    model_name=unlist(lapply(output,function(x) x$model_name))
  )
  if(any(distr=="polyweibull")) {
    misc$km=lapply(formula,function(f) make_KM(f,data))
  } else {
    misc$km=make_KM(formula,data)
  }

  # Names the elements of the list
  names(models) <- labs

  # Formats all output in a list
  res <- list(models=models,model.fitting=model.fitting,method=method,misc=misc)
  # And sets its class attribute to "survHE"
  class(res) <- "survHE"
  return(res)
}

manipulate_distributions <- function(x){
  # selected model checks -----
  matchTable = list(
    "exp" = c("exponential", "exp"),
    "wei" = c("weibull", "weibullaft", "weiaft", "waft", "weibullaf", "weiaf", "waf", "wei"),
    "wph" = c("weibullph", "weiph", "wph"),
    "gam" = c("gamma", "gam", "gma"),
    "lno" = c("lognormal", "lnormal", "lnorm", "lognorm", "lno"),
    "llo" = c("loglogistic", "loglog", "llogistic", "llogis", "llo", "llogist"),
    "gga" = c("generalisedgamma", "generalizedgamma", "ggamma", "gengamma", "gga", "ggam"),
    "ggo" = c("gengamma.orig", "ggo"),
    "gef" = c("generalisedf", "generalizedf", "genf", "gef"),
    "gof" = c("genf.orig", "gof"),
    "gom" = c("gompertz", "gpz", "gomp", "gompz", "gom"),
    "rps" = c("roystonparmar", "roystonparmarsplines", "roystonparmarspline", "spline", "splines", "rps"),
    "pow" = c("polyweibull","pow","PolyWeibull")
  )
  # Human readable label
  labelTable = c(
    "exp" = "Exponential",
    "wei" = "Weibull (AFT)",
    "wph" = "Weibull (PH)",
    "gam" = "Gamma",
    "lno" = "log-Normal", 
    "llo" = "log-Logistic",
    "gga" = "Gen. Gamma", "ggo" = "Gen. Gamma (orig parametrisation)",
    "gef" = "Gen. F", "gof" = "Gen. F (orig parametrisation)",
    "gom" = "Gompertz",
    "rps" = "Royston-Parmar",
    "pow" = "Poly-Weibull")
  # Labels used by R to define p..., r... and d... commands
  labelR = c(
    "exp" = "exp",
    "wei" = "weibull", 
    "wph" = "weibullPH", 
    "gam" = "gamma",
    "lno" = "lnorm", 
    "llo" = "llogis",
    "gga" = "gengamma", 
    "ggo" = "gengamma.orig",
    "gef" = "genf",
    "gof" = "genf.orig",
    "gom" = "gompertz",
    "rps" = "survspline",
    "pow" = "polyweibull"
  )
  
  distr = gsub("[ ]*[-]*", "", tolower(x))
  isDistrUnmatched = which(!sapply(
    1:length(distr),
    '%in%',
    unname(unlist(sapply(matchTable, match, distr)))))
  if (length(isDistrUnmatched) > 0) {
    stop(paste0("Distribution ", paste(distr[isDistrUnmatched], collapse = ", "), " could not be matched."))
  }
  
  distr3 <- numeric()
  for (i in 1:length(distr)) {
    distr3[i] <- names(which(unlist(lapply(matchTable,function(x) distr[i]%in%x))))  
  }
  labs <- unname(labelTable[distr3])
  distr <- unname(labelR[distr3])
  
  list(distr=distr,distr3=distr3,labs=labs)
}

make_KM <- function(formula,data) {
  km.formula <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  # Computes the Kaplan Meier curve using the package "rms"
  ObjSurvfit <- rms::npsurv(      # Uses the function "npsurv" from the package "rms"
    formula = km.formula,         # to fit the model specified in the "formula" object
    data = data                   # to the dataset named "data"
  )
  return(ObjSurvfit)
}

make.surv1 <- function(fit, mod=1, t=times, newdata=NULL, nsim=1,...) {
  ## Creates the survival curves for the fitted model(s)
  # fit = the result of the call to the fit.models function, containing the model fitting (and other relevant information)
  # mod = the index of the model. Default value is 1, but the user can choose which model fit to visualise, 
  #     if the call to fit.models has a vector argument for distr (so many models are fitted & stored in the same object)
  # t = the time framework to be used for the estimation of the survival curve
  # newdata = a list (of lists), specifying the values of the covariates at which the computation is performed. For example
  #           'list(list(arm=0),list(arm=1))' will create two survival curves, one obtained by setting the covariate 'arm'
  #           to the value 0 and the other by setting it to the value 1. In line with 'flexsurv' notation, the user needs
  #           to either specify the value for *all* the covariates or for none (in which case, 'newdata=NULL', which is the
  #           default). If some value is specified and at least one of the covariates is continuous, then a single survival
  #           curve will be computed in correspondence of the average values of all the covariates (including the factors, 
  #           which in this case are expanded into indicators). The order of the variables in the list *must* be the same
  #           as in the formula used for the model
  # nsim = the number of simulations from the distribution of the survival curves. Default at nsim=1, in which case
  #          uses the point estimate for the relevant distributional parameters and computes the resulting survival curve
  # ... = additional options
  
  # Defines list with optional parameters
  exArgs <- list(...)
  
  # Extracts the model object and the data from the survHE output
  m <- fit$models[[mod]]
  data <- fit$misc$data
  
  # By default uses the mean to compute summary statistics (eg for the case when nsim=1)
  # but the user can specify the median, which works better for very skewed parameters
  # which happens for example when trying to fit a model with too many parameters, which
  # results in a huge uncertainty, blowing up the estimates and making even the mean survival
  # curve impossible to compute
  if (exists("summary_stat",exArgs)){
    summary_stat <- exArgs$summary_stat
  } else {summary_stat <- "mean"}
  
  # Makes sure the distribution name(s) vector is in a usable format
  dist <- fit$misc$model_name[mod]
  
  # Now creates the profile of covariates for which to compute the survival curves
  X <- make_profile_surv(fit$misc$formula, data, newdata)
  
  # This is needed to rescale correctly the INLA models (which are fitted
  # on a range [0-1] for numerical stability)
  time_max <- max(times)
  
  # Draws a sample of nsim simulations from the distribution of the model parameters
  sim <- do.call(paste0("make_sim_",fit$method),
                 args=list(m=m,t=t,X=X,nsim=nsim,newdata=newdata,dist=dist,data=data,
                           formula=fit$misc$formula,summary_stat=summary_stat,time_max=time_max)
  )
  # Computes the survival curves - first in matrix form with all the simulations
  # Needs to add more inputs for the case of hmc/rps
  if (fit$method=="hmc" && dist=="rps") {
    exArgs$data.stan <- fit$misc$data.stan[[mod]]
    t[t==0] <- min(0.00001,min(t[t>0]))
  }
  if (fit$method=="mle" && dist=="rps") {
    exArgs$knots <- fit$models[[mod]]$knots
  }
  
  mat <- do.call(compute_surv_curve,
                 args=list(sim=sim,exArgs=exArgs,nsim=nsim,
                           dist=dist,t=t,method=fit$method,X=X) 
  )
  
  
  # Finally computes the actual survival curves, in summary forms
  if (nsim == 1) {
    # If nsim=1 then only save the point estimates of the survival curve
    S <- lapply(mat, function(x) {
      rowwise(x, time) |>
        summarise(S = mean(c_across(contains("S")))) |> ungroup()
    })
  } else {
    # If nsim>1 then also give the lower and upper quartile of the underlying distribution
    # avoids no visibile binding for global variable
    
    S <- lapply(mat, function(x) {
      rowwise(x, time) |>
        summarise(S = mean(c_across(contains("S"))),
                  low = quantile(c_across(contains("S")), 0.025, na.rm = T),
                  upp = quantile(c_across(contains("S")), 0.975, na.rm = T)) |> ungroup()
    })
  }
  
  list(
    S = S,
    sim = sim,
    nsim = nsim,
    mat = mat,
    des.mat = X,
    times = t)
}

compute_surv_curve <- function(sim,exArgs,nsim,dist,t,method,X) {  
  # Computes the survival curves
  args <- args_surv()
  distr <- manipulate_distributions(dist)$distr 
  
  if (dist=="rps") {
    # RPS-related options
    if (exists("scale",where=exArgs)) {scale <- exArgs$scale} else {scale <- "hazard"}
    if (exists("timescale",where=exArgs)) {timescale <- exArgs$timescale} else {timescale <- "log"}
    if (exists("log",where=exArgs)) {log <- exArgs$log} else {log <- FALSE}
    
    if (method=="hmc") {
      knots <- exArgs$data.stan$knots
      mat <-
        lapply(sim, function(x) {
          gamma=as_tibble(x) %>% select(contains("gamma"))
          offset=as_tibble(x) %>% select(offset)
          matrix(
            unlist(
              lapply(1:nsim,function(i){
                1-do.call(psurvspline,args=list(
                  q=t,
                  gamma=as.numeric(gamma %>% slice(i)),
                  beta=0,
                  X=0,
                  knots=knots,
                  scale=scale,
                  timescale=timescale,
                  offset=as.numeric(offset%>% slice(i)),
                  log=log
                ))
              })
            ),nrow=length(t),ncol=nsim,byrow=FALSE
          )
        })
    } 
    if (method=="mle") {
      # First needs to fiddle with the matrix of covariates profile
      if ("(Intercept)" %in% colnames(X)){
        X <- X %>% as_tibble() %>% select(-"(Intercept)")
      } else {
        X <- X %>% as_tibble()
      }
      ###      if(exists("offset",where=exArgs)) {offset=exArgs$offset} else {offset=0}
      mat <- 
        lapply(sim,function(x) {
          gamma=as_tibble(x) |> select(contains("gamma"))
          matrix(unlist(
            lapply(1:nsim,function(i) {
              1-do.call(psurvspline,args=list(
                q=t,
                gamma=as.numeric(gamma |> slice(i)),
                beta=0,
                X=0,
                knots=exArgs$knots,
                scale=scale,
                timescale=timescale,
                offset=0,
                log=log
              ))
            })
          ),nrow=length(t),ncol=nsim,byrow=FALSE
          )
        })
    }
  } else {
    mat <- lapply(sim, function(x) {
      matrix(
        unlist(
          lapply(1:nsim, function(i) {
            1 - do.call(paste0("p",distr),args=eval(parse(text=args[[distr]])))
          })
        ),nrow=length(t),ncol=nsim,byrow=FALSE
      )
    })
  }
  for (i in seq_along(mat)) {
    colnames(mat[[i]]) <- paste0("S_", 1:nsim)
  }
  mat <- lapply(mat, function(x) bind_cols(tibble(time = t),
                                           as_tibble(x)))
  
  return(mat)
}

args_surv <- function() {
  list(
    exp='list(t,rate=x[,"rate"][i])',
    weibull='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
    weibullPH='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
    gamma='list(t,shape=x[,"shape"][i],rate=x[,"rate"][i])',
    gengamma='list(t,mu=x[,"mu"][i],sigma=x[,"sigma"][i],Q=x[,"Q"][i])',
    gompertz='list(t,shape=x[,"shape"][i],rate=x[,"rate"][i])',
    lnorm='list(t,meanlog=x[,"meanlog"][i],sdlog=x[,"sdlog"][i])',
    llogis='list(t,shape=x[,"shape"][i],scale=x[,"scale"][i])',
    genf='list(t,mu=x[,"mu"][i],sigma=x[,"sigma"][i],Q=x[,"Q"][i],P=x[,"P"][i])',
    rps='list(t,gamma=x[,"gamma"][i],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)',
    survspline='list(t,gamma=x[,"gamma"][i],knots=m$knots,scale=scale,timescale=timescale,offset=offset,log=log)'
  )
}

make_profile_surv <- function(formula,data,newdata) {
  # Checks how many elements are given in 'newdata'
  n.elements <- ifelse(is.null(newdata),0,length(newdata))
  n.provided <- unlist(lapply(newdata,function(x) length(x)))
  
  # Temporarily re-writes the model formula to avoid issues with naming
  formula_temp <- update(formula,paste(all.vars(formula,data)[1],"~",all.vars(formula,data)[2],"+."))
  # Creates a tibble with all the covariates in their original format
  covs <- data %>% model.frame(formula_temp,.) %>% as_tibble(.) %>%  select(-c(1:2)) %>% 
    rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
    rename_if(is.factor,.funs=~gsub("[( )]","",.x)) %>% 
    bind_cols(as_tibble(model.matrix(formula_temp,data)) %>% select(contains("Intercept"))) 
  if("(Intercept)"%in% names(covs)) {
    covs=covs %>% select(`(Intercept)`,everything())
  }
  
  ncovs <- covs %>% select(-contains("Intercept")) %>% with(ncol(.))
  # Selects the subset of categorical covariates
  is.fac <- covs %>% select(where(is.factor))
  fac.levels=lapply(is.fac,levels)
  nfacts <- covs %>% select(where(is.factor)) %>% with(ncol(.))
  
  # If formula is in 'inla' terms now change it back to 'flexsurv' terms
  formula_temp <- as.formula(gsub("inla.surv","Surv",deparse(formula)))
  # Computes the "average" profile of the covariates
  X <- data %>% model.matrix(formula_temp,.) %>% as_tibble(.) %>% summarise_all(mean) 
  # If there's at least one factor with more than 2 levels, then do *not* rename the columns to the simpler version
  # which only has the name of the variable (rather than the combination, eg 'groupMedium' that R produces)
  if(all(unlist(lapply(fac.levels,length))<=2)) {colnames(X)=colnames(covs)}
  
  # The way the object X *must* be formatted depends on which way it's been generated.
  # The point is that it *always* has to be a matrix for other functions to process
  if(n.elements==0){
    # If all the covariates are factors, then get survival curves for *all* the combinations
    if(nfacts==ncovs & nfacts>0) {
      X=unique(model.matrix(formula,data))
      # If there's at least one factor with more than 2 levels, then do *not* rename the columns to the simpler version
      # which only has the name of the variable (rather than the combination, eg 'groupMedium' that R produces)
      if(all(unlist(lapply(fac.levels,length))<2)) {colnames(X)=colnames(covs)}
    } else {
      X <- as.matrix(X,nrow=nrow(X),ncol=ncol(X))
    }
  }
  
  # If 'newdata' provides a specific (set of) profile(s), then use that
  if (n.elements>=1) {
    # This means that n.provided will also be > 0 but needs to check it's the right number and if not, stop
    if (!all(n.provided==ncovs)) {
      stop("You need to provide data for *all* the covariates specified in the model, in the list 'newdata'")
    } else {
      # Turns the list 'newdata' into a data.frame & ensures the factors are indeed factors
      nd=do.call(rbind.data.frame,newdata) %>% as_tibble() %>% 
        mutate(across(names(is.fac),factor))
      # Augments the original data with the values supplied in 'newdata'. To make things work, 
      # needs to change the type of variables that are 'factors' in the analysis as well as remove
      # the NAs and replace with 0s
      aug_data=data %>% mutate(across(names(is.fac),factor)) %>% add_row(nd) %>% replace(is.na(.),0)
      # Creates a model frame with IDs
      mf=model.frame(formula,aug_data) %>% as_tibble() %>% select(-1) %>% 
        mutate(id=row_number()) %>% rename_if(is.factor,.funs=~gsub("as.factor[( )]","",.x)) %>% 
        rename_if(is.factor,.funs=~gsub("[( )]","",.x))
      # Now creates the 'model matrix' with the combination of all the factors
      mm=model.matrix(formula,aug_data) %>% as_tibble() %>% mutate(id=row_number())
      mf=suppressMessages(mf %>% right_join(nd))
      # And selects only the rows that match with the profile selected in 'newdata'
      X=as.matrix(mm %>% filter(id %in% mf$id) %>% select(-id) %>% unique,drop=FALSE)
    }
  }
  
  # Returns the output (the matrix with the covariates profile)
  return(X)
}

make_sim_mle <- function(m,t,X,nsim,newdata,dist,summary_stat,...) {
  # Simulates from the distribution of the model parameters - takes 100000 bootstrap samples
  nboot=100000
  B=ifelse(nsim<nboot,nboot,nsim)
  if(is.null(newdata)) {
    #####X=X %>% as_tibble()
    # NB: 'flexsurv' needs to exclude the intercept
    if(grep("Intercept",colnames(X))>0) {
      # If the intercept is part of the design matrix X then remove it (to make normboot work!)
      X=matrix(X[,-grep("Intercept",colnames(X))],nrow=nrow(X))
    } 
    # If X has only one row, needs to create a list, with length equal to the number of profiles (=nrow(X))
    if(nrow(X)==1) {
      sim=list(flexsurv::normboot.flexsurvreg(m,B=B,X=as.matrix(X)))
    } else {
      # Otherwise normboot will take care of it with the proper length for the automatically created list
      sim=flexsurv::normboot.flexsurvreg(m,B=B,X=as.matrix(X))
    }
  } else {
    # If there are newdata, then create the list of sims using it
    sim <- lapply(1:nrow(X),function(i) flexsurv::normboot.flexsurvreg(m,B=B,newdata=newdata[[i]]))
  }
  # Then if 'nsim'=1, then take the average over the bootstrap samples. 
  if(nsim==1) {
    sim=lapply(sim,function(x) x %>% as_tibble() %>% summarise_all(summary_stat) %>% as.matrix(.,ncol=ncol(X)))
  } 
  # If nsim<=5000 (number of bootstrap samples), then samples only 'nsim' of them
  if(nsim>1 & nsim<nboot) {
    sim=lapply(sim,function(x) x %>% as_tibble %>% sample_n(ifelse(nsim<nboot,nsim,B),replace=FALSE) %>% 
                 as.matrix(.,nrow=nsim,ncol=ncol(X)))
  }
  return(sim)
}

runMLE1 <- function(x,exArgs) {
  
  # Loads the model formula & data
  formula <- exArgs$formula
  data=exArgs$data
  bhazard = exArgs$bhazard
  
  # Loads in the available models in each method
  availables <- load_availables()
  # Uses the helper 'manipulated_distributions' to create the vectors distr, distr3 and labs
  d3 <- manipulate_distributions(x)$distr3
  x <- manipulate_distributions(x)$distr
  
  tic <- proc.time()
  
  # If it's one of the other available models under MLE, then simply runs flexsurv::flexsurvreg
  # But allows to use weight and subset
  if(exists("weights",where=exArgs)) {weights <- exArgs$weights} else {weights=NULL}
  if(exists("subset",where=exArgs)) {subset <- exArgs$subset} else {subset <- NULL}
  # if(exists("bhazard",where=exArgs)) {bhazard <- exArgs$bhazard} else {bhazard <-NULL}    
  ###model <- flexsurv::flexsurvreg(formula=formula,data=data,dist=x,weights=weights)
  model <- tryCatch(do.call("flexsurvreg",args=list(
    formula=formula,
    data=quote(data),
    dist=x,
    weights=weights,
    subset=subset,
    bhazard=bhazard
  )), error = function(e) NA, warning = function(w) list(NA))
  # Fix the 'call' by using parts of the original call to fit.models
  # if(!is.null(weights)) {model$call$weights=exArgs$call$weights} else {model$call$weights=NULL}
  # if(!is.null(subset)) {model$call$subset=exArgs$call$subset} else {model$call$subset=NULL}
  # if(!is.null(bhazard)) {model$call$bhazard=exArgs$call$bhazard} else {model$call$bhazard=NULL}
  
  toc <- proc.time()-tic
  
  # Replaces a field used in 'make.surv' to indicate the model used (standardised across models)
  model_name <- d3
  
  # Finally returns the output
  
  if (is_na(model)) {
    res<- list(
      model=NA,
      aic=NA,
      bic=NA,
      aicc=NA,
      bicc=NA,
      time2run=toc[3],
      model_name=model_name
    )
  }  else{
    res <- list(
      model=model,
      aic=model$AIC,
      bic=-2*model$loglik+model$npars*log(model$N),
      aicc = ifelse(anyNA(model), NA, - 2 * model$loglik + model$npars * 2 * model$N / (model$N - model$npars - 1)),
      bicc = ifelse(anyNA(model), NA, - 2 * model$loglik + model$npars * log(model$N) * model$N / (model$N - model$npars - 1)),
      time2run=toc[3],
      model_name=model_name
    )
  }
  return(res)
}

Probs <- function(M_t, v_Ts, t, p_HS, p_HD, p_SD) { 
  # Arguments:
  # M_t: health state occupied by at cycle t (character variable)
  # v_Ts: vector with the duration of being Sick
  # t:     current cycle 
  # Returns: 
  # transition probabilities for that cycle
  
  # create matrix of state transition probabilities
  m_p_t           <- matrix(0, nrow = n_states, ncol = n_i) 
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states                              
  
  # update m_p_t with the appropriate probabilities   
  # transition probabilities when Healthy 
  m_p_t[, M_t == "pfs"]     <- rbind( (1 - p_HD[t]) * (1 - p_HS[t]),
                                      (1 - p_HD[t]) *      p_HS[t],
                                      p_HD[t])     
  # transition probabilities when Sick 
  m_p_t[, M_t == "pd"]      <- rbind(0,
                                     1 - p_SD[v_Ts],
                                     p_SD[v_Ts])  
  # transition probabilities when Dead     
  m_p_t[, M_t == "dead"]    <- rbind(0, 0, 1)                            
  return(t(m_p_t))
}       

trans_prob <- function(surv){
  # t.p <- 1- surv[-1]/(surv[-length(surv)])
  d_surv <- surv[-1]/(surv[-length(surv)])
  t.p <- 1 - d_surv
  # if (sum(t.p < 0) > 0) {
  #   message("Negative transition probabilities were set to 0.")
  # }
  d_surv[d_surv > 1] <- 1
  t.p <- 1 - d_surv
  return(t.p = t.p)
}

partsurv_bhaz <- function(pfs_survHE = NULL, os_survHE = NULL, l_d.data = NULL, l_vc.data = NULL, par = FALSE, chol = FALSE,
                          choose_PFS = NULL, choose_OS = NULL, time = times, v_names_states, PA = FALSE, n_sim = 1000){
    
  deter <- ifelse(PA == 1, 0, 1) # determine if analysis is deterministic or probabilistic
  dist_PFS <- choose_PFS
  dist_OS  <- choose_OS
  chosen_models <- paste0("PFS: ", dist_PFS, ", ", "OS: ", dist_OS) # chosen model names
  
  # Calculate survival probabilities
  if (deter == 0) { # probabilistic
    # use survival models
    # Model-setup
    # model objects
    pfs_survHE <- pfs_survHE$model.objects
    os_survHE <-  os_survHE$model.objects
    # model names
    mod.pfs <- names(pfs_survHE$models)
    mod.os <- names(os_survHE$models)
    # chosen model index based on name
    mod.pfs.chosen <- which(mod.pfs == dist_PFS)
    mod.os.chosen  <- which(mod.os == dist_OS)
    fit_PFS <- make.surv1(pfs_survHE,
                         mod = mod.pfs.chosen,
                         nsim = n_sim,
                         t = times)
    fit_OS  <- make.surv1(os_survHE,
                         mod = mod.os.chosen,
                         nsim = n_sim,
                         t = times)
    pfs.surv <- surv_prob(fit_PFS, PA = TRUE)
    os.surv.rel  <- surv_prob( fit_OS, PA = TRUE)
  }
  else { # deterministic
    # use survival models
    # Model-setup
    # model objects
    pfs_survHE <- pfs_survHE$model.objects
    os_survHE <-  os_survHE$model.objects
    # model names
    mod.pfs <- names(pfs_survHE$models)
    mod.os <- names(os_survHE$models)
    pfs.surv <- surv_prob(pfs_survHE$models[[which(mod.pfs == dist_PFS)]], time = times)
    os.surv.rel  <- surv_prob(os_survHE$models[[which(mod.os  ==  dist_OS)]], time = times)
  }
  
  # if PFS > OS, make PFS = OS
 # check_PFS_OS(os.surv.rel * os.surv.genpop - pfs.surv)          # print warning message if PFS > OS
  if (deter == 0) { # probabilistic
    pfs.surv <- as.matrix(pfs.surv)
    os.surv.rel  <- as.matrix(os.surv.rel)
    os.surv <- os.surv.rel
    # 
    for (i in 1:ncol(pfs.surv)) {
      os.surv[,i] <- os.surv.rel[,i]*os.surv.genpop
      pfs.surv[,i][pfs.surv[,i] > os.surv[,i]] <- os.surv[,i][pfs.surv[,i] > os.surv[,i]]
    }
  } else { # deterministic
    os.surv <- os.surv.rel * os.surv.genpop
    pfs.surv[pfs.surv > os.surv] <- os.surv[pfs.surv > os.surv]
  }
  
  # Calculate state occupation proportions
  Sick                 <- os.surv - pfs.surv    # estimate the probability of remaining in the progressed state
  Healthy              <- pfs.surv              # probability of remaining stable
  Dead                 <- 1 - os.surv           # probability of being Dead
  trace <- abind(Healthy,
                 Sick,
                 Dead, rev.along = 0)
  # Calculate area under survival curve (expected survival)
  if (deter == 0) { # probabilistic
    pfs.expected.surv <- mean(apply(pfs.surv, 2, function(x){expected_surv(time = times, surv = x)}))
    os.expected.surv  <- mean(apply(os.surv,  2, function(x){expected_surv(time = times, surv = x)}))
  } else { # deterministic
    pfs.expected.surv <- expected_surv(time = times, surv = pfs.surv)
    os.expected.surv  <- expected_surv(time = times, surv = os.surv)
  }
  
  # Calculate transition probabilities
  if (deter == 0) { # probabilistic
    pfs.trans.prob <- apply(pfs.surv, 2, trans_prob)
    os.trans.prob  <- apply( os.surv, 2, trans_prob)
  } else { # deterministic
    pfs.trans.prob <- trans_prob(pfs.surv)
    os.trans.prob  <- trans_prob(os.surv)
  }
  
  # Calculate Markov trace
  if (deter == 0){ # probabilistic
    trace      <- aperm(trace, perm = c(1,3,2)) # Markov trace of all simulations
    mean.trace <- apply(trace, 1:2, mean, na.rm = TRUE)
    CI <- apply(trace, 1:2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    CI <- aperm(CI, perm = c(2, 3, 1))
    dimnames(mean.trace)[[2]] <- v_names_states
    dimnames(CI)[[3]] <- c("low", "high")
    dimnames(CI)[[2]] <- v_names_states
    dimnames(trace)[[2]] <- v_names_states
    
    # List of items to return
    res <- list(trace = trace, CI = CI, Mean = mean.trace, # Markov trace 
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  } else { # deterministic
    dimnames(trace)[[2]] <- v_names_states
    
    # List of items to return
    res <- list(trace = trace, # Markov trace
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  }
  
  return(res)
  
}

partsurv_nobhaz <- function(pfs_survHE = NULL, os_survHE = NULL, l_d.data = NULL, l_vc.data = NULL, par = FALSE, chol = FALSE,
                          choose_PFS = NULL, choose_OS = NULL, time = times, v_names_states, PA = FALSE, n_sim = 1000){
  
  deter <- ifelse(PA == 1, 0, 1) # determine if analysis is deterministic or probabilistic
  dist_PFS <- choose_PFS
  dist_OS  <- choose_OS
  chosen_models <- paste0("PFS: ", dist_PFS, ", ", "OS: ", dist_OS) # chosen model names
  
  # Calculate survival probabilities
  if (deter == 0) { # probabilistic
    # use survival models
    # Model-setup
    # model objects
    pfs_survHE <- pfs_survHE$model.objects
    os_survHE <-  os_survHE$model.objects
    # model names
    mod.pfs <- names(pfs_survHE$models)
    mod.os <- names(os_survHE$models)
    # chosen model index based on name
    mod.pfs.chosen <- which(mod.pfs == dist_PFS)
    mod.os.chosen  <- which(mod.os == dist_OS)
    fit_PFS <- make.surv1(pfs_survHE,
                         mod = mod.pfs.chosen,
                         nsim = n_sim,
                         t = times)
    fit_OS  <- make.surv1(os_survHE,
                         mod = mod.os.chosen,
                         nsim = n_sim,
                         t = times)
    pfs.surv <- surv_prob(fit_PFS, PA = TRUE)
    os.surv  <- surv_prob( fit_OS, PA = TRUE)
  }
  else { # deterministic
    # use survival models
    # Model-setup
    # model objects
    pfs_survHE <- pfs_survHE$model.objects
    os_survHE <-  os_survHE$model.objects
    # model names
    mod.pfs <- names(pfs_survHE$models)
    mod.os <- names(os_survHE$models)
    pfs.surv <- surv_prob(pfs_survHE$models[[which(mod.pfs == dist_PFS)]], time = times)
    os.surv  <- surv_prob(os_survHE$models[[which(mod.os  ==  dist_OS)]], time = times)
  }
  
  # if PFS > OS, make PFS = OS
  # check_PFS_OS(os.surv - pfs.surv)          # print warning message if PFS > OS
  if (deter == 0) { # probabilistic
    pfs.surv <- as.matrix(pfs.surv)
    os.surv  <- as.matrix(os.surv)
    # 
    for (i in 1:ncol(pfs.surv)) {
      pfs.surv[,i][pfs.surv[,i] > os.surv[,i]] <- os.surv[,i][pfs.surv[,i] > os.surv[,i]]
    }
  } else { # deterministic
    pfs.surv[pfs.surv > os.surv] <- os.surv[pfs.surv > os.surv]
  }
  
  # Calculate state occupation proportions
  Sick                 <- os.surv - pfs.surv    # estimate the probability of remaining in the progressed state
  Healthy              <- pfs.surv              # probability of remaining stable
  Dead                 <- 1 - os.surv           # probability of being Dead
  trace <- abind(Healthy,
                 Sick,
                 Dead, rev.along = 0)
  # Calculate area under survival curve (expected survival)
  if (deter == 0) { # probabilistic
    pfs.expected.surv <- mean(apply(pfs.surv, 2, function(x){expected_surv(time = times, surv = x)}))
    os.expected.surv  <- mean(apply(os.surv,  2, function(x){expected_surv(time = times, surv = x)}))
  } else { # deterministic
    pfs.expected.surv <- expected_surv(time = times, surv = pfs.surv)
    os.expected.surv  <- expected_surv(time = times, surv = os.surv)
  }
  
  # Calculate transition probabilities
  if (deter == 0) { # probabilistic
    pfs.trans.prob <- apply(pfs.surv, 2, trans_prob)
    os.trans.prob  <- apply( os.surv, 2, trans_prob)
  } else { # deterministic
    pfs.trans.prob <- trans_prob(pfs.surv)
    os.trans.prob  <- trans_prob(os.surv)
  }
  
  # Calculate Markov trace
  if (deter == 0){ # probabilistic
    trace      <- aperm(trace, perm = c(1,3,2)) # Markov trace of all simulations
    mean.trace <- apply(trace, 1:2, mean)
    CI <- apply(trace, 1:2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    CI <- aperm(CI, perm = c(2, 3, 1))
    dimnames(mean.trace)[[2]] <- v_names_states
    dimnames(CI)[[3]] <- c("low", "high")
    dimnames(CI)[[2]] <- v_names_states
    dimnames(trace)[[2]] <- v_names_states
    
    # List of items to return
    res <- list(trace = trace, CI = CI, Mean = mean.trace, # Markov trace 
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  } else { # deterministic
    dimnames(trace)[[2]] <- v_names_states
    
    # List of items to return
    res <- list(trace = trace, # Markov trace
                pfs.expected.surv = pfs.expected.surv, os.expected.surv = os.expected.surv, # expected survival
                pfs.surv = pfs.surv, os.surv = os.surv, # survival probabilities
                pfs.trans.prob = pfs.trans.prob, os.trans.prob = os.trans.prob, # transition probabilities
                chosen_models = chosen_models # chosen model names
    )
  }
  
  return(res)
  
}


fit_function_bhaz <- function(data) {
  
  fit_pfs  <- list(model.objects = fit.models1(
    formula = Surv(time_pfs2, status_pfs2) ~ 1,
    data = data,
    distr = distribs,
    bhazard = NULL
    ))
  
  bhaz <- tryCatch(flexsurv::hgompertz(data$time_os2, shape = mort_reg$coefficients[2], rate = mort_reg$coefficients[1]),
                   error = function(e) list(NA), warning = function(w) list(NA))
  
  fit_os_rel  <- list(model.objects = fit.models1(
    formula = Surv(time_os2, status_os2) ~ 1,
    data = data,    
    distr = distribs,
    bhazard = bhaz))
  
  pfs.events <- sum(data$status_pfs2, na.rm = TRUE)
  os.events <- sum(data$status_os2, na.rm = TRUE)
  f.up = max(data$time_pfs2, na.rm = TRUE)
  
  
  list(setNames(sapply(seq_along(c("aicc", "bicc")), function(i){
    
    fit_pfs$bestfit <- tryCatch(
      names(fit_pfs$model.objects$models[which.min(fit_pfs[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    fit_os_rel$bestfit <- tryCatch(
      names(fit_os_rel$model.objects$models[which.min(fit_os_rel[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
if (is_na(fit_os_rel$bestfit)) {
  res1 <- list(pfs = NA, 
               pfs_ci = NA,
               pd = NA,
               pd_ci = NA,
               os = NA, 
               os_ci = NA,
               chosen_dist = NA,
               pfs_dist = fit_pfs$bestfit,
               os_dist = fit_os_rel$bestfit,
               pfs_events = pfs.events,
               os_events = os.events,
               f_up = f.up,
               pfs_obs = NA,
               pd_obs = NA,
               os_obs = NA)
  
} else {
    
    parts1 <- tryCatch(partsurv_bhaz(
      pfs_survHE = fit_pfs,
      os_survHE = fit_os_rel,
      choose_PFS = fit_pfs$bestfit,
      choose_OS = fit_os_rel$bestfit,
      time = times,
      PA = FALSE,
      #seed = seed,
      v_names_states = v_names_states),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    
      partsa <- tryCatch(partsurv_bhaz(
      pfs_survHE = fit_pfs,
      os_survHE = fit_os_rel,
      choose_PFS = fit_pfs$bestfit,
      choose_OS = fit_os_rel$bestfit,
      time = times,
      PA = TRUE,
      n_sim = 1000,
      #seed = seed,
      v_names_states = v_names_states),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    if (is_na(partsa)){
      
      res1 <- list(pfs = NA, 
                   pfs_ci = NA,
                   pd = NA,
                   pd_ci = NA,
                   os = NA, 
                   os_ci = NA,
                   chosen_dist = NA,
                   pfs_dist = fit_pfs$bestfit,
                   os_dist = fit_os_rel$bestfit,
                   pfs_events = pfs.events,
                   os_events = os.events,
                   f_up = f.up,
                   pfs_obs = NA,
                   pd_obs = NA,
                   os_obs = NA,
                   pfsdet = parts1$pfs.expected.surv, 
                   osdet = parts1$os.expected.surv)
      
    } else {
    ## Calculate mean time in each state and related CI
    pd.surv <- apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times, surv = x)}) - apply(partsa$pfs.surv, 2, function(x){
        expected_surv(time = times, surv = x)})
    
    pd.expected.surv <- mean(pd.surv)
    
    pfs.expected.surv.ci <- quantile(apply(partsa$pfs.surv, 2, function(x){
      expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
    
    pd.expected.surv.ci <- quantile(pd.surv, c(0.025, 0.975), na.rm = TRUE)
    
    os.expected.surv.ci <- quantile(apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
    
    ## Calculated mean time in each state to end of trial followup
    
    pd.obs.surv <- mean(apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times[times < f.up], surv = x)}) - apply(partsa$pfs.surv, 2, function(x){
        expected_surv(time = times[times < f.up], surv = x)}))
    
    pfs.obs.surv <- mean(apply(partsa$pfs.surv, 2, function(x){
      expected_surv(time = times[times < f.up], surv = x)}))
    
    os.obs.surv <- mean(apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times[times < f.up], surv = x)}))
    
       res1 <- list(pfs = partsa$pfs.expected.surv, 
                 pfs_ci = pfs.expected.surv.ci,
                 pd = pd.expected.surv,
                 pd_ci = pd.expected.surv.ci,
                 os = partsa$os.expected.surv, 
                 os_ci = os.expected.surv.ci,
                 chosen_dist = partsa$chosen_models,
                 pfs_dist = fit_pfs$bestfit,
                 os_dist = fit_os_rel$bestfit,
                 pfs_events = pfs.events,
                 os_events = os.events,
                 f_up = f.up,
                 pfs_obs = pfs.obs.surv,
                 pd_obs = pd.obs.surv,
                 os_obs = os.obs.surv,
                 pfsdet = parts1$pfs.expected.surv, 
                 osdet = parts1$os.expected.surv
                 )
    }
      
}
    return(res1)          
  }, simplify = FALSE), c("aicc", "bicc")))
}

fit_function_nobhaz <- function(data) {
  
  fit_pfs  <- list(model.objects = fit.models1(
    formula = Surv(time_pfs2, status_pfs2) ~ 1,
    distr = distribs,
    data = data))
  
  fit_os  <- list(model.objects = fit.models1(
    formula = Surv(time_os2, status_os2) ~ 1,
    data = data,
    distr = distribs))
  
  pfs.events <- sum(data$status_pfs2, na.rm = TRUE)
  os.events <- sum(data$status_os2, na.rm = TRUE)
  f.up = max(data$time_pfs2, na.rm = TRUE)
  
  
  list(setNames(sapply(seq_along(c("aicc", "bicc")), function(i){
    
    fit_pfs$bestfit <- tryCatch(
      names(fit_pfs$model.objects$models[which.min(fit_pfs[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    fit_os$bestfit <- tryCatch(
      names(fit_os$model.objects$models[which.min(fit_os[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    
    if (is_na(fit_os$bestfit)) {
      res1 <- list(pfs = NA, 
                   pfs_ci = NA,
                   pd = NA,
                   pd_ci = NA,
                   os = NA, 
                   os_ci = NA,
                   chosen_dist = NA,
                   pfs_dist = fit_pfs$bestfit,
                   os_dist = fit_os_rel$bestfit,
                   pfs_events = pfs.events,
                   os_events = os.events,
                   f_up = f.up,
                   pfs_obs = NA,
                   pd_obs = NA,
                   os_obs = NA)
      
    } else {
      
      parts1 <- tryCatch(partsurv_bhaz(
        pfs_survHE = fit_pfs,
        os_survHE = fit_os_rel,
        choose_PFS = fit_pfs$bestfit,
        choose_OS = fit_os_rel$bestfit,
        time = times,
        PA = FALSE,
        #seed = seed,
        v_names_states = v_names_states),
        error = function(e) list(NA), warning = function(w) list(NA))
    
    partsa <- tryCatch(partsurv_nobhaz(
      pfs_survHE = fit_pfs,
      os_survHE = fit_os,
      choose_PFS = fit_pfs$bestfit,
      choose_OS = fit_os$bestfit,
      time = times,
      PA = TRUE,
      n_sim = 1000,
      #seed = seed,
      v_names_states = v_names_states),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    if (is_na(partsa)){
      
      res1 <- list(pfs = NA, 
                   pfs_ci = NA,
                   pd = NA,
                   pd_ci = NA,
                   os = NA, 
                   os_ci = NA,
                   chosen_dist = NA,
                   pfs_dist = fit_pfs$bestfit,
                   os_dist = fit_os$bestfit,
                   pfs_events = pfs.events,
                   os_events = os.events,
                   f_up = f.up,
                   pfs_obs = NA,
                   pd_obs = NA,
                   os_obs = NA)
      
    } else {
    
    ## Calculate mean time in each state and related CI
    pd.surv <- apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times, surv = x)}) - apply(partsa$pfs.surv, 2, function(x){
        expected_surv(time = times, surv = x)})
    
    pd.expected.surv <- mean(pd.surv)
    
    pfs.expected.surv.ci <- quantile(apply(partsa$pfs.surv, 2, function(x){
      expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
    
    pd.expected.surv.ci <- quantile(pd.surv, c(0.025, 0.975), na.rm = TRUE)
    
    os.expected.surv.ci <- quantile(apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
    
    ## Calculated mean time in each state to end of trial followup
    
    pd.obs.surv <- mean(apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times[times < f.up], surv = x)}) - apply(partsa$pfs.surv, 2, function(x){
        expected_surv(time = times[times < f.up], surv = x)}))
    
    pfs.obs.surv <- mean(apply(partsa$pfs.surv, 2, function(x){
      expected_surv(time = times[times < f.up], surv = x)}))
    
    os.obs.surv <- mean(apply(partsa$os.surv, 2, function(x){
      expected_surv(time = times[times < f.up], surv = x)}))
    
    res1 <- list(pfs = partsa$pfs.expected.surv, 
                 pfs_ci = pfs.expected.surv.ci,
                 pd = pd.expected.surv,
                 pd_ci = pd.expected.surv.ci,
                 os = partsa$os.expected.surv, 
                 os_ci = os.expected.surv.ci,
                 chosen_dist = partsa$chosen_models,
                 pfs_dist = fit_pfs$bestfit,
                 os_dist = fit_os$bestfit,
                 pfs_events = pfs.events,
                 os_events = os.events,
                 f_up = f.up,
                 pfs_obs = pfs.obs.surv,
                 pd_obs = pd.obs.surv,
                 os_obs = os.obs.surv
    )
    }
    }
    return(res1)          
  }, simplify = FALSE), c("aicc", "bicc")))
}

fit_function_part <- function(data, PA, n_sim) {
  
  fit_pfs  <- list(model.objects = fit.models1(
    formula = Surv(time_pfs2, status_pfs2) ~ 1,
    data = data,
    distr = distribs,
    bhazard = NULL
  ))
  
  bhaz <- tryCatch(flexsurv::hgompertz(data$time_os2, shape = mort_reg$coefficients[2], rate = mort_reg$coefficients[1]),
                   error = function(e) list(NA), warning = function(w) list(NA))
  
  fit_os_rel  <- list(model.objects = fit.models1(
    formula = Surv(time_os2, status_os2) ~ 1,
    data = data,    
    distr = distribs,
    bhazard = bhaz))
  
    fit_os  <- list(model.objects = fit.models1(
    formula = Surv(time_os2, status_os2) ~ 1,
    data = data,
    distr = distribs))
  
  pfs.events <- sum(data$status_pfs2, na.rm = TRUE)
  os.events <- sum(data$status_os2, na.rm = TRUE)
  f.up = max(data$time_pfs2, na.rm = TRUE)
  
  
  list(setNames(sapply(seq_along(c("aicc", "bicc")), function(i){
    
    fit_pfs$bestfit <- tryCatch(
      names(fit_pfs$model.objects$models[which.min(fit_pfs[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    fit_os_rel$bestfit <- tryCatch(
      names(fit_os_rel$model.objects$models[which.min(fit_os_rel[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    fit_os$bestfit <- tryCatch(
      names(fit_os$model.objects$models[which.min(fit_os[["model.objects"]][["model.fitting"]][[i]])]),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    
    partsa <- tryCatch(partsurv_bhaz(
      pfs_survHE = fit_pfs,
      os_survHE = fit_os_rel,
      choose_PFS = fit_pfs$bestfit,
      choose_OS = fit_os_rel$bestfit,
      time = times,
      PA = PA,
      n_sim = n_sim,
      #seed = seed,
      v_names_states = v_names_states),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    partsa_s <- tryCatch(partsurv_nobhaz(
      pfs_survHE = fit_pfs,
      os_survHE = fit_os,
      choose_PFS = fit_pfs$bestfit,
      choose_OS = fit_os$bestfit,
      time = times,
      PA = PA,
      n_sim = n_sim,
      #seed = seed,
      v_names_states = v_names_states),
      error = function(e) list(NA), warning = function(w) list(NA))
    
    if (is_na(partsa)) {
      pd.surv <-  NA
      pd.expected.surv <- NA
      pd.expected.surv.ci <- NA
      pfs.expected.surv <- NA
      pfs.expected.surv.ci <- NA
      os.expected.surv <- NA
      os.expected.surv.ci <-  NA
      chosen_dist <- NA
      
    } else {
      
      ## Calculate mean time in each state and related CI from simplified analysis
      chosen_dist <- partsa$chosen_models
      pfs.expected.surv <- partsa$pfs.expected.surv  
      os.expected.surv <- partsa$os.expected.surv
      
      if (PA == TRUE) {
        
        pd.surv <- apply(partsa$os.surv, 2, function(x){
          expected_surv(time = times, surv = x)}) - apply(partsa$pfs.surv, 2, function(x){
            expected_surv(time = times, surv = x)})
        
        pd.expected.surv <- mean(pd.surv)
        
        pfs.expected.surv.ci <- quantile(apply(partsa$pfs.surv, 2, function(x){
          expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
        
        pd.expected.surv.ci <- quantile(pd.surv, c(0.025, 0.975), na.rm = TRUE)
        
        os.expected.surv.ci <- quantile(apply(partsa$os.surv, 2, function(x){
          expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
        
      } else {
        
        pd.surv <-  NA
        
        pd.expected.surv <- partsa$os.expected.surv - partsa$pfs.expected.surv
        
        pfs.expected.surv.ci <-  NA
        
        pd.expected.surv.ci <- NA
        
        os.expected.surv.ci <- NA
        
      }
    }
    
  if (is_na(partsa_s)) {
      pd.surv_s <-  NA
      pd.expected.surv_s <- NA
      pd.expected.surv.ci_s <- NA
      pfs.expected.surv_s <- NA
      pfs.expected.surv.ci_s <- NA
      os.expected.surv_s <- NA
      os.expected.surv.ci_s <-  NA
      chosen_dist_s <- NA
      
    } else {  
    
        ## Calculate mean time in each state and related CI from simplified analysis
      chosen_dist_s <- partsa_s$chosen_models
      pfs.expected.surv_s <- partsa_s$pfs.expected.surv  
      os.expected.surv_s <- partsa_s$os.expected.surv
      
      if (PA == TRUE) {
        
      pd.surv_s <- apply(partsa_s$os.surv, 2, function(x){
            expected_surv(time = times, surv = x)}) - apply(partsa_s$pfs.surv, 2, function(x){
            expected_surv(time = times, surv = x)})
      
      pd.expected.surv_s <- mean(pd.surv_s)
      
      pfs.expected.surv.ci_s <- quantile(apply(partsa_s$pfs.surv, 2, function(x){
        expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
      
      pd.expected.surv.ci_s <- quantile(pd.surv_s, c(0.025, 0.975), na.rm = TRUE)
      
      os.expected.surv.ci_s <- quantile(apply(partsa_s$os.surv, 2, function(x){
        expected_surv(time = times, surv = x)}), c(0.025, 0.975), na.rm = TRUE)
        
      } else {
     
      pd.surv_s <-  NA
      
      pd.expected.surv_s <- partsa_s$os.expected.surv - partsa_s$pfs.expected.surv
        
      pfs.expected.surv.ci_s <-  NA
     
      pd.expected.surv.ci_s <- NA
        
      os.expected.surv.ci_s <- NA
        
      }
    }
      

      ## Combine

        res <- list(pfs = pfs.expected.surv, 
                    pfs_ci = pfs.expected.surv.ci, 
                    pd = pd.expected.surv,
                    pd_ci = pd.expected.surv.ci, 
                    os = os.expected.surv, 
                    os_ci = os.expected.surv.ci,
                    chosen_dist = chosen_dist,
                    pfs_s = pfs.expected.surv_s, 
                    pfs_ci_s = pfs.expected.surv.ci_s,
                    pd_s = pd.expected.surv_s,
                    pd_ci_s = pd.expected.surv.ci_s,
                    os_s = os.expected.surv_s, 
                    os_ci_s = os.expected.surv.ci_s,
                    chosen_dist_s = chosen_dist_s,                 
                    pfs_events = pfs.events,
                    os_events = os.events,
                    f_up = f.up,
                    pfs_dist = fit_pfs$bestfit,
                    os_dist = fit_os_rel$bestfit,
                    os_dist_s = fit_os$bestfit)
    
    # ## Calculated mean time in each state to end of trial followup
    # 
    # pd.obs.surv <- mean(apply(partsa$os.surv, 2, function(x){
    #   expected_surv(time = times[times < f.up], surv = x)}) - apply(partsa$pfs.surv, 2, function(x){
    #     expected_surv(time = times[times < f.up], surv = x)}))
    # 
    # pfs.obs.surv <- mean(apply(partsa$pfs.surv, 2, function(x){
    #   expected_surv(time = times[times < f.up], surv = x)}))
    # 
    # os.obs.surv <- mean(apply(partsa$os.surv, 2, function(x){
    #   expected_surv(time = times[times < f.up], surv = x)}))
    
    return(res)          
  }, simplify = FALSE), c("aicc", "bicc")))
}

fit_function_des <- function(data, B1, B2) {
  
  fit_HS <- list(models = map(seq_along(distribs),function(i){
    tryCatch(flexsurvreg(Surv(time2_2, status2_2) ~ 1, dist = distribs[[i]], data = data),
             error = function(e) NA, warning = function(w) NA)}))
  
  fit_HS$model.fitting = list(
    aic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HS$models[i]), NA, fit_HS$models[[i]]$AIC)})),
    bic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HS$models[i]), NA, glance(fit_HS$models[[i]])$BIC)})),
    aicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HS$models[i]), NA, 
             - 2 * fit_HS$models[[i]]$loglik + fit_HS$models[[i]]$npars * 2 * fit_HS$models[[i]]$N / 
               (fit_HS$models[[i]]$N - fit_HS$models[[i]]$npars - 1))})),        
    bicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HS$models[i]), NA,
             - 2 * fit_HS$models[[i]]$loglik + fit_HS$models[[i]]$npars * log(fit_HS$models[[i]]$N) * fit_HS$models[[i]]$N /
               (fit_HS$models[[i]]$N - fit_HS$models[[i]]$npars - 1))}))
  ) 
  
  
  bhazHD <- flexsurv::hgompertz(data$time2_3, shape = mort_reg$coefficients[2], rate = mort_reg$coefficients[1])  
  
  
  fit_HD <- list(models = map(seq_along(distribs),function(i){
    tryCatch(flexsurvreg(Surv(time2_3, status2_3) ~ 1, dist = distribs[[i]], data = data),
             error = function(e) NA, warning = function(w) NA)}))
  
  fit_HD$model.fitting = list(
    aic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HD$models[i]), NA, fit_HD$models[[i]]$AIC)})),
    bic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HD$models[i]), NA, glance(fit_HD$models[[i]])$BIC)})),
    aicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HD$models[i]), NA, 
             - 2 * fit_HD$models[[i]]$loglik + fit_HD$models[[i]]$npars * 2 * fit_HD$models[[i]]$N / 
               (fit_HD$models[[i]]$N - fit_HD$models[[i]]$npars - 1))})),     
    bicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HD$models[i]), NA,
             - 2 * fit_HD$models[[i]]$loglik + fit_HD$models[[i]]$npars * log(fit_HD$models[[i]]$N) * fit_HD$models[[i]]$N /
               (fit_HD$models[[i]]$N - fit_HD$models[[i]]$npars - 1))}))
  )
  
  
  fit_HDbhaz <- list(models = map(seq_along(distribs),function(i){
    tryCatch(flexsurvreg(Surv(time2_3, status2_3) ~ 1, dist = distribs[[i]], data = data, bhazard = bhazHD),
             error = function(e) NA, warning = function(w) NA)}))
  
  fit_HDbhaz$model.fitting = list(
    aic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HDbhaz$models[i]), NA, fit_HDbhaz$models[[i]]$AIC)})),
    bic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HDbhaz$models[i]), NA, glance(fit_HDbhaz$models[[i]])$BIC)})),
    aicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HDbhaz$models[i]), NA, 
             - 2 * fit_HDbhaz$models[[i]]$loglik + fit_HDbhaz$models[[i]]$npars * 2 * fit_HDbhaz$models[[i]]$N / 
               (fit_HDbhaz$models[[i]]$N - fit_HDbhaz$models[[i]]$npars - 1))})),     
    bicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_HDbhaz$models[i]), NA,
             - 2 * fit_HDbhaz$models[[i]]$loglik + fit_HDbhaz$models[[i]]$npars * log(fit_HDbhaz$models[[i]]$N) * fit_HDbhaz$models[[i]]$N /
               (fit_HDbhaz$models[[i]]$N - fit_HDbhaz$models[[i]]$npars - 1))}))
  )
  
  
  
  bhazSD <- flexsurv::hgompertz(data$time2_4, shape = mort_reg$coefficients[2], rate = mort_reg$coefficients[1])  
  
  fit_SD <- list(models = map(seq_along(distribs),function(i){
    tryCatch(flexsurvreg(Surv(time2_4, status2_4) ~ 1, dist = distribs[[i]], bhazard = bhazSD, data = data),
             error = function(e) NA, warning = function(w) NA)}))
  
  fit_SD$model.fitting <- list(
    aic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_SD$models[i]), NA, fit_SD$models[[i]]$AIC)})),
    bic = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_SD$models[i]), NA, glance(fit_SD$models[[i]])$BIC)})),
    aicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_SD$models[i]), NA, 
             - 2 * fit_SD$models[[i]]$loglik + fit_SD$models[[i]]$npars * 2 * fit_SD$models[[i]]$N / 
               (fit_SD$models[[i]]$N - fit_SD$models[[i]]$npars - 1))})),     
    bicc = c(plyr::laply(seq_along(distribs), function(i){
      ifelse(is.na(fit_SD$models[i]), NA,
             - 2 * fit_SD$models[[i]]$loglik + fit_SD$models[[i]]$npars * log(fit_SD$models[[i]]$N) * fit_SD$models[[i]]$N /
               (fit_SD$models[[i]]$N - fit_SD$models[[i]]$npars - 1))}))
  )  
  
  os.events <- sum(data$status2_3, data$status2_4, na.rm = TRUE)
  
  list(setNames(sapply(seq_along(c("aicc", "bicc")), function(i){
    
    # Extract best fitting models
    
    bestfit_HS_dist <- ifelse(is.na(which.min(fit_HS$model.fitting[[i]])), NA,
                              distribs[[which.min(fit_HS$model.fitting[[i]])]])
    
    bestfit_HD_dist <- ifelse(is.na(which.min(fit_HD$model.fitting[[i]])), NA,
                              distribs[[which.min(fit_HD$model.fitting[[i]])]])
    
    bestfit_HDbhaz_dist <- ifelse(is.na(which.min(fit_HDbhaz$model.fitting[[i]])), NA,
                                  distribs[[which.min(fit_HDbhaz$model.fitting[[i]])]])
    
    bestfit_SD_dist <- ifelse(is.na(which.min(fit_SD$model.fitting[[i]])), NA,
                              distribs[[which.min(fit_SD$model.fitting[[i]])]])
    
    # Extract best fitting models
    bestfit_HS <- tryCatch(fit_HS$models[[which.min(fit_HS$model.fitting[[i]])]],
                           error = function(e) NA, warning = function(w) NA)
    
    bestfit_HD <- tryCatch(fit_HD$models[[which.min(fit_HD$model.fitting[[i]])]],
                           error = function(e) NA, warning = function(w) NA)
    
    bestfit_HDbhaz <- tryCatch(fit_HDbhaz$models[[which.min(fit_HDbhaz$model.fitting[[i]])]],
                               error = function(e) NA, warning = function(w) NA)
    
    
    bestfit_SD <- tryCatch(fit_SD$models[[which.min(fit_SD$model.fitting[[i]])]],
                           error = function(e) NA, warning = function(w) NA)
    
   
     # Fit partly naive MSM
    
    list_mod <- list(bestfit_HS, bestfit_HD, genpop_reg, bestfit_SD)
    
    DES_los <-  tryCatch(totlos.simfs(list_mod, trans = m_P_diag_des, t = max(times), M=10000, B = B1, ci = TRUE, group=c(1,2,3,3)),
                         error = function(e) NA, warning = function(w) NA)
    
    DES_los <- as.data.frame(DES_los)

    
    #Fit MSM with GPM
    
    list_mod_bhaz <- list(bestfit_HS, genpop_reg, bestfit_HDbhaz, genpop_reg, bestfit_SD)
    
    DES_los_bhaz <-  tryCatch(totlos.simfs(list_mod_bhaz, trans = m_P_diag_desbhaz, t = max(times), M=10000, B = 100, ci = TRUE, group=c(1,2,3,3)),
                              error = function(e) NA, warning = function(w) NA)
    
    DES_los_bhaz <- as.data.frame(DES_los_bhaz)
    
    simp = 0
    
    if (is_na(DES_los_bhaz)) {   
      
      list_mod_simp <- list(bestfit_HS, genpop_reg, genpop_reg, bestfit_SD)
      
      DES_los_simp <-  tryCatch(totlos.simfs(list_mod_simp, trans = m_P_diag_des, t = max(times), M=10000, B = 100, ci = TRUE, group=c(1,2,3,3)),
                                error = function(e) NA, warning = function(w) NA)
      
      DES_los_bhaz <- as.data.frame(DES_los_simp)
      
      simp = 1
      
    }
    
    
   # Store results
    
    res  <- list(pfs = DES_los$est[[1]],
                 pfs_ci = c(DES_los$L[[1]], DES_los$U[[1]]),
                 pd = DES_los$est[[2]],
                 pd_ci = c(DES_los$L[[2]], DES_los$U[[2]]),
                 dead = DES_los$est[[3]],
                 dead_ci = c(DES_los$L[[3]], DES_los$U[[3]]),
                 os = max(times) - DES_los$est[[3]],
                 os_check = sum(DES_los$est[[1]], DES_los$est[[2]]),
                 os_ci = c(max(times) - DES_los$U[[3]], max(times) - DES_los$L[[3]]),
                 bestfit_HS_dist = bestfit_HS_dist,
                 bestfit_HD_dist = bestfit_HD_dist,
                 bestfit_HDbhaz_dist = bestfit_HDbhaz_dist,
                 bestfit_SD_dist = bestfit_SD_dist,
                 pfs_b = DES_los_bhaz$est[[1]],
                 pfs_ci_b = c(DES_los_bhaz$L[[1]], DES_los_bhaz$U[[1]]),
                 pd_b = DES_los_bhaz$est[[2]],
                 pd_ci_b = c(DES_los_bhaz$L[[2]], DES_los_bhaz$U[[2]]),
                 dead_b = DES_los_bhaz$est[[3]],
                 dead_ci_b = c(DES_los_bhaz$L[[3]], DES_los_bhaz$U[[3]]),
                 os_b = max(times) - DES_los_bhaz$est[[3]],
                 os_check_b = sum(DES_los_bhaz$est[[1]], DES_los_bhaz$est[[2]]),
                 os_ci_b = c(max(times) - DES_los_bhaz$U[[3]], max(times) - DES_los_bhaz$L[[3]]),
                 os_events = os.events,
                simp = simp
                 )
    
    return(res)
  }, simplify = FALSE), c("aicc", "bicc")))
}

