





.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This package relies on R INLA, you should install it manually by running: \n\n\tinstall.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/stable\"), dep=TRUE) \n\nFull instructions here:\n\n\thttp://www.r-inla.org/download \n\n")
}

# .onLoad <- function(libname, pkgname) {
#   tryCatch(library(INLA), error = function(e) {
#     # packageStartupMessage("Package INLA not found, attempting install")
#     install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#   })
# }


#' @importFrom utils flush.console write.table install.packages
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline points
#' @importFrom stats as.formula cor formula lm plogis sd


#### FUNCTION #1
#' @title Make a mesh, build an spde model, and a projector matrix
#' @description Construct a nonconvex boundary for a set of points. Create a triangle mesh based on initial point locations, specified or automatic boundaries, and mesh quality parameters. Construct observation/prediction weight matrices for models. Create an inla.spde2 model object for a Matern model.
#' @param response A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).
#' @param control Parameters controlling inla.nonconvex.hull and inla.mesh.2d, including convex, concave, max_edge, and cutoff.
#' @param plot_mesh A logical argument indicating if mesh is to be plotted, Default: FALSE.
#' @return OUTPUT_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
makeMeshSPDE <- function(response, control = list(convex = NULL, concave = NULL, max_edge = NULL, cutoff = NULL), plot_mesh = FALSE) {
  coords = cbind(response$longitude, response$latitude)
  boundary = INLA::inla.nonconvex.hull(points = coords, convex = control$convex, concave = control$concave)

  # Make the mesh
  mesh = INLA::inla.mesh.2d(loc = coords, boundary = boundary, max.edge = control$max_edge, cutoff =  control$cutoff)

  # Mapping between meshes and continuous space
  A = INLA::inla.spde.make.A(mesh = mesh, loc = coords)

  # SPDE model construction
  spde = INLA::inla.spde2.matern(mesh, alpha = 2)

  if(plot_mesh) {
    plot(mesh)
    points(coords, pch = 20, col = "red", cex = 0.2)
  }

  return(list(spde = spde, A_mat = A, n_points = mesh$n, mesh = mesh))
}



#### FUNCTION #2
#' @title Calculate variance inflation factor (VIF) and remove variables with VIF > threshold value
#' @description The function checks for collinearity between variables and performs stepwise VIF selection. Source: https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/.
#' @author Marcus W. Beck
#' @param response A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).
#' @param raster_stack A collection of RasterLayer objects with the same spatial extent and resolution (These are the environmental covariates).
#' @param thresh The threshold value used for retaining variables, Default: 10.
#' @param trace A logical argument indicating if text output is returned as the stepwise selection progresses, Default: T.
#' @param ... Additional arguments passed to `lm'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
calculateVIFandRemoveVariables <- function(response, raster_stack, thresh = 10, trace = T,...){

  coords = cbind(response$longitude, response$latitude)
  covariate_DF = data.frame(raster::extract(raster_stack, coords))

  if(class(covariate_DF) != 'data.frame') covariate_DF <- data.frame(covariate_DF)
  # get initial vif value for all comparisons of variables
  vif_init <- NULL
  var_names <- names(covariate_DF)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init <- rbind(vif_init, c(val, 1/(1 - summary(lm(form_in, data = covariate_DF, ...))$r.squared)))
  }
  vif_max <- max(as.numeric(vif_init[ ,2]))
  if(vif_max < thresh){
    if(trace == T){ # print output of each iteration
      prmatrix(vif_init, collab = c('var', 'vif'), rowlab=rep('', nrow(vif_init)), quote = F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh, ', max VIF ', round(vif_max,2), sep = ''), '\n\n')
    }
    return(var_names)
  }
  else{
    in_dat <- covariate_DF
    # backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      vif_vals <- NULL
      var_names <- names(in_dat)
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add <- 1/(1 - summary(lm(form_in, data = in_dat, ...))$r.squared)
        vif_vals <- rbind(vif_vals, c(val, vif_add))
      }
      max_row <- which(vif_vals[ , 2] == max(as.numeric(vif_vals[ , 2])))[1]
      vif_max <- as.numeric(vif_vals[max_row, 2])
      if(vif_max < thresh) break
      if(trace == T){ # print output of each iteration
        prmatrix(vif_vals, collab = c('var', 'vif'), rowlab = rep('', nrow(vif_vals)), quote = F)
        cat('\n')
        cat('removed: ', vif_vals[max_row,1], vif_max,'\n\n')
        flush.console()
      }
      in_dat <- in_dat[ , !names(in_dat) %in% vif_vals[max_row, 1]]
    }
    return(names(in_dat))
  }
}



#### FUNCTION #3
#' @title Find the best spatial model, i.e. the model with the smallest DIC
#' @description The function performs forward and backward elimination of variables in order to find the spatial model with a set of covariates that result in the smallest DIC.
#' @param response A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).
#' @param raster_stack A collection of RasterLayer objects with the same spatial extent and resolution.
#' @param A_mat An observation/prediction weight matrix returned from makeMeshSPDE.
#' @param spde An inla.spde2 model object for a Matern model returned from makeMeshSPDE.
#' @param family A string indicating the likelihood family, Default: `binomial'.
#' @param save_output A logical argument indicating if output to be saved to a csv file, Default: TRUE.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
findModelWithSmallestDIC <- function(response, raster_stack, A_mat, spde, family = "binomial", save_output = TRUE){

  coords = cbind(response$longitude, response$latitude)
  covariate_DF = data.frame(raster::extract(raster_stack, coords))
  keep = calculateVIFandRemoveVariables(response, raster_stack, thresh = 10)
  covars_keep = covariate_DF[ , keep]

  # Stacking data
  data.stk = INLA::inla.stack(data = list(n_positive = as.vector(response$n_positive),
                                          examined = as.vector(response$examined)),
                              A = list(A_mat, 1),
                              effects=list(field = 1:spde$n.spde,
                                           list(intercept = rep(1, nrow(response)),
                                                covariate = covars_keep)))

  # Model fitting
  options(useFancyQuotes = FALSE)
  formula = as.formula(paste("n_positive ~ -1 + intercept + f(field, model = spde) +", paste(names(covars_keep), collapse="+")))

  result = INLA::inla(formula, family = family,
                      Ntrials = INLA::inla.stack.data(data.stk)$examined,
                      data = INLA::inla.stack.data(data.stk), verbose=FALSE,
                      control.predictor = list(compute=TRUE,
                                               A = INLA::inla.stack.A(data.stk)),
                      control.fixed = list(expand.factor.strategy = 'inla'),
                      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                      control.inla = list(int.strategy = "eb"))

  # Assess the importance of the variable one at a time
  nm = noquote(names(covars_keep))
  DIC = rep(NA, length(nm)+1)
  DIC[1] = result$dic$dic

  nmFinal <- nm

  for (i in 1:length(nm)) {
    j = nm[i]

    covars = nmFinal[-which(nmFinal == j)]

    formula = as.formula(paste("n_positive ~ -1 + intercept + f(field, model = spde) +", paste(covars, collapse="+")))

    output = INLA::inla(formula, family = family,
                        Ntrials = INLA::inla.stack.data(data.stk)$examined,
                        data = INLA::inla.stack.data(data.stk), verbose = FALSE,
                        control.predictor = list(compute = TRUE,
                                                 A = INLA::inla.stack.A(data.stk)),
                        control.fixed = list(expand.factor.strategy = 'inla'),
                        control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                        control.inla = list(int.strategy = "eb"))

    if (output$dic$dic < min(DIC[1:i])){
      nmFinal = covars}

    DIC[i+1] = output$dic$dic

  }

  # Add variables back to the model one at a time
  outnm = nm[which(!nm %in% nmFinal)]

  if (length(outnm) > 0){
    nmFinal = nmFinal
    DIC2 = rep(NA, length(outnm)+1)
    DIC2[1] = min(DIC)

    for (i in 1:length(outnm)) {
      j = outnm[i]

      covars = c(nmFinal, j)

      formula = as.formula(paste("n_positive ~ -1 + intercept + f(field, model = spde) +", paste(covars, collapse="+")))

      output = INLA::inla(formula, family = family,
                          Ntrials = INLA::inla.stack.data(data.stk)$examined,
                          data = INLA::inla.stack.data(data.stk), verbose=FALSE,
                          control.predictor = list(compute = TRUE,
                                                   A = INLA::inla.stack.A(data.stk)),
                          control.fixed = list(expand.factor.strategy = 'inla'),
                          control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                          control.inla = list(int.strategy = "eb"))

      if (output$dic$dic < min(DIC2[1:i])){
        nmFinal = covars}

      DIC2[i+1] = output$dic$dic

    }}

  # nmFinal is the final set of cavariates with the smallest DIC
  covars = nmFinal

  formula = as.formula(paste("n_positive ~ -1 + intercept + f(field, model = spde) +", paste(covars,  collapse="+")))

  output = INLA::inla(formula, family = family,
                      Ntrials = INLA::inla.stack.data(data.stk)$examined,
                      data = INLA::inla.stack.data(data.stk), verbose = FALSE,
                      control.predictor = list(compute = TRUE,
                                               A = INLA::inla.stack.A(data.stk)),
                      control.fixed = list(expand.factor.strategy='inla'),
                      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                      control.inla = list(int.strategy = "eb"))
  if(save_output) {
    write.table(output$summary.fixed, "Fixed_effects.txt", row.names = TRUE, col.names = TRUE)
  }
  return(output$summary.fixed)

}



#### FUNCTION #4
#' @title Prediction of the response on a grid/raster
#' @description The function predicts the response on target locations where data are not observed using posterior distributions.
#' @param response A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).
#' @param finalmodel An object returned from the function findModelWithSmallestDIC.
#' @param A_mat An observation/prediction weight matrix returned from makeMeshSPDE.
#' @param spde An inla.spde2 model object for a Matern model returned from makeMeshSPDE,
#' @param mesh A triangle mesh created based on initial point locations and returned from makeMeshSPDE.
#' @param family A string indicating the likelihood family, Default: `binomial'.
#' @param raster_stack A collection of RasterLayer objects with the same spatial extent and resolution.
#' @param nsamp Number of samples to draw from an approximated posterior of a fitted model. Make nsamp >= 100 in order to compute mean, sd, IQR, and 95\% CI of posterior samples.
#' @param int.strategy Character. The integration strategy to use; one of `auto', `ccd', `grid', `eb' (empirical bayes), `user' or `user.std'.
#' @param write_posterior A logical argument indicating if posterior realizations are to be written into raster files, Default: TRUE.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
predictionOnAGrid <- function(response, finalmodel, A_mat, spde, mesh, family = "binomial", raster_stack, nsamp, int.strategy, write_posterior = TRUE){

  covar_names = rownames(finalmodel)[-1]
  coords = cbind(response$longitude, response$latitude)
  covariate_DF = data.frame(raster::extract(raster_stack, coords))
  covars_final = subset(covariate_DF, select = as.vector(covar_names))

  # Stacking data
  data.stk = INLA::inla.stack(data = list(n_positive = as.vector(response$n_positive),
                                          examined = as.vector(response$examined)),
                              A = list(A_mat, 1),
                              effects=list(field = 1:spde$n.spde,
                                           list(intercept = rep(1, nrow(response)),
                                                covariate = covars_final)))

  options(useFancyQuotes = FALSE)
  formula = as.formula(paste("n_positive ~ -1 + intercept + f(field, model=spde) +", paste(names(covars_final), collapse="+")))

  result = INLA::inla(formula, family = family,
                      Ntrials = INLA::inla.stack.data(data.stk)$examined,
                      data = INLA::inla.stack.data(data.stk), verbose = F,
                      control.predictor = list(compute = TRUE,
                                               A = INLA::inla.stack.A(data.stk)),
                      control.fixed = list(expand.factor.strategy = 'inla'),
                      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                      control.inla = list(int.strategy = int.strategy))


  # Locations for predictions
  mask <- raster_stack[[1]]
  NAvalue(mask) <- -9999
  pred_val <- raster::getValues(mask)
  w <- is.na(pred_val)
  index <- 1:length(w)
  index <- index[!w]
  pred_locs <- raster::xyFromCell(mask, 1:ncell(mask))
  pred_locs <- pred_locs[!w, ]
  colnames(pred_locs) <- c('longitude', 'latitude')
  locs_pred <- pred_locs

  A.pred <- INLA::inla.spde.make.A(mesh = mesh, loc = locs_pred)

  covariates_pred = data.frame(raster::extract(raster_stack, locs_pred))
  covar_pred = subset(covariates_pred, select = as.vector(covar_names))

  ## Draw from posterior samples
  nsamp = nsamp
  samp = INLA::inla.posterior.sample(nsamp, result)
  k = dim(covar_pred)[2] # number of final covariates
  test = raster_stack[[1]]
  pred = matrix(NA, nrow = dim(A.pred)[1], ncol = 100)

  for (i in 1:nsamp){

    field = samp[[i]]$latent[grep('field', rownames(samp[[i]]$latent)), ]
    intercept = samp[[i]]$latent[grep('intercept', rownames(samp[[i]]$latent)), ]
    beta = NULL
    for (j in 1:k){
      beta[j] = samp[[i]]$latent[grep(names(covar_pred)[j], rownames(samp[[i]]$latent)), ]
    }

    lp = as.matrix(intercept + c(as.matrix(covar_pred)%*%beta) + drop(A.pred%*%field))

    # Predicted values
    if(i <= 100) {
      pred[,i] = plogis(lp)
    }
    pred_prev = plogis(lp)

    ## Write the posterior surface into a raster file
    if(write_posterior) {
      pred_val[!w] <- round(pred_prev,3)
      out_mean = raster::setValues(test, pred_val)
      out_mean = raster::writeRaster(out_mean, paste0("Prevalence_posterior_realization_", i, ".tif"), overwrite = TRUE)
    }

  }


  ## Use 100 posterior samples and compute mean, median, CI, IQR etc.
  pred_mean = rowMeans(pred, na.rm = T)
  pred_sd = apply(pred, 1, sd, na.rm = T)
  pred_quants <- apply(pred, 1, quantile, probs = c(0.025, 0.25, 0.75, 0.975), na.rm = T)
  pred_IQR = pred_quants[3,] - pred_quants[2,]
  pred_95CI = pred_quants[4,] - pred_quants[1,]


  # Plot prevalence surfaces
  pdf("Prevalence_maps.pdf")

  # Write the output into a raster file
  pred_val[!w] <- round(pred_mean, 3)
  out_mean = raster::setValues(test, pred_val)
  out_mean = raster::writeRaster(out_mean, "Prevalence_mean.tif", overwrite = TRUE)
  plot(out_mean, main = "Prevalence mean")

  pred_val[!w] <- round(pred_sd, 3)
  out_sd = raster::setValues(test, pred_val)
  out_sd = raster::writeRaster(out_sd, "Prevalence_sd.tif", overwrite = TRUE)
  plot(out_sd, main = "Prevalence sd")

  pred_val[!w] <- round(pred_95CI, 3)
  out_CI = raster::setValues(test, pred_val)
  out_CI = raster::writeRaster(out_CI, "Prevalence_95CI.tif", overwrite = TRUE)
  plot(out_CI, main = "Prevalence 95%_CI")

  pred_val[!w] <- round(pred_IQR, 3)
  out_IQR = raster::setValues(test, pred_val)
  out_IQR = raster::writeRaster(out_IQR, "Prevalence_IQR.tif", overwrite = TRUE)
  plot(out_IQR, main = "Prevalence interquartile range")
  
  pred_val[!w] <- round(pred_quants[4,], 3)
  out_975 = raster::setValues(test, pred_val)
  out_975 = raster::writeRaster(out_975, "Prevalence_upperbound_975.tif", overwrite = TRUE)
  plot(out_975, main = "Prevalence (upper bound)")
  
  pred_val[!w] <- round(pred_quants[1,], 3)
  out_025 = raster::setValues(test, pred_val)
  out_025 = raster::writeRaster(out_025, "Prevalence_lowerbound_025.tif", overwrite = TRUE)
  plot(out_025, main = "Prevalence (lower bound)")

  dev.off()


  # Validation plot
  pdf("In_sample_fit.pdf")

  PR = response$n_positive / response$examined
  coords = cbind(response$longitude, response$latitude)
  pred = raster::extract(out_mean, coords)
  PR = PR[which(is.na(pred) == FALSE)]
  pred = pred[which(is.na(pred) == FALSE)]
  plot(PR, pred, main=paste("Cor = ", round(cor(PR, pred),2)), xlab = "Observed", ylab = "Predicted", cex = 0.8, mgp = c(2,1,0), xlim = c(0,1), ylim = c(0,1))
  abline(0, 1, col = "red")

  dev.off()

  return(invisible())

}




#### FUNCTION #5
#' @title Fit the final spatial model using a subset of data
#' @description Fit the final spatial model using a subset of data for the purpose of cross-validation.
#' @param response A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).
#' @param val_DAT Subset of data for cross-validation.
#' @param finalmodel An object returned from the function findModelWithSmallestDIC.
#' @param A_mat An observation/prediction weight matrix returned from makeMeshSPDE.
#' @param spde An inla.spde2 model object for a Matern model returned from makeMeshSPDE.
#' @param family A string indicating the likelihood family, Default: `binomial'.
#' @param raster_stack A collection of RasterLayer objects with the same spatial extent and resolution.
#' @param int.strategy Character. The integration strategy to use; one of `auto', `ccd', `grid', `eb' (empirical bayes), `user' or `user.std'.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname fitFinalModel
#'
#' @import raster
fitFinalModel = function(response, val_DAT, finalmodel, A_mat, spde, family = "binomial", raster_stack, int.strategy){

  aa = which(is.na(val_DAT$n_positive) == TRUE)
  covar_names = rownames(finalmodel)[-1]
  coords = cbind(val_DAT$longitude, val_DAT$latitude)
  covariate_DF = data.frame(raster::extract(raster_stack, coords))
  covars_final = subset(covariate_DF, select = as.vector(covar_names))

  # Stacking data
  data.stk = INLA::inla.stack(data = list(n_positive = as.vector(val_DAT$n_positive),
                                          examined = as.vector(val_DAT$examined)),
                              A = list(A_mat, 1),
                              effects=list(field = 1:spde$n.spde,
                                           list(intercept = rep(1, nrow(val_DAT)),
                                                covariate = covars_final)))

  options(useFancyQuotes = FALSE)
  formula = as.formula(paste("n_positive ~ -1 + intercept + f(field, model=spde) +", paste(names(covars_final), collapse="+")))

  result = INLA::inla(formula, family = family,
                      Ntrials = INLA::inla.stack.data(data.stk)$examined,
                      data = INLA::inla.stack.data(data.stk), verbose = F,
                      control.predictor = list(link=1, compute = TRUE, A = INLA::inla.stack.A(data.stk)),
                      control.fixed = list(expand.factor.strategy = 'inla'),
                      control.compute = list(dic = TRUE, waic = TRUE, config = TRUE),
                      control.inla = list(int.strategy = int.strategy))

  fitted = result$summary.fitted.values$mean[1:nrow(val_DAT)]
  obs = response$n_positive / response$examined
  out = c(round(cor(fitted, obs), 4), round(cor(fitted, obs)^2, 4))

  return(out)
}




#### FUNCTION #6
#' @title Performs cross-validation
#' @description Performs cross-validation for a subset of data to assess predictive performance of the final model.
#' @param response A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).
#' @param finalmodel An object returned from the function findModelWithSmallestDIC.
#' @param A_mat An observation/prediction weight matrix returned from makeMeshSPDE.
#' @param spde An inla.spde2 model object for a Matern model returned from makeMeshSPDE.
#' @param family A string indicating the likelihood family, Default: `binomial'.
#' @param raster_stack A collection of RasterLayer objects with the same spatial extent and resolution.
#' @param int.strategy Character. The integration strategy to use; one of `auto', `ccd', `grid', `eb' (empirical bayes), `user' or `user.std'.
#' @param n_reps Number of replicates of subsets of data for cross-validation, Default: `100'.
#' @param pct_out Percentage of data to be used as a test set.
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
crossValidation <- function(response, finalmodel, A_mat, spde, family = "binomial", raster_stack, int.strategy, n_reps = 100, pct_out){
  i = 1:n_reps

  cor_mat <- parallel::mclapply(X = i, FUN = function(...) {
    requireNamespace("INLA")
    requireNamespace("raster")
    aa = sample(1:nrow(response), size = floor(nrow(response)*pct_out), replace = FALSE)
    DAT = response
    DAT$n_positive[aa] = NA
    fitFinalModel(response = response, val_DAT = DAT, finalmodel = finalmodel, A_mat = A_mat, spde = spde, family = family, raster_stack = raster_stack, int.strategy = int.strategy)
  }, mc.cores = ifelse(tolower(Sys.info()[["sysname"]] == "windows"), 1, parallel::detectCores()))
  
  cor_mat <- t(as.data.frame(cor_mat))
  colnames(cor_mat) = c("Correlation mean", "Cross-validated R-squared")
  rownames(cor_mat) <- NULL

  return(cor_mat)

}




