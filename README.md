# PrevalenceMapping
The package performs Bayesian geostatistical spatial mapping of disease prevalence using prevalence survey data and environmental covariate raster files. Using an integrated nested Laplace approximation (INLA) framework, the outputs include mapped prevalence in the form of raster files, and measures of in-sample-fit and cross-validation.

#### Install the package direct from github using:

``` r
devtools::install_github(repo = "suyunkang/PrevalenceMapping")
```

#### Install with the example and mock data (Longer download time):

``` r
devtools::install_github(repo = "suyunkang/PrevalenceMapping", ref = "Documentation")
```


----------

## Worked example

Below is an example using mock data and covariates which are included in the _"Documentation"_ branch of the package. The full script may be found on GitHub at  `inst/doc/Spatial_model_mock_data.R`. Here I present a brief bare-bones example

 1. Load a response dataset
 2. Create a mesh for it
 3. Load in covariate raster files
 4. Fit a model with parameter selection (minimum DIC wins)
 5. Perform cross-validation
 6. Visualise the cross-validation results 


First, specify working directory and load packages. We'll make use of a few `raster` functions at the top level, so it's good to load it as we go as well. You may choose to use your own working directory, we're going to store the model plots and output in a new folder `Outputs`.
``` r
dir.create("Outputs", showWarnings = FALSE)
setwd("Outputs")
library(PrevalenceMapping)
library(raster)
```

Load the response data and do any pre-procsessing.
``` r
Resp = read.csv(file.path(.libPaths(), "PrevalenceMapping/doc/ResponseData/MOZ_mock_data.csv"))
Resp$examined = round(Resp$tested, 0)
Resp$n_positive = round(Resp$positive, 0)
```

The mock data has columns `latitude`, `longitude`, `tested` and `positive`, eg, the first five rows are

| latitude | longitude | tested | positive|
|-|-|-|-|
| -14.24 | 39.82 | 37.58 | 15.37 |
| -24.93 | 31.82 | 6.11  | 0.00  |
| -11.00 | 35.42 | 5.59  | 2.93  |
| -15.95 | 33.48 | 20.54 | 11.56 |
| -25.52 | 33.11 | 13.56 | 0.00  |
| ... |  |  |  |

With the response we make a mesh and build an SPDE model. The parameters in `control`, `"convex"`, `"concave"`, `"max_edge"`, and `"cutoff"` may require some tuning, it's common to try a broad range of meshes before moving on to prediction later.

``` r
out_mesh = makeMeshSPDE(response = Resp, 
                        control = list(convex = -0.1,
                                       concave = -0.1,
                                       max_edge = c(0.2, 0.4),
                                       cutoff = 0.3),
                        plot_mesh = TRUE)
```

We have our dataset `Resp` and our mesh `out_mesh`,  now stack the environmental covariate rasters using the example covariates `.tif` files. We set any `NA`'s to be `-9999`. The covariates in this data are useful malaria prevalence predictors like `distance to water`,  and several `land surface temperature` among others.
``` r
folder = file.path(.libPaths(), "PrevalenceMapping/doc/Covariates")
lsf = grep("*.tif$", list.files(folder, full.names = TRUE), value = T)
covariate_stack = stack(lsf)
NAvalue(covariate_stack) = -9999
```

With the response, the mesh and the covariates loaded we fit and, with parameter selection, re-fit the model, choosing the "best" model as the one with the smallest deviance information criterion (DIC). This function produces `"Fixed_effects.txt"` -- a table of coefficients for the final set of covariates.

``` r
mod_final = findModelWithSmallestDIC(response = Resp, 
                                     raster_stack = covariate_stack,
                                     A_mat = out_mesh$A_mat,
                                     spde = out_mesh$spde,
                                     family = "binomial")
```

Using the model we predict prevalence on a grid/raster and compute posterior samples. Use the `nsamp` parameter to choose number of desired realisations of the posterior distribution and `write_posterior` to choose whether to write the posterior samples into `.tif` files.

``` r
nsamp = 1000
mod_final = read.table("Fixed_effects.txt")
predict = predictionOnAGrid(response = Resp,
                            finalmodel = mod_final,
                            A_mat = out_mesh$A_mat,
                            spde = out_mesh$spde,
                            mesh = out_mesh$mesh,
                            family = "binomial",
                            raster_stack = covariate_stack,
                            nsamp = nsamp,
                            int.strategy = "eb",
                            write_posterior = TRUE)
```

Cross-validation statistics are simple to produce, use parameter `n_reps` to choose number of desired replicates and inspect the output `validation_result`.
``` r
n_reps = 100
test_pct = c(0.1, 0.2, 0.3, 0.4)
validation_result = vector("list", 4)
names(validation_result) = paste0(test_pct*100, "pct")
for (k in 1:length(test_pct)){
  fn <- paste0("Cross_validation_", test_pct[k]*100, "pct.txt")
  if (file.exists(fn)) file.remove(fn)
  validation_result[[k]] = crossValidation(response = Resp, 
                                           finalmodel = mod_final, 
                                           A_mat = out_mesh$A_mat, 
                                           spde = out_mesh$spde, 
                                           family = "binomial", 
                                           raster_stack = covariate_stack, 
                                           int.strategy = "eb", 
                                           n_reps = n_reps, 
                                           pct_out = test_pct[k])
}
validation_result
```

Finally, we may also choose to visualise the validation results.
``` r
# pdf("Cross_validation_result.pdf")
par(mfrow = c(2,1), mar = c(3.5,4,3,1), mgp = c(2,1,0))
## Pearson correlation 
DF = data.frame(Cor10 = validation_result$`10pct`[,1],
                Cor20 = validation_result$`20pct`[,1],
                Cor30 = validation_result$`30pct`[,1],
                Cor40 = validation_result$`40pct`[,1])
boxplot(DF,
        names = c("10%", "20%", "30%", "40%"),
        ylab = expression(rho),
        xlab = "Percentage of validation set",
        col = 2:6,
        main = "(A) Pearson correlation coefficient",
        cex.main = 1)
## Cross-validated R-squared
DF2 = data.frame(Cor10 = validation_result$`10pct`[,2],
                 Cor20 = validation_result$`20pct`[,2],
                 Cor30 = validation_result$`30pct`[,2],
                 Cor40 = validation_result$`40pct`[,2])
boxplot(DF2,
        names = c("10%", "20%", "30%", "40%"), 
        ylab = expression(R^2), 
        xlab = "Percentage of validation set", 
        col = 2:6,
        main = "(B) Cross-validated R-squared",
        cex.main = 1)
# dev.off()
```

