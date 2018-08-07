

#################################################################
#### Prevalence mapping using a mock dataset : spatial model ####
#################################################################

rm(list = ls())

#Specify working directory
workd = ""  
setwd(workd)

## Install "PrevalenceMapping"
# devtools::install_github(repo = "suyunkang/PrevalenceMapping")
devtools::install_github(repo = "suyunkang/PrevalenceMapping", ref = "Documentation")

library(PrevalenceMapping)
library(raster)



##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##

## Response data
Resp = read.csv(file.path(.libPaths(), "PrevalenceMapping/doc/ResponseData/MOZ_mock_data.csv"))
Resp$examined = round(Resp$tested, 0)
Resp$n_positive = round(Resp$positive, 0)
Resp$longitude = Resp$longitude
Resp$latitude = Resp$latitude

## Stack covariates (a mix of static and dynamic covariates)
folder = file.path(.libPaths(), "PrevalenceMapping/doc/Covariates")
lsf = grep("*.tif$", list.files(folder, full.names = TRUE), value = T)
covariate_stack = stack(lsf)
NAvalue(covariate_stack) = -9999


## Store results here
setwd(paste0(workd,"/Outputs"))


##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##~~##

## Make a mesh and build an SPDE model
out_mesh = makeMeshSPDE(response = Resp, control = list(convex = -0.1, concave = -0.1, max_edge = c(0.2, 0.4), cutoff = 0.3), plot_mesh = TRUE)


## Find the best model with the smallest DIC (This function produces "Fixed_effects.txt")
mod_final = findModelWithSmallestDIC(response = Resp, raster_stack = covariate_stack, A_mat = out_mesh$A_mat, spde = out_mesh$spde, family = "binomial")


## Predict prevalence on a grid/raster and compute posterior samples 
## Use "nsamp" to choose number of desired realizations of the posterior distribution
nsamp = 100
mod_final = read.table("Fixed_effects.txt")
predict = predictionOnAGrid(response = Resp, finalmodel = mod_final, A_mat = out_mesh$A_mat, spde = out_mesh$spde, mesh = out_mesh$mesh, family = "binomial", raster_stack = covariate_stack, nsamp = nsamp, int.strategy = "eb", write_posterior = TRUE)


## Cross-validation statistics
## Use "n_reps" to choose number of desired replicates
n_reps = 100
test_pct = c(0.1, 0.2, 0.3, 0.4)
validation_result = vector("list", 4)
names(validation_result) = paste0(test_pct*100, "pct")
for (k in 1:length(test_pct)){
  fn <- paste0("Cross_validation_", test_pct[k]*100, "pct.txt")
  if (file.exists(fn)) file.remove(fn)
  validation_result[[k]] = crossValidation(response = Resp, finalmodel = mod_final, A_mat = out_mesh$A_mat, spde = out_mesh$spde, family = "binomial", raster_stack = covariate_stack, int.strategy = "eb", n_reps = n_reps, pct_out = test_pct[k])
}
validation_result


## Plot validation results
pdf("Cross_validation_result.pdf")

par(mfrow=c(2,1), mar=c(3.5,4,3,1), mgp=c(2,1,0))
## Pearson correlation 
DF = data.frame(Cor_10=validation_result$`10pct`[,1], Cor20=validation_result$`20pct`[,1], Cor30=validation_result$`30pct`[,1], Cor40=validation_result$`40pct`[,1])
boxplot(DF, names = c("10%", "20%", "30%", "40%"), ylab=expression(rho), xlab="Percentage of validation set", col=2:6, main="(A) Pearson correlation coefficient", cex.main=1)
## Cross-validated R-squared
DF2 = data.frame(Cor_10=validation_result$`10pct`[,2], Cor20=validation_result$`20pct`[,2], Cor30=validation_result$`30pct`[,2], Cor40=validation_result$`40pct`[,2])
boxplot(DF2, names = c("10%", "20%", "30%", "40%"), ylab=expression(R^2), xlab="Percentage of validation set", col=2:6, main="(B) Cross-validated R-squared", cex.main=1)

dev.off()


