% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllTheFunctions.R
\name{fitFinalModel}
\alias{fitFinalModel}
\title{Fit the final spatial model using a subset of data}
\usage{
fitFinalModel(response, val_DAT, finalmodel, A_mat, spde, family = "binomial",
  raster_stack, int.strategy)
}
\arguments{
\item{response}{A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).}

\item{val_DAT}{Subset of data for cross-validation.}

\item{finalmodel}{An object returned from the function findModelWithSmallestDIC.}

\item{A_mat}{An observation/prediction weight matrix returned from makeMeshSPDE.}

\item{spde}{An inla.spde2 model object for a Matern model returned from makeMeshSPDE.}

\item{family}{A string indicating the likelihood family, Default: `binomial'.}

\item{raster_stack}{A collection of RasterLayer objects with the same spatial extent and resolution.}

\item{int.strategy}{Character. The integration strategy to use; one of `auto', `ccd', `grid', `eb' (empirical bayes), `user' or `user.std'.}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
Fit the final spatial model using a subset of data for the purpose of cross-validation.
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
