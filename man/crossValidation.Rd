% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllTheFunctions.R
\name{crossValidation}
\alias{crossValidation}
\title{Performs cross-validation}
\usage{
crossValidation(response, finalmodel, A_mat, spde, family = "binomial",
  raster_stack, int.strategy, n_reps = 100, pct_out)
}
\arguments{
\item{response}{A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).}

\item{finalmodel}{An object returned from the function findModelWithSmallestDIC.}

\item{A_mat}{An observation/prediction weight matrix returned from makeMeshSPDE.}

\item{spde}{An inla.spde2 model object for a Matern model returned from makeMeshSPDE.}

\item{family}{A string indicating the likelihood family, Default: `binomial'.}

\item{raster_stack}{A collection of RasterLayer objects with the same spatial extent and resolution.}

\item{int.strategy}{Character. The integration strategy to use; one of `auto', `ccd', `grid', `eb' (empirical bayes), `user' or `user.std'.}

\item{n_reps}{Number of replicates of subsets of data for cross-validation, Default: `100'.}

\item{pct_out}{Percentage of data to be used as a test set.}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
Performs cross-validation for a subset of data to assess predictive performance of the final model.
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
