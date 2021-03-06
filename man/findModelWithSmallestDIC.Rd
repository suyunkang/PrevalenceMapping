% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllTheFunctions.R
\name{findModelWithSmallestDIC}
\alias{findModelWithSmallestDIC}
\title{Find the best spatial model, i.e. the model with the smallest DIC}
\usage{
findModelWithSmallestDIC(response, raster_stack, A_mat, spde,
  family = "binomial", save_output = TRUE)
}
\arguments{
\item{response}{A data frame containing the response data including number of positive cases (n_positive), number of individuals examined (examined), and point locations (longitude and latitude).}

\item{raster_stack}{A collection of RasterLayer objects with the same spatial extent and resolution.}

\item{A_mat}{An observation/prediction weight matrix returned from makeMeshSPDE.}

\item{spde}{An inla.spde2 model object for a Matern model returned from makeMeshSPDE.}

\item{family}{A string indicating the likelihood family, Default: `binomial'.}

\item{save_output}{A logical argument indicating if output to be saved to a csv file, Default: TRUE.}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
The function performs forward and backward elimination of variables in order to find the spatial model with a set of covariates that result in the smallest DIC.
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
