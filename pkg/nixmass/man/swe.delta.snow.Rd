\name{swe.delta.snow}
\alias{swe.delta.snow}


\title{SWE modeling from daily snow height differences}

\description{
  Model daily values of Snow Water Equivalent (SWE) solely from daily differences of snow depth.
}
\usage{
swe.delta.snow(data, rho.max=401.2588, rho.null=81.19417, c.ov=0.0005104722, 
               k.ov=0.37856737,k=0.02993175, tau=0.02362476, eta.null=8523356,
               timestep=24, verbose=FALSE)
}


\arguments{
  \item{data}{A data.frame with at least two columns named \code{date} and \code{hs}. They should contain date and corresponding daily observations of snow depth \eqn{hs \ge 0} measured at one site. The unit must be meters (m). No gaps or NA are allowed.}
  \item{rho.max}{Maximum density of an individual snow layer produced by the deltasnow model [kg/m3], \eqn{rho.max > 0}}
  \item{rho.null}{Fresh snow density for a newly created layer [kg/m3], \eqn{rho.null > 0}}
  \item{c.ov}{Overburden factor due to fresh snow [-], \eqn{c.ov > 0}}
  \item{k.ov}{Defines the impact of the individual layer density on the compaction due to overburden [-], \eqn{k.ov \in [0,1]}.}
  \item{k}{Exponent of the exponential-law compaction [m3/kg], \eqn{k > 0}.}
  \item{tau}{Uncertainty bound [m], \eqn{tau > 0}.}
  \item{eta.null}{Effective compactive viscosity of snow for "zero-density" [Pa s].}
  \item{timestep}{Timestep between snow depth observations in hours. Default is 24 hours, i.e. daily snow depth observations.}
  \item{verbose}{Should additional information be given during runtime? Can be \code{TRUE} or \code{FALSE}.}
}

\details{
\code{swe.delta.snow}{ computes SWE solely from daily changes of snow depth at an observation site. \cr
Compression of a snow layer without additional load on top is computed on the basis of Sturm and Holmgren (1998), who regard snow as a viscous fluid: \cr
  \deqn{\rho_i(t_{i+1}) = \rho_i(t_i)*(1+(SWE*g)/\eta_0 * exp^{-k_2*\rho_i(t_i)})}
  with \eqn{\rho_i(t_{i+1}) and \rho_i(t_i)} being tomorrow's and today's respective density of layer i, the gravitational acceleration \eqn{g = 9.8ms^{-2}}, viscosity \eqn{\eta_0} [Pa] and factor \eqn{k2 [m^3kg^{-1}}], determining the importance of today's for tomorrow's density.

}
}

\value{
  A vector with daily SWE values in mm.
}

\author{Harald Schellander, Michael Winkler}

\references{
 	Gruber, S. (2014) "Modelling snow water equivalent based on daily snow depths", Masterthesis, Institute for Atmospheric and Cryospheric Sciences, University of Innsbruck.
 	\cr\cr
 	Martinec, J., Rango, A. (1991) "Indirect evaluation of snow reserves in mountain basins". Snow, Hydrology and Forests in High Alpine Areas. pp. 111-120.
 	\cr\cr
 	Sturm, M., Holmgren, J. (1998) "Differences in compaction behavior of three climate classes of snow". Annals of Glaciology 26, 125-130.
 	\cr\cr
 	Winkler, M., Schellander, H., and Gruber, S. (submitted) "Snow Water Equivalents exclusively from Snow Heights and their temporal Changes: The delta.snow.model". Hydrology and Earth System Sciences Discussions.
}


\concept{swe}
\concept{snow depth}
\concept{modeling}
