\name{nixmass}
\alias{nixmass}

\title{SWE modeling with the delta.snow process based model and several emprical regression models.}

\description{
  Snow Water Equivalent (SWE) is modeled either exclusiveley from daily snow depth changes or statistically depending on snow depth, altitude, date and climate class.
}
\usage{
nixmass(data, model = c("delta.snow","jo09","pi16","st10","gu19"), 
       alt, region.jo09, region.gu19, snowclass.st10, verbose = FALSE)
}

\arguments{
 \item{data}{A data.frame of daily observations with two columns named \emph{date} and \emph{hs} referring to day and snow depth at that day. Values in the date column must be of class character with format \code{YYYY-MM-DD}. Values in the hs column must be snow depth values of class numeric \eqn{\ge 0} in m. No gaps or NA are allowed.}
  \item{model}{Defines model for SWE computation. Can be one, several or all of "delta.snow","jo09","pi16","st10","gu19". If no model is given, all models are computed.}
    \item{alt}{Must be given if one of model is \code{"jo09"}. Station elevation in meters}
  \item{region.jo09}{Must be given if one of model is \code{"jo09"}. This must be an integer number between 1 and 7 of the Swiss region where the station belongs to, according to Fig. 1 in the original reference. }
  \item{region.gu19}{If model contains \code{"gu19"} this must be one of "italy", "southwest", "central" or "southeast" as described in the original reference.}
  \item{snowclass.st10}{Must be given if one of model is \code{"st10"}. Must be one of the following character strings: "alpine","maritime","prairie","tundra","taiga" as outlined in the original reference. }
  \item{verbose}{Logical. Should additional information be given during runtime?}
}

\details{
\code{nixmass}{ This function is a wrapper for the computation of SWE with different models. The process based model \code{\link[=swe.delta.snow]{delta.snow}} can be chosen, as well as different empirical regression models of \code{\link[=swe.jo09]{Jonas},\link[=swe.pi16]{Pistocchi}, \link[=swe.st10]{Sturm}} and \link[=swe.gu19]{Guyennon}}.
For the "delta.snow" model and the ones of "Pistocchi" and "Guyennon", the needed parameters and coefficients from the original references are set as default. They can however be changed according to results from other datasets. For the other models of "Jonas" and "Sturm" regression coefficients are fixed. 

In its current implementation, matrices holding the snowpack and all its layers are predefined. This makes the iteration over the whole provided period of snow height values slow, if many years have to be modeled at once. It is therefore strongly recommended to split computation of SWE values into single years, separated by zero values for snow height. In a future implementation, this will be speeded up. By computing only two days at once, pre-allocation of matrices can be restricted, but may not be completely avoided.
}

\value{
  A list of class \code{"nixmass"} with components:
  \item{swe}{ Contains a list of numerical vectors. Each entry refers to SWE values computed with the selected model(s). }
  \item{date}{Character vector of date strings in the format \code{YYYY-MM-DD}.}
  \item{hs}{Vector of snow height values used to compute SWE. }
}

\author{Harald Schellander, Michael Winkler}

\references{
     Guyennon, N., Valt, M., Salerno, F., Petrangeli, A., Romano, E. (2019) 'Estimating the snow water equivalent from snow depth measurements in the Italian Alps', Cold Regions Science and Technology. Elsevier, 167 (August), p. 102859. doi: 10.1016/j.coldregions.2019.102859.
\cr\cr     
 	Jonas, T., Marty, C. and Magnusson, J. (2009) "Estimating the snow water equivalent from snow depth measurements in the Swiss Alps"", Journal of Hydrology, 378(1 - 2), pp. 161 - 167. doi: 10.1016/j.jhydrol.2009.09.021.
\cr\cr 	
 	Pistocchi, A. (2016) "Simple estimation of snow density in an Alpine region", Journal of Hydrology: Regional Studies. Elsevier B.V., 6(Supplement C), pp. 82 - 89. doi: 10.1016/j.ejrh.2016.03.004.
\cr\cr 	
 	Sturm, M. et al. (2010) "Estimating Snow Water Equivalent Using Snow Depth Data and Climate Classes", Journal of Hydrometeorology, 11(6), pp. 1380 - 1394. doi: 10.1175/2010JHM1202.1.
\cr\cr
Winkler, M., Schellander, H., and Gruber, S. (submitted) "Snow Water Equivalents exclusively from Snow Heights and their temporal Changes: The delta.snow.model". Hydrology and Earth System Sciences Discussions.
}


\examples{
## Load example data with realistic snow height values 
## from a station at 600 meters in the northern Alps
## Note that the winter season is set to an arbitrary date 
## to mask its origin
data("hsdata")
o <- nixmass(hsdata, model="delta.snow",verbose=TRUE)
plot(o)
     
o1 <- nixmass(hsdata, alt=600, region.jo09=6, region.gu19 = "central",
             snowclass.st10 = "alpine", verbose = FALSE)
plot(o1)
summary(o1)
}


\concept{swe}
\concept{snow depth}
\concept{statistical}
\concept{modeling}

