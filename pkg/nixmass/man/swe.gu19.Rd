\name{swe.gu19}
\alias{swe.gu19}

\title{Statistical SWE modeling based on a quadratic dependance on the day-of-year}

\description{
  This model parameterizes bulk snow density with day-of-the-year as the only input similar to \code{\link{swe.pi16}} but adds a quadratic dependance. It was calibrated for the regions of the whole Italian alps, and the subregions South-West, Central and South-East. By setting the cofficients of the empirical regression it can however be used with results from other datasets.
}
\usage{
swe.gu19(data, region.gu19, n0=NA ,n1=NA, n2=NA)
}

\arguments{
  \item{data}{A data.frame of daily observations with two columns named \emph{date} and \emph{hs} referring to day and snow depth at that day. The date column must be a character string with the format \code{YYYY-MM-DD}. The hs column must be snow depth values \eqn{\ge 0} in m. }
  \item{region.gu19}{ Must be one of the italian subalpine regions \emph{italy}, \emph{southwest}, \emph{central} or \emph{southeast}, defined in the original reference (see details), or \emph{myregion}, in which case the coefficients n0, n1 and n2 have to be set. }
  \item{n0,n1,n2}{Values \code{no}, \code{n1} and \code{n2} refer to Intercept, slope and quadratic dependence of an empirical regression between densities and the day-of-year (see details) 
  
  
  }
}

\details{
\code{swe.gu19}{ Similar to the model of Pistocchi (2016), this function uses only the day-of-year (DOY) as parameterization for bulk snow density and hence SWE. In contrast to the latter, here, a quadratic term for DOY was added, to reflect non-linearity in the snow bulk density variability. The datums in the input data.frame are converted to DOY as days spent since November 1st. Regression coefficients depend on regions defined in Guyennon et al. (2019), which are \emph{italy} for the Italian Alps, \emph{southwest} for the South-western Italian Alps, \emph{central} for the Central Italian Alpes or \emph{southeast} for the South-western Italian Alps. 

If \code{region.gu19} is set to \emph{myregion}, the coefficients \code{no}, \code{n1} and \code{n2} must be set to values, obtained from a regression between densities and day-of-year from another dataset. It has to have the form density ~ DOY + DOY^2, where DOY is the day-of-year as defined in the original reference. 
Non computable values are returned as NA.}
}

\value{
  A vector with daily SWE values in mm.
}

\references{
     Guyennon, N., Valt, M., Salerno, F., Petrangeli, A., Romano, E. (2019) 'Estimating the snow water equivalent from snow depth measurements in the Italian Alps', Cold Regions Science and Technology. Elsevier, 167 (August), p. 102859. doi: 10.1016/j.coldregions.2019.102859.

 	Pistocchi, A. (2016) 'Simple estimation of snow density in an Alpine region', Journal of Hydrology: Regional Studies. Elsevier B.V., 6 (Supplement C), pp. 82 - 89. doi: 10.1016/j.ejrh.2016.03.004.
}


\concept{swe}
\concept{snow depth}
\concept{statistical}
\concept{modeling}
\concept{doy}
