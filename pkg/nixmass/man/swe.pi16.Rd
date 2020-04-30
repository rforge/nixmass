\name{swe.pi16}
\alias{swe.pi16}

\title{Statistical SWE modeling depending on the day-of-year}

\description{
  This model parameterizes bulk snow density with day-of-the-year as the only input. It was calibrated for the region of South Tyrol, Italy, and is therefore called ST model in the original reference.
}
\usage{
swe.pi16(data, rho_0=200, K=1)
}

\arguments{
  \item{data}{A data.frame of daily observations with two columns named \emph{date} and \emph{hs} referring to day and snow depth at that day. The date column must be a character string with the format \code{YYYY-MM-DD}. The hs column must be snow depth values \eqn{\ge 0} in m. }
  \item{rho_0}{Intercept of the linear regression between observed snow depths and SWE values. \code{rho_0} is set to 200 as default, which is the value from the original reference. It can however be set to any value according to regression modeling with other datasets.}
  \item{K}{Slope of the linear regression between observed densities and the day-of-year as defined in the original reference. \code{K} is set to 1 as default, which is the value from the original reference. It can however be set to any value according to regression modeling with other datasets. }
}

\details{
\code{swe.pi16}{ This function uses only the day-of-year (DOY) as parameterization for bulk snow density and hence SWE. Here, the datums in the input data.frame are converted to DOY as defined in the original reference: negative values between 1.10. and 31.12. DOY=-92 at 1.10. In leap years 31.12. has DOY = 0, in non-leap years 31.12. has DOY = -1 with no day being 0. Non computable values are returned as NA.}
}

\value{
  A vector with daily SWE values in mm.
}

\references{
 	Pistocchi, A. (2016) 'Simple estimation of snow density in an Alpine region', Journal of Hydrology: Regional Studies. Elsevier B.V., 6(Supplement C), pp. 82 - 89. doi: 10.1016/j.ejrh.2016.03.004.
}


\concept{swe}
\concept{snow depth}
\concept{statistical}
\concept{modeling}
\concept{doy}
