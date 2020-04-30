\name{swe.jo09}
\alias{swe.jo09}

\title{Statistical SWE modeling depending on month and climatic region in Switzerland}

\description{
  Snow Water Equivalent (SWE) is modeled statistically depending on snow depth, altitude, date and region in Switzerland.
}
\usage{
swe.jo09(data, alt, region.jo09)
}

\arguments{
  \item{data}{A data.frame of daily observations with two columns named \emph{date} and \emph{hs} referring to day and snow depth at that day. The date column must be a character string with the format \code{YYYY-MM-DD}. The hs column must be snow depth values \eqn{\ge 0} in m. }
  \item{alt}{station elevation in meters}
  \item{region.jo09}{Integer number of the Swiss region where the station belongs to, according to Fig. 1 in the original reference. Must be one of 1,2,3,4,5,6,7.}
}

\details{
\code{swe.jo09}{ This model parametrizes bulk snow density using snow depth, season (i.e. month), site altitude and site location. The location is implemented by a density offset according to the region in Switzerland, where the station belongs to. Non computable values are returned as NA.}
}

\value{
  A list with SWE values for each region in mm.
}

\references{
 	Jonas, T., Marty, C. and Magnusson, J. (2009) 'Estimating the snow water equivalent from snow depth measurements in the Swiss Alps', Journal of Hydrology, 378(1 - 2), pp. 161 - 167. doi: 10.1016/j.jhydrol.2009.09.021.
}


\concept{swe}
\concept{snow depth}
\concept{statistical}
\concept{modeling}
