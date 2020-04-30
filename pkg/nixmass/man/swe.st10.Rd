\name{swe.st10}
\alias{swe.st10}

\title{Statistical SWE modeling depending on day of year and a climatic region}

\description{ The \emph{Sturm} model parametrizes bulk snow density with day of the year and a \emph{snowclass.st10}. It was trained on historical snow height - density - SWE data from the United States, Canada, and Switzerland.
}
\usage{
swe.st10(data,snowclass.st10)
}

\arguments{
  \item{data}{A data.frame of daily observations with two columns named \emph{date} and \emph{hs} referring to day and snow height at that day. The date column must be a character string with the format \code{YYYY-MM-DD}. The hs column must be snow height values \eqn{\ge 0} in m. }
  \item{snowclass.st10}{Must be one of the following character strings: "alpine","maritime","prairie","tundra","taiga".  }
}

\details{
\code{swe.st10}{ This model converts snow height to SWE using snow height, day of year and station location (from which a climate class of snow can be inferred. The day of year (DOY) is the day-number of in the season 1.10. - 30.6. The 1.10. refers to DOY = -92. The 1.2. would be DOY = 32, while 15.11. would be DOY = -47. The \emph{snowclass.st10} must be one out of the character strings "alpine","maritime","prairie","tundra" and "taiga". For the Alps probably "alpine" would be the most appropriate climate classification.
Non computable values are returned as NA.}
}

\value{
  A vector with daily SWE values in mm.
}

\references{
 	Sturm, M. et al. (2010) 'Estimating Snow Water Equivalent Using Snow Depth Data and Climate Classes', Journal of Hydrometeorology, 11(6), pp. 1380 - 1394. doi: 10.1175/2010JHM1202.1.
}


\concept{swe}
\concept{snow depth}
\concept{statistical}
\concept{modeling}
