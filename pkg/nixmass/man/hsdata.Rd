\name{hsdata}
\alias{hsdata}
\docType{data}

\title{Daily snow height data for a northern alpine station}

\description{
  Gapless daily snow height observations for a winter season from 1.8. - 31.7. from a station situated in the northern earstern alps at an altitude of 600 m. For anonymization the years are intentionally set to 1900 - 1901. 
  
This data series is free of gaps with a minimum of 0 and a maximum of 1.3 meters. It is intended to be used as is as input data for the package \code{\link{nixmass}} to calculate snow water equivalent and bulk snow density with the delta.snow method and several empirical regression models from the literature. 
}
\usage{
data("hsdata")
}

\format{
A \code{data.frame} named \code{data} with columns \code{date} and \code{hs}. The date column contains character strings of the format "YYYY-MM-DD" and is of class \code{character}. The hs column holds daily observed snow heights in meters and is of class \code{numeric}. 
}


\examples{
## Load example data 
data("hsdata")

## explore dataset
head(hsdata)
plot(hsdata$hs, type="o")

## compute snow water equivalents
o <- nixmass(hsdata, model="delta.snow",verbose=TRUE)
plot(o)
     
o1 <- nixmass(hsdata, alt=600, region.jo09 = 6, region.gu19 = "central",
             snowclass.st10 = "alpine", verbose = FALSE)
plot(o1)
summary(o1)

}

