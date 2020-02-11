\name{growth}
\alias{growth}
\title{
  Berkeley Growth Study data.
}
\description{
  A list containing the heights of 39 boys and 54 girls from age 1 to 18
  and the ages at which they were collected.
}
\format{
    This list contains the following components:
    \describe{
	\item{hgtm}{
	  a 31 by 39 numeric matrix giving the heights in centimeters of
	  39 boys at 31 ages.
	}
	\item{hgtf}{
	  a 31 by 54 numeric matrix giving the heights in centimeters of
	  54 girls at 31 ages.
	}
	\item{age}{
	  a numeric vector of length 31 giving the ages at which the
	  heights were measured.
	}
    }
}
\details{
These data are included here for convenience as they are also available
in the "fda" R-package.
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2006,ISBN 978-0-387-22751-1), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York.

  Ramsay, James O., and Silverman, Bernard W. (2002, ISBN 978-0-387-22465-7), \emph{Applied
    Functional Data Analysis}, Springer, New York, ch. 6.

  Tuddenham, R. D., and Snyder, M. M. (1954) "Physical growth of
  California boys and girls from birth to age 18", \emph{University of
  California Publications in Child Development}, 1, 183-364.

  Ramsay, James O., Wickham, Hadley, Graves, Spencer, and Hooker, Giles (2018)<https://CRAN.R-project.org/package=fda>,
    \emph{fda: Functional Data Analysis},R package version 2.4.8,

}
%\seealso{\code{\link{eval.penalty}}}
%\examples{
%}
% docclass is function
%\keyword{datasets}
