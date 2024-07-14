\name{slopes}
\alias{slopes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
slopes
}
\description{
This function calculates all slopes of water oxygen content versus time, using a rolling window.
Also returned are the start times for each slope and r-sqaured values. Blank slope is also calculated
and subtracted from all rolling window slopes. Data for up to four animals can be processed at once.
}
\usage{
slopes(data, exclude_time = x, window_duration = x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exclude_time}{
time (in seconds) to be exlcuded from the beginning of the dataset
}
\item{window_duration}{
time duration (in seconds) of the rolling window for calculating slopes
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns all slopes, and their corresponding start times and r-squared values. Also given is the minimum slope (the steepest decline in O2 over time), the start time for th specific slope, and it's r-squared.

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Shaun Killen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(123)
time <- seq(0, 600, by = 2)
MO2_fish1 <- runif(length(time), 7.5, 11.0)
MO2_fish2 <- runif(length(time), 7.5, 11.0)
MO2_fish3 <- runif(length(time), 7.5, 11.0)
MO2_fish4 <- runif(length(time), 7,5, 11.0)
MO2_blank <- runif(length(time), 10.8, 11.0)

slopes(dataset, exclude_time = 30, window_duration = 200)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
