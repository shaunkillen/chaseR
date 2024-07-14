\name{window_effect}
\alias{window_effect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
window_effect
}
\description{
A function to calculate rolling slopes for numerous rolling window sizes. Returns a plot of steepest estimated decline in water oxygen content
versus size of rolling window, colour-scaled to the mean r-squared value for each window size. Also displays the minimum rolling window size to achieve threshold mean r-squared values. Data for up to four animals can be analysed at once.
}
\usage{
window_effect(data, exclude_time = x, window_durations = seq(lower limit, upper limit, increment size), r_squared_threshold = x)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exclude_time}{
Time (in seconds) to be exlcuded from the beginning of the dataset.
}
 \item{window_durations}{
Range of times (in seconds) of rolling window sizes to be included in the analysis. User specifies the lower and upper limits of times to be
included, as well as the increments for all times to be tested in between these limits.
 \item{r_squared_threshold}{
The minimum r-sqaured value desired for selecting a rolling window duration.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For each animal, returns a plot of the minimum slope after subtraction of the blank (the steepest decline in water oxygen content) versus window duration. Points are scaled accroding to the
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

window_effect(dataset, exclude_time = 30, window_durations = seq(30, 600, 30), r_squared_threshold = 0.90)
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
