\name{calculate_variables}
\alias{calculate_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate_variables
}
\description{
Uses slopes from the "slopes" function to calculate MMR, the corresponding r-squared for the slope used to calculate MMR, the time at which MMR occurred, the exponential rate of decline in oxygen uptake throughout the data set, the aymptote for oxygen uptake during recovery, and optionally (if values for SMR are supplied) the excess post-expercise oxygen consumption (EPOC; calculated as the area under the curve of MO2 vs time). Returns plots of MO2 vs time, and EPOC for each animal. Up to four animals for be analysed at once.
}
\usage{
calculate_variables(data, tubing_diameters, tubing_lengths, V_chambers, fish_masses, exclude_time = 10, window_duration = 120, smr_values)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tubing_diameters}{
Diameters of tubing used in the mixing circuit of the respirometry setup (in cm).
}
  \item{tubing_lengths}{
Lengths of tubing used in the mixing circuit of the respirometry setup (in cm).
}
  \item{V_chambers}{
Volumes of respirometry chambers (in mL).
  \item{fish_masses}{
Masses of fish (in g).
  \item{exclude_time}{
time (in seconds) to be exlcuded from the beginning of the dataset
}
  \item{window_duration}{
time duration (in seconds) of the rolling window for calculating slopes
}
  \item{smr_values}{
Standard metabolic rate values for each fish (in mg O2/h). This is optional but required if calculating EPOC.
}
}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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

 # Example input values for each fish column
  tubing_diameters <- c(0.03, 0.03, 0.03, 0.03)  # in cm
  tubing_lengths <- c(50, 50, 50, 50)  # in cm
  V_chambers <- c(50, 50, 50, 50)  # in mL
  fish_masses <- c(1.429, 0.924, 0.879, 1.324)  # in grams
  smr_values <- c(0.187, 0.095, 0.124, 0.085)  # Example SMR values for each MO2_fish column in mg O2/h

calculate_variables(data, tubing_diameters, tubing_lengths, V_chambers, fish_masses, exclude_time = 10, window_duration = 120, smr_values)
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
