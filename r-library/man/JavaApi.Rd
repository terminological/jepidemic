% Generated by r6-generator-maven-plugin: do not edit by hand

\name{JavaApi}

\alias{JavaApi}

\title{Java Epidemic}

\section{Usage}{
	\if{html}{\out{<div class="r">}}
	\preformatted{
 J = jepidemic::JavaApi$get(logLevel)
 	}
  \if{html}{\out{</div>}}
}

\arguments{
	\if{html}{\out{<div class="arguments">}}
	\item{logLevel}{optional - the slf4j log level as a string - one of OFF (most specific, no logging), 
	FATAL (most specific, little data), ERROR, WARN, INFO, DEBUG, 
	TRACE (least specific, a lot of data), ALL (least specific, all data)}
	\if{html}{\out{</div>}}
}

\description{
Some work in progress R functions related to the covid 19 epidemic

Version: 0.03

Generated: 2022-05-24T22:12:30.918083
}

\author{\email{rc538@exeter.ac.uk}}

\keyword{java api}

\section{Static methods and constructors}{
	\itemize{
		\item \code{JavaApi$get()}
		\item \href{#method-GrowthRateEstimator-new}{\code{J$GrowthRateEstimator$new(minWindow, maxWindow)}}
		\item \href{#method-GrowthRateEstimator-strictGrowthRateAndRtEstimator}{\code{J$GrowthRateEstimator$strictGrowthRateAndRtEstimator(minWindow, maxWindow, infectivityProfiles)}}
		\item \href{#method-GrowthRateEstimator-defaultGrowthRateAndRtEstimator}{\code{J$GrowthRateEstimator$defaultGrowthRateAndRtEstimator(minWindow, maxWindow, initialIncidence)}}
		\item \href{#method-GrowthRateEstimator-defaultGrowthRateEstimator}{\code{J$GrowthRateEstimator$defaultGrowthRateEstimator(minWindow, maxWindow)}}
		\item \href{#method-GrowthRateEstimator-basicGrowthRateEstimator}{\code{J$GrowthRateEstimator$basicGrowthRateEstimator()}}
		\item \href{#method-Serialiser-new}{\code{J$Serialiser$new()}}
		\item \href{#method-Serialiser-serialiseDataframe}{\code{J$Serialiser$serialiseDataframe(dataframe, filename)}}
		\item \href{#method-Serialiser-deserialiseDataframe}{\code{J$Serialiser$deserialiseDataframe(filename)}}
		\item \href{#method-Serialiser-serialiseList}{\code{J$Serialiser$serialiseList(dataframe, filename)}}
		\item \href{#method-Serialiser-deserialiseList}{\code{J$Serialiser$deserialiseList(filename)}}
		\item \href{#method-Serialiser-serialiseNamedList}{\code{J$Serialiser$serialiseNamedList(dataframe, filename)}}
		\item \href{#method-Serialiser-deserialiseNamedList}{\code{J$Serialiser$deserialiseNamedList(filename)}}
		\item \href{#method-InfectivityProfile-new}{\code{J$InfectivityProfile$new(discretePdf, id)}}
		\item \href{#method-CoriEstimator-new}{\code{J$CoriEstimator$new(r0Mean, r0SD, maxWindow)}}
		\item \href{#method-CoriEstimator-defaultEpiEstim}{\code{J$CoriEstimator$defaultEpiEstim(infectivityProfile, meanR0, sdR0, window)}}
		\item \href{#method-CoriEstimator-defaultUncertainEpiEstim}{\code{J$CoriEstimator$defaultUncertainEpiEstim(infectivityProfiles, meanR0, sdR0, window, n2)}}
	}
	

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-GrowthRateEstimator-new"></a>}}
	
	\subsection{Method \code{GrowthRateEstimator$new()}}{
		Produces a GrowthRateEstimator that is configured to behave in a sensible fashion, with sane defaults. 
In that it estimates an initial poisson rate prior from the first entry and returns a fixed window size and aggregates uncertain SI 
distributions with a random resampling strategy<br><br>
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$GrowthRateEstimator$new(minWindow, maxWindow)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{minWindow - the smallest amount of data considered in estimating the growth rate, 4 is a sensible default unless you want bang up to date estimate.}{ - (java expects a int)}
					\item{maxWindow - the largest amount of data considered in estimating the growth rate, 14 is a sensible amount.}{ - (java expects a int)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 GrowthRateEstimator object: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-GrowthRateEstimator-strictGrowthRateAndRtEstimator"></a>}}
	
	\subsection{Method \code{GrowthRateEstimator$strictGrowthRateAndRtEstimator()}}{
		Produces a GrowthRateEstimator that will estimate growth rate and $R_t$, with sane defaults. The $R_t$ estimation is 
done using bootstrapping, and combining the results using a mixture distribution. This will be relatively slow.
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$GrowthRateEstimator$strictGrowthRateAndRtEstimator(minWindow, maxWindow, infectivityProfiles)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.}{ - (java expects a int)}
					\item{maxWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount.}{ - (java expects a int)}
					\item{infectivityProfiles - a numeric matrix with each column containing a discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)}{ - (java expects a RNumericArray)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 GrowthRateEstimator object: a strict growth rate and \(R_t\) estimator that accounts for uncertainty in the infectivity profile.
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-GrowthRateEstimator-defaultGrowthRateAndRtEstimator"></a>}}
	
	\subsection{Method \code{GrowthRateEstimator$defaultGrowthRateAndRtEstimator()}}{
		Produces a GrowthRateEstimator that will estimate growth rate and $R_t$, with sane defaults. The $R_t$ estimation is 
done using a single infectivity profile and does not account for infectivity profile uncertainty. This is comparatively fast.<br><br>
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$GrowthRateEstimator$defaultGrowthRateAndRtEstimator(minWindow, maxWindow, initialIncidence)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.}{ - (java expects a int)}
					\item{maxWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount.}{ - (java expects a int)}
					\item{initialIncidence - a guess at the initial incidence of the time series}{ - (java expects a double)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 GrowthRateEstimator object: a growth rate and \(R_t\) estimator that does not account for uncertainty in the infectivity profile.
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-GrowthRateEstimator-defaultGrowthRateEstimator"></a>}}
	
	\subsection{Method \code{GrowthRateEstimator$defaultGrowthRateEstimator()}}{
		Produces a GrowthRateEstimator that will estimate growth rate only, with sane defaults, unless an infectivity profile is lated added<br><br>
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$GrowthRateEstimator$defaultGrowthRateEstimator(minWindow, maxWindow)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.}{ - (java expects a int)}
					\item{maxWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount.}{ - (java expects a int)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 GrowthRateEstimator object: a growth rate and \(R_t\) estimator that does not account for uncertainty in the infectivity profile.
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-GrowthRateEstimator-basicGrowthRateEstimator"></a>}}
	
	\subsection{Method \code{GrowthRateEstimator$basicGrowthRateEstimator()}}{
		Produces a GrowthRateEstimator that will estimate growth rate only, with sane defaults, and estimating using windows between 4 and 14 days<br><br>
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$GrowthRateEstimator$basicGrowthRateEstimator()
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{none}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 GrowthRateEstimator object: a growth rate and \(R_t\) estimator that does not account for uncertainty in the infectivity profile.
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-new"></a>}}
	
	\subsection{Method \code{Serialiser$new()}}{
		the default no-args constructor
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$new()
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{none}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 Serialiser object: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-serialiseDataframe"></a>}}
	
	\subsection{Method \code{Serialiser$serialiseDataframe()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$serialiseDataframe(dataframe, filename)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{dataframe}{ - (java expects a RDataframe)}
					\item{filename}{ - (java expects a String)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			void: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-deserialiseDataframe"></a>}}
	
	\subsection{Method \code{Serialiser$deserialiseDataframe()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$deserialiseDataframe(filename)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{filename}{ - (java expects a String)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			RDataframe: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-serialiseList"></a>}}
	
	\subsection{Method \code{Serialiser$serialiseList()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$serialiseList(dataframe, filename)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{dataframe}{ - (java expects a RList)}
					\item{filename}{ - (java expects a String)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			void: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-deserialiseList"></a>}}
	
	\subsection{Method \code{Serialiser$deserialiseList()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$deserialiseList(filename)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{filename}{ - (java expects a String)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			RList: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-serialiseNamedList"></a>}}
	
	\subsection{Method \code{Serialiser$serialiseNamedList()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$serialiseNamedList(dataframe, filename)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{dataframe}{ - (java expects a RNamedList)}
					\item{filename}{ - (java expects a String)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			void: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-Serialiser-deserialiseNamedList"></a>}}
	
	\subsection{Method \code{Serialiser$deserialiseNamedList()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$Serialiser$deserialiseNamedList(filename)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{filename}{ - (java expects a String)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			RNamedList: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-InfectivityProfile-new"></a>}}
	
	\subsection{Method \code{InfectivityProfile$new()}}{
		no description
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$InfectivityProfile$new(discretePdf, id)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{discretePdf - A discrete pdf of probabilities over time of secondary infection given priamry infection where probability at time 0 is 0}{ - (java expects a RNumericVector)}
					\item{id - a numeric index or id for this profile}{ - (java expects a int)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 InfectivityProfile object: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-CoriEstimator-new"></a>}}
	
	\subsection{Method \code{CoriEstimator$new()}}{
		The Estimator is configured with the baseline \(R_0\) assumptions which are used as prior values for \(R_t\) when no other 
values are available, or the {@link CoriEstimator#withDefaultPrior()} option is used.
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$CoriEstimator$new(r0Mean, r0SD, maxWindow)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{r0Mean - The mean of the gamma distribution of the baseline prior assumption about R0 - the theorerical reproduction number in an infinite pool of evenly mixing susceptible individuals}{ - (java expects a double)}
					\item{r0SD - The SD of the prior gamma distribution}{ - (java expects a double)}
					\item{maxWindow - The maximum sliding window over which to estimate \(R_t\).}{ - (java expects a int)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 CoriEstimator object: 
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-CoriEstimator-defaultEpiEstim"></a>}}
	
	\subsection{Method \code{CoriEstimator$defaultEpiEstim()}}{
		Produces a CoriEstimator that is configured to behave like the default settings of EpiEstim. 
In that it uses a fixed prior for Rt and returns a fixed window size<br><br>
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$CoriEstimator$defaultEpiEstim(infectivityProfile, meanR0, sdR0, window)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{infectivityProfile - a vector with the discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)}{ - (java expects a RNumericVector)}
					\item{meanR0 - the mean prior of R0}{ - (java expects a double)}
					\item{sdR0 - the sd prior of R0}{ - (java expects a double)}
					\item{window - the window size}{ - (java expects a int)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 CoriEstimator object: a default cori estimator
		}
	}

	\if{html}{\out{<hr>}}
	\if{html}{\out{<a id="method-CoriEstimator-defaultUncertainEpiEstim"></a>}}
	
	\subsection{Method \code{CoriEstimator$defaultUncertainEpiEstim()}}{
		Produces a CoriEstimator that is configured to behave like the default settings of EpiEstim. 
In that it uses a fixed prior for Rt and returns a fixed window size and aggregates uncertain SI 
distributions with a random resampling strategy<br><br>
	
		\subsection{Usage}{
			\if{html}{\out{<div class="r">}}
			\preformatted{
J = jepidemic::JavaApi$get()
J$CoriEstimator$defaultUncertainEpiEstim(infectivityProfiles, meanR0, sdR0, window, n2)
	  		}
			\if{html}{\out{</div>}}
		}
	
		\subsection{Arguments}{
			\if{html}{\out{<div class="arguments">}}
			\describe{
				\itemize{
					\item{infectivityProfiles - a numeric matrix with the discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)}{ - (java expects a RNumericArray)}
					\item{meanR0 - the mean prior of R0}{ - (java expects a double)}
					\item{sdR0 - the sd prior of R0}{ - (java expects a double)}
					\item{window - the window size}{ - (java expects a int)}
					\item{n2 - the number of samples to use for estimating confidence limits on bootstrapped infectivity profiles}{ - (java expects a int)}
				}
			}
			\if{html}{\out{</div>}}
		}

		\subsection{Returns}{
			R6 CoriEstimator object: a default cori estimator
		}
	}

}