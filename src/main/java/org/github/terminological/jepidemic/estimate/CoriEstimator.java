package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.github.terminological.jepidemic.IncompleteTimeseriesException;
import org.github.terminological.jepidemic.InfectivityProfile;
import org.github.terminological.jepidemic.ParameterOutOfRangeException;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.distributions.Summary;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.co.terminological.rjava.RClass;
import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.RFunctions;
import uk.co.terminological.rjava.RMethod;
import uk.co.terminological.rjava.UnconvertableTypeException;
import uk.co.terminological.rjava.UnexpectedNaValueException;
import uk.co.terminological.rjava.ZeroDimensionalArrayException;
import uk.co.terminological.rjava.types.RCharacter;
import uk.co.terminological.rjava.types.RDataframe;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;
import uk.co.terminological.rjava.types.RNumericArray;
import uk.co.terminological.rjava.types.RNumericVector;

/**
 * An estimator of the reproduction number (\(R_t\)) using a renewal equation method. 
 * The renewal equation method depends on a time series of infections, and on the infectivity profile - 
 * a measure of the probability that a secondary infection occurred on a specific day after the primary case, 
 * given a secondary infection occurred. A Bayesian framework is then used to update a prior probabilistic 
 * estimate of \(R_t\) on any given day with both information gained from the time series of infections in the 
 * epidemic to date and the infectivity profile to produce a posterior estimate of \(R_t\). This is based on the 
 * article:<br>
 * 
 * <br>A. Cori, N. M. Ferguson, C. Fraser, and S. Cauchemez, ‘A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics’, Am J Epidemiol, vol. 178, no. 9, pp. 1505–1512, Nov. 2013, doi: 10.1093/aje/kwt133.<br>
 * 
 * <br> and the implementation of it in the R Library EpiEstim. It assumes that infections (or other type of observation) are a 
 * Poisson distributed process and that the reproduction number is a gamma distributed quantity. The estimates of \(R_t\) may be 
 * aggregated over a specific "window" to provide a more stable estimate.  
 */
@RClass
public class CoriEstimator {

	static Logger log = LoggerFactory.getLogger(CoriEstimator.class);
	
	List<InfectivityProfile> infProf = new ArrayList<>();
	GammaParameters rZero;
	private int maxTau = 14;
	private boolean atStart = false;
	private boolean summarise = false;
	private boolean epiEstimMode = false;
	
	Function<CoriEstimationResultEntry,List<DatedRtEstimate>> priorSelectionStrategy;
	Function<CoriEstimationResultEntry,DatedRtEstimate> posteriorSelectionStrategy; 
	Function<CoriEstimationSummaryEntry,Summary> combiningStrategy;

	
	
	/**
	 * The Estimator is configured with the baseline \(R_0\) assumptions which are used as prior values for \(R_t\) when no other 
	 * values are available, or the {@link CoriEstimator#withDefaultPrior()} option is used.
	 * @param r0Mean - The mean of the gamma distribution of the baseline prior assumption about R0 - the theorerical reproduction number in an infinite pool of evenly mixing susceptible individuals 
	 * @param r0SD - The SD of the prior gamma distribution 
	 * @param maxWindow - The maximum sliding window over which to estimate \(R_t\).
	 */
	@RMethod 
	public CoriEstimator(double r0Mean, double r0SD, int maxWindow) {
		if (r0Mean < 0) throw new ParameterOutOfRangeException("The mean of the prior Rt must be positive");
		if (r0SD <= 0) throw new ParameterOutOfRangeException("The SD of the prior Rt must be greater than zero");
		if (maxWindow <= 1) throw new ParameterOutOfRangeException("The maxWindow parameter must be 2 or more. Probably it should be >=7 unless working with very large number of timeseries");
		this.rZero = new GammaMoments(r0Mean,r0SD).convert();
		this.maxTau = maxWindow;
		this.withDefaultPrior();
		this.selectSpecificWindow(maxWindow);
		this.collectFirst();
	}
	
	/**
	 * @param on - a boolean flag to switch on (TRUE) or off (FALSE) eipestim naming convention  
	 * @return the estimator itself (a fluent method)
	 */
	@RMethod
	public CoriEstimator legacySupport(boolean on) {
		this.epiEstimMode = on;
		return this;
	}
	
	

	/** 
	 * Produces a CoriEstimator that is configured to behave like the default settings of EpiEstim. 
	 * In that it uses a fixed prior for Rt and returns a fixed window size<br><br>
	 * 
	 * @param infectivityProfile - a vector with the discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)
	 * @param meanR0 - the mean prior of R0
	 * @param sdR0 - the sd prior of R0
	 * @param window - the window size
	 * @return a default cori estimator
	 */
	@RMethod 
	public static CoriEstimator defaultEpiEstim(RNumericVector infectivityProfile, double meanR0, double sdR0, int window) {
		CoriEstimator out = new CoriEstimator(meanR0,sdR0,window);
		
		return out.withDefaultPrior()
				.selectSpecificWindow(window)
				.atStartOfTimeseries()
				.withInfectivityProfile(infectivityProfile, true)
				.collectFirst()
				.legacySupport(true);
		
	}
	
	/** 
	 * Produces a CoriEstimator that is configured to behave like the default settings of EpiEstim. 
	 * In that it uses a fixed prior for Rt and returns a fixed window size and aggregates uncertain SI 
	 * distributions with a random resampling strategy<br><br>
	 * 
	 * 
	 * @param infectivityProfiles - a numeric matrix with the discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)
	 * @param meanR0 - the mean prior of R0
	 * @param sdR0 - the sd prior of R0
	 * @param window - the window size
	 * @param n2 - the number of samples to use for estimating confidence limits on bootstrapped infectivity profiles
	 * @return a default cori estimator
	 */
	@RMethod 
	public static CoriEstimator defaultUncertainEpiEstim(RNumericArray infectivityProfiles, double meanR0, double sdR0, int window, int n2) {
		CoriEstimator out = new CoriEstimator(meanR0,sdR0,window);
		out.withDefaultPrior()
			.selectSpecificWindow(window)
			.atStartOfTimeseries()
			.collectResampledQuantiles(n2)
			.legacySupport(true)
			.withInfectivityProfileMatrix(infectivityProfiles);
		return out;
		
	}
	
	
	
	/**
	 * The select methods define how the posterior estimate of \(R_t\) is selected. In this case a set value for the window is used. This must be less or equal to the 
	 * value of maxWindow supplied in the CoriEstimator::new constructor. The window defines how many data points prior to the time of the estimate are used to synthesise the posterior
	 * estimate of \(R_t\), and so tends to reduce the uncertainty around the value of \(R_t\), but longer values of window introduces temporal uncertainty into
	 * the estimates. This is as implemented in EpiEstim (with a default value of 7 for the window).
	 * @param window - a number of days
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator selectSpecificWindow(int window) {
		if (window-1 > this.maxTau) throw new ParameterOutOfRangeException("The window must be smaller than maxWindow (from the constructor)");
		this.posteriorSelectionStrategy = (res) -> res.forWindow(window);
		return this;
	}
	
	public CoriEstimator selectAll() {
		this.posteriorSelectionStrategy = null;
		return this;
	}
	
	/**
	 * The adaptive windowing strategy selects a posterior estimate of \(R_t\) by selecting the shortest window length (&gt; minWindow) which spans at least incidenceSum infection 
	 * (or observation) counts. This ensures that where numbers are low the estimation window is lengthened (up ot the value of maxWindow provided to the CoriEstimator::new constructor) to 
	 * include more data at the cost of temporal precision. This prevents stochastic noise from dominating the estimates of \(R_t\). 
	 * @param incidenceSum - the number of cases that must be observed in a window. It there are fewer that this the window will be lengthened. Larger values lead to smoother estiamtes at the cost of reduced temporal precision.
	 * @param minWindow - the minimum size of the window that will be selected. In times of very high case numbers a minimum window prevents spurious outliers from completely dominating that days estimate of \(R_t\).
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator selectAdaptiveWindow(double incidenceSum, int minWindow) {
		if (incidenceSum <= 0) throw new ParameterOutOfRangeException("incidenceSum limit must be > 0");
		if (minWindow-1 > this.maxTau) throw new ParameterOutOfRangeException("The minWindow must be smaller than maxWindow (from the constructor)");
		this.posteriorSelectionStrategy = (res) -> res.selectWindowByCases(incidenceSum, minWindow);
		return this;
	}
	
	/**
	 * The minimum uncertainty strategy aims to select a window for the posterior \(R_t\) estimate based on minimising a combined value for uncertainty. The uncertainty in this case is the 
	 * product of the (timeVsRt)th root of SD of the posterior estimate and the window size to the power of timeVsRt. <br><br>
	 * 
	 * $$window^{timeVsRt} * SD_{posterior}^\frac{1}{timeVsRt}$$
	 * 
	 * @param timeVsRt - the relative importance parameter. Values above one favour temporal certainty; values under one favour \(R_t\) estimate certainty 
	 * @param minWindow - the minimum size of the window that will be selected.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator selectMinimumUncertainty(double timeVsRt, int minWindow) {
		if (minWindow-1 > this.maxTau) throw new ParameterOutOfRangeException("The minWindow must be smaller than maxWindow (from the constructor)");
		this.posteriorSelectionStrategy = (res) -> res.selectWindowByMinimumUncertainty(timeVsRt, minWindow);
		return this;
	}
	
	/**
	 * The mixture combination strategy for selecting a posterior \(R_t\) involves finding all the estimates that have a window that span the estimate date (this may involve longer windowed estimate from
	 * dates in the future) and estimating a gamma distribution with the same mean and SD as the mixture of all the estimates. This combined estimate involves estimates of all different windows. 
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod
	public CoriEstimator selectMixtureCombination() {
		this.posteriorSelectionStrategy = (res) -> res.selectMixtureCombination();
		return this;
	}
	
	/**
	 * The mixture combination strategy for selecting a posterior \(R_t\) involves finding all the estimates that have a window that span the estimate date (this may involve longer windowed estimate from
	 * dates in the future) and estimating a gamma distribution with the same mean and SD as the mixture of all the estimates. This combined estimate involves estimates of all different windows.
	 * @param weights - the weighting of the different windows according to their length as a vector which much be the same length as the longest window, maxTau.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod
	public CoriEstimator selectWeightedMixtureCombination(RNumericVector weights) {
		double[] wts = weights.javaPrimitive(0);
		if (wts.length != maxTau) throw new ParameterOutOfRangeException("The weights array must be the same length as the maximum window: "+maxTau);
		this.posteriorSelectionStrategy = (res) -> res.selectMixtureCombination(
				(datedRt) -> wts[datedRt.getWindow()-1]
				);
		return this;
	}
	
	/**
	 * This sets the selection strategy for prior estimate of the reproduction number to be the fixed R0 value provided to the CoriEstimator::new constructor. All estimates use the same prior
	 * and this is equivalent to the strategy implemented in EpiEstim (which has a default value of 5 for the mean and 5 for the SD of the R0 prior). This is an uninformed prior strategy.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator withDefaultPrior() {
		this.priorSelectionStrategy = (res) -> defaultPriors(res.dateValue(), res.getInfProfId());
		return this;
	}
	
	/**
	 * This sets the selection strategy for prior estimate of the reproduction number to be fixed value. All estimates use the same value as an uninformed prior, but enables us to compare the
	 * impact of different priors.
	 * @param mean - the mean of the fixed prior \(R_t\)
	 * @param sd - the sd  of the fixed prior \(R_t\) 
	 * @return The estimator itself (a fluent method)
	 */	
	@RMethod 
	public CoriEstimator withFixedPrior(double mean, double sd) {
		if (mean < 0) throw new ParameterOutOfRangeException("The mean of the prior Rt must be positive");
		if (sd <= 0) throw new ParameterOutOfRangeException("The SD of the prior Rt must be greater than zero");
		this.priorSelectionStrategy = (res) -> fixedPriors(mean,sd,res.dateValue(), res.getInfProfId());
		return this;
	}
	
	/**
	 * This sets the prior to be the posterior for the previous timestep, but with the SD of the previous posterior widened to allow the equivalent of a random walk. The widening factor defines how much 
	 * change we expect to see in each time step
	 * @param factor - must be above 1 (or else the estimate will become locked on a specific value)
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator withAdaptivePrior(double factor) {
		if (factor<1) throw new ParameterOutOfRangeException("The factor must be greater than or equal to 1"); 
		this.priorSelectionStrategy = (res) -> res.adaptiveNextPrior(factor);
		return this;
	}
	
	private CoriEstimator collectFirst() {
		this.combiningStrategy = (l) -> {
			if (l.posteriorsForProfiles.size() > 1) throw new RuntimeException("If there are more than one infectivity profiles a combining strategy must be defined");
			return l.posteriorsForProfiles.get(0);
		};
		return this;
	}
	
	/**
	 * This configures how the estimator combines multiple estimates of \(R_t\) using different infectivity profiles. In this case the estimator will 
	 * randomly sample from the selected posterior gamma distribution produced for each infectivity profile. These samples are combined and a mean, sd and empirical 
	 * quantiles calculated.    
	 * @param sampleSize - the number of samples to draw for each posterior estimate of \(R_t\). This number will be multiplied by the number of infectivity profile estimates and the number time points, so should not be set too high or it will lead to memory consumption issues. This is the equivalent of the $n2$ parameter in EpiEstim 
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator collectResampledQuantiles(int sampleSize) {
		this.combiningStrategy = (cere) -> GammaParameters.resamplingCombination(cere.posteriorsForProfiles, sampleSize);
		this.summarise = true;
		return this;
	}
	
	/**
	 * This configures how the estimator combines multiple estimates of \(R_t\) using different infectivity profiles. In this case it calculates a mixture distribution 
	 * of all the posterior estimates, and combined values for mean, and standard deviation, of the mixture plus the quantiles using a brent solver.  
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator collectMixtureQuantiles() {
		this.combiningStrategy = (cere) -> GammaParameters.mixtureDistribution(cere.posteriorsForProfiles);
		this.summarise = true;
		return this;
	}
	
	/**
	 * This configures how the estimator combines multiple estimates of \(R_t\) using different infectivity profiles. In this case it calculates a mixture distribution 
	 * of all the posterior estimates, and then estimates a single gamma distribution with the same mean and SD as the mixture. This will only be a valid approximation
	 * if all the posterior estimates are close to each other, which will be the case if the infectivitiy profiles are fairly similar. It shoudl tend to produce estimates
	 * with less confidence. It is however the quickest method.  
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator collectMixtureApproximation() {
		this.combiningStrategy = (cere) -> GammaParameters.mixtureApproximation(cere.posteriorsForProfiles);
		this.summarise = true;
		return this;
	}
	
	/**
	 * Configure the estimator to output full details for all infectivity profiles.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator detailedOutput() {
		this.selectAll();
		this.combiningStrategy = null;
		this.summarise = false;
		return this;
	}
	
	/**
	 * Configure the estimator with an infectivity profile. At least one profile must be added
	 * @param infectivityProfile - A numeric vector of the daily discrete probability of a secondary infection event occurring given an infected individual. It is assumed that this vector starts 
	 * at time zero with a probability of zero. 
	 * @param replace -  replace the existing infectivity profile with this new one? if false then the profile is added to a list.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod
	public CoriEstimator withInfectivityProfile(RNumericVector infectivityProfile, boolean replace) {
		/* Messing with the settings at this point turns out to be a bad idea
		 * because we assume defniition has occurred in some ordered way.
		if (infProf.size() == 0) {
			// set the do nothing collector if there is only one;
			this.collectFirst();
		} else {
			if (!this.summarise) {
				// set the mixture quantiles as the default collector for multiple
				// infectivity profiles. This will typically be overridden.
				this.collectMixtureQuantiles();
				this.summarise=true;
			}
			// if summarise is already true collector will have been set.
		}
		*/
		if(replace) infProf = new ArrayList<>();
		infProf.add(new InfectivityProfile(infectivityProfile, infProf.size()));
		return this;
	}
	
	/**
	 * Configure the estimator with an infectivity profile by specifying it as a matrix. This will replace all previously set profiles.
	 * @param infectivityProfiles - A numeric array of the daily discrete probability of a secondary infection event occurring given an infected individual. It is assumed that this vector starts 
	 * at time zero with a probability of zero. 
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod
	public CoriEstimator withInfectivityProfileMatrix(RNumericArray infectivityProfiles) {
		infProf = new ArrayList<>();
		if (infectivityProfiles.getDimensionality() != 2) throw new ParameterOutOfRangeException("Infectivity profiles should be defined in a matrix");
		try {
			infectivityProfiles.get().forEach(a -> {
				this.withInfectivityProfile(a.getVector(), false);
			});
		} catch (ZeroDimensionalArrayException e) {
			e.printStackTrace();
		}
		return this;
	}
	
	/**
	 * If the timeseries is at the start then we can safely impute zero values prior to the start of the timeseries. If the timeseries starts in the middle of the epidemic we must infer the 
	 * infections prior to the beginning of the time series if we are to prevent sudden jumps in the estimates over the initial part of the timeseries.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator atStartOfTimeseries() {
		atStart = true;
		return this;
	}
	
	/**
	 * If the timeseries starts in the middle of the epidemic we must infer the 
	 * infections prior to the beginning of the time series if we are to prevent sudden jumps in the 
	 * estimates over the initial part of the timeseries. This is done using the supplied value of R0 to 
	 * infer a growth rate.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod 
	public CoriEstimator inMiddleOfTimeseries() {
		atStart = false;
		return this;
	}
	
	/** 
	 * Estimates \(R_t\) for a one or more timeseries of infections in an epidemic, where the supplied dataframe is grouped
	 * @param incidence - a dataframe containing at least 2 columns (one containing a date, the other a numeric count of cases). The dataframe may be grouped, in which case grouping is preserved, and each group is treated as a seperate time series 
	 * @param dateColName - the name of the date column
	 * @param incidenceColName - the name of the incidence column
	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
	 */
	@RMethod
	public RDataframe estimateRt(RDataframe incidence, String dateColName, String incidenceColName) {
		
		
		RDataframe out = incidence
			.groupModify((d,g) -> {
				try {
					return estimateRtSingle(d.withColIfAbsent("errors", RCharacter.NA), dateColName, incidenceColName);
				} catch (IncompleteTimeseriesException e) {
					return d.mergeWithCol("errors", RConverter.convert("Incomplete timeseries: "+e.getMessage()), (s1,s2) -> RFunctions.paste("; ",s1,s2));
				} catch (UnexpectedNaValueException e) {
					return d.mergeWithCol("errors", RConverter.convert("Unexpected NA value"), (s1,s2) -> RFunctions.paste("; ",s1,s2));
				}
			});
		System.gc();
		return out;
		
	}
	
	/**
	 * Estimates \(R_t\) for a single timeseries of infections in an epidemic.
	 * @param incidence - a dataframe containing at least 2 columns (one containing a date, the other a numeric count of cases). The dataframe may be grouped, in which case grouping is preserved 
	 * @param dateColName - the name of the date column
	 * @param incidenceColName - the name of the incidence column
	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
	 */
	@RMethod
	public RDataframe estimateRtSingle(RDataframe incidence, String dateColName, String incidenceColName) throws IncompleteTimeseriesException {
		
		if (!incidence.containsKey(dateColName)) throw new ParameterOutOfRangeException("The dataframe does not contain a column: "+dateColName);
		if (!incidence.getTypeOfColumn(dateColName).equals(RDate.class)) throw new ParameterOutOfRangeException("The date column: "+dateColName+" is not of type RDate");
		if (!incidence.containsKey(incidenceColName)) throw new ParameterOutOfRangeException("The dataframe does not contain a column: "+incidenceColName);
		if (!incidence.getTypeOfColumn(incidenceColName).equals(RNumeric.class)) {
			if (incidence.getTypeOfColumn(incidenceColName).equals(RInteger.class)) {
				incidence.mutate(incidenceColName, RInteger.class, i -> RFunctions.asNumeric(i));
			} else {
				throw new ParameterOutOfRangeException("The incidence column: "+incidenceColName+" is not of type RNumeric (or RInteger)");
			}
		}
		if (!RFunctions.all(RFunctions::isFinite, incidence.pull(incidenceColName,RNumericVector.class))) {
			throw new ParameterOutOfRangeException("The incidence column contains non-finite values");
		}
		
		CoriTimeseries rtSeries = new CoriTimeseries(this); //, initialGrowth);
		
		try {
			incidence
				.rename("date", dateColName)
				.rename("I", incidenceColName)
				.stream(TimeseriesEntry.Incidence.class)
				.forEach(x -> rtSeries.add(new CoriTimeseriesEntry(x, rtSeries)));
		} catch (UnconvertableTypeException e) {
			throw new RuntimeException(e);
		}
		
		CoriEstimationResult allResults = estimateRtSingle(rtSeries);
		CoriEstimationSummary filteredResults = allResults.selectPosterior();
		
		//Maybe stream here and combine later. 
		
		if (!summarise) {
			// TODO: sometimes it woudl be good to apply a posteriorSelectionStrategy independently of a bootstrap collection.
			// at the moment collect bootstraps applies the posterior selection strategy first
			// doing this independently would mena the following has to change as woudl be a summay entry rather than a result entry
			
			RDataframe out = filteredResults.stream().collect(
				CoriEstimationSummaryEntry.collector(dateColName,incidenceColName,epiEstimMode));
			return out;
		
		} else { 
		
			RDataframe out2 = filteredResults.stream().collect(
				CoriEstimationSummaryEntry.summaryCollector(dateColName,incidenceColName,epiEstimMode ));
			return(out2);
		}
		
	}
	
	
	
	/**
	 * This converts the input timeseries into a EstimateResult time series. This is the native Java methd in the api. 
	 * @param rtSeries - a CoriTimeseries 
	 * @return a CoriEstimateResult.
	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
	 */
	public CoriEstimationResult estimateRtSingle(CoriTimeseries rtSeries) throws IncompleteTimeseriesException {
		if (infProf.isEmpty()) throw new RuntimeException("Needs an infectivity profile to be defined before we can estimate Rt");
		
		rtSeries.checkComplete();
		
		CoriEstimationResult allResults = new CoriEstimationResult(this);
		infProf.parallelStream().map(prof -> {
			
			CoriEstimationResult singleResult = new CoriEstimationResult(this); 
			CoriTimeseries rtWorking = new CoriTimeseries(rtSeries);
			GrowthRateUtils conv = new GrowthRateUtils(prof);
			//TODO: this is a bit of a mess. Probably should embed the InfectionProfile into rtWorking
			//which woudl trigger the inferStart
			//and cache the estimation componenets
			double initialGrowth = conv.RToGrowth(rZero.convert().getMean());
			rtWorking.inferStart(Optional.ofNullable(atStart ? 0 : null), prof.length(), initialGrowth);
			rtWorking.estimateForProfile(prof.profile());
			CoriEstimationResultEntry last = null;
			Optional<CoriTimeseriesEntry> tsEntry = rtWorking.start();
			while (tsEntry.isPresent()) {
				
				List<DatedRtEstimate> priors;
				if (last == null) {
					priors = defaultPriors(tsEntry.get().dateValue(), prof.getId());
				} else {
					priors = priorSelectionStrategy.apply(last);
				}
				List<DatedRtEstimate> results = tsEntry.get().results(priors);
				
				CoriEstimationResultEntry current = new CoriEstimationResultEntry(prof.getId(), tsEntry.get(), results, singleResult);
				singleResult.add(current);
				last = current;
				tsEntry = tsEntry.get().next();
			};
			
			return singleResult;
		
		}).forEach(allResults::addAll);
		System.gc();
		return allResults;
 	}

	/** 
	 * Estimates \(R_t\) for a one or more timeseries of infections in an epidemic, where the supplied dataframe is grouped
	 * @param rates - a dataframe containing at least 2 columns (one containing a date, the other a numeric poission rate of cases). The dataframe may be grouped, in which case grouping is preserved, and each group is treated as a seperate time series 
	 * @param dateColName - the name of the date column
	 * @param rateColName - the name of the rate column
	 * @param samplesPerProfile - the number of times to bootstrap per infectivity profile
	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
	 */
	@RMethod
	public RDataframe estimateRtFromRates(RDataframe rates, String dateColName, String rateColName, int samplesPerProfile) {
		
		
		RDataframe out = rates
			.groupModify((d,g) -> {
				try {
					return estimateRtFromSingleRate(d.withColIfAbsent("errors", RCharacter.NA), dateColName, rateColName, samplesPerProfile);
				} catch (IncompleteTimeseriesException e) {
					return d.mergeWithCol("errors", RConverter.convert("Incomplete timeseries: "+e.getMessage()), (s1,s2) -> RFunctions.paste("; ",s1,s2));
				} catch (UnexpectedNaValueException e) {
					return d.mergeWithCol("errors", RConverter.convert("Unexpected NA value"), (s1,s2) -> RFunctions.paste("; ",s1,s2));
				}
			});
		System.gc();
		return out;
		
	}
	
	/**
	 * Estimates \(R_t\) for a single timeseries of infections in an epidemic.
	 * @param rate - a dataframe containing at least 2 columns (one containing a date, the other a numeric estimate of possion rate of cases). The dataframe may be grouped, in which case grouping is preserved 
	 * @param dateColName - the name of the date column
	 * @param rateColName - the name of the estimated poisson rate
	 * @param samplesPerProfile - the number of times to bootstrap per infectivity profile
	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
	 */
	@RMethod
	public RDataframe estimateRtFromSingleRate(RDataframe rate, String dateColName, String rateColName, int samplesPerProfile) throws IncompleteTimeseriesException {
		
		if (!rate.containsKey(dateColName)) throw new ParameterOutOfRangeException("The dataframe does not contain a column: "+dateColName);
		if (!rate.getTypeOfColumn(dateColName).equals(RDate.class)) throw new ParameterOutOfRangeException("The date column: "+dateColName+" is not of type RDate");
		if (!rate.containsKey(rateColName)) throw new ParameterOutOfRangeException("The dataframe does not contain a column: "+rateColName);
		if (!rate.getTypeOfColumn(rateColName).equals(RNumeric.class)) {
			if (!rate.getTypeOfColumn(rateColName).equals(RInteger.class)) {
				rate.mutate(rateColName, RInteger.class, i -> RFunctions.asNumeric(i));
			} else {
				throw new ParameterOutOfRangeException("The incidence column: "+rateColName+" is not of type RNumeric (or RInteger)");
			}
		}
		
		CoriBootstrappedTimeseries rtSeries = new CoriBootstrappedTimeseries(this); //, initialGrowth);
		
		try {
			rate
				.rename("date", dateColName)
				.rename("pois", rateColName)
				.stream(TimeseriesEntry.Rate.class)
				.forEach(x -> rtSeries.add(new CoriBootstrappedTimeseriesEntry(x, rtSeries)));
		} catch (UnconvertableTypeException e) {
			throw new RuntimeException(e);
		}
		
		CoriEstimationResult allResults = estimateRtBootstrapped(rtSeries, samplesPerProfile);
		CoriEstimationSummary summaryResults = allResults.selectPosterior();
		
		//Maybe stream here and combine later. 
		
		if (!summarise) {
			RDataframe out = summaryResults.stream().collect(
				CoriEstimationSummaryEntry.poissonCollector(dateColName,rateColName,epiEstimMode));
			return out;
		
		} else { 
		
			RDataframe out2 = summaryResults.stream().collect(
				CoriEstimationSummaryEntry.poissonSummaryCollector(dateColName,rateColName,epiEstimMode ));
			return out2;
		}
		
	}
	
//	private static int lcm(int number1, int number2) {
//	    if (number1 == 0 || number2 == 0) {
//	        return 0;
//	    }
//	    int absNumber1 = Math.abs(number1);
//	    int absNumber2 = Math.abs(number2);
//	    int absHigherNumber = Math.max(absNumber1, absNumber2);
//	    int absLowerNumber = Math.min(absNumber1, absNumber2);
//	    int lcm = absHigherNumber;
//	    while (lcm % absLowerNumber != 0) {
//	        lcm += absHigherNumber;
//	    }
//	    return lcm;
//	}
	
	/**
	 * This converts the input timeseries into a EstimateResult time series. This is the native Java methd in the api. 
	 * @param rtSeries - a CoriTimeseries
	 * @param samplesPerProfile - the number of times to bootstrap per infectivity profile 
	 * @return a CoriEstimateResult.
	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
	 */
	public CoriEstimationResult estimateRtBootstrapped(CoriBootstrappedTimeseries rtSeries, int samplesPerProfile) throws IncompleteTimeseriesException {
		if (infProf.isEmpty()) throw new RuntimeException("Needs an infectivity profile to be defined before we can estimate Rt");
		
		rtSeries.checkComplete();
		
		// Make sure infProf has eough entries to allow for bootstrapping 
		// of timeseries
//		int infProfLen = lcm(samples, infProf.size());
//		int repeats = infProfLen / (infProf.size());
//		if (infProfLen > 10*samples) infProfLen = 10*samples;
		List<InfectivityProfile> infProfRep = new ArrayList<>();
		for (int i=0; i<samplesPerProfile; i++) {
			infProfRep.addAll(infProf); 
		}
				
		CoriEstimationResult allResults = new CoriEstimationResult(this);
		infProfRep.parallelStream().map(prof -> {
			
			CoriEstimationResult singleResult = new CoriEstimationResult(this);
			CoriTimeseries rtWorking = rtSeries.sample();
			GrowthRateUtils conv = new GrowthRateUtils(prof);
			double initialGrowth = conv.RToGrowth(rZero.convert().getMean());
			rtWorking.inferStart(Optional.ofNullable(atStart ? 0 : null), prof.length(), initialGrowth);
			rtWorking.estimateForProfile(prof.profile());
			CoriEstimationResultEntry last = null;
			Optional<CoriTimeseriesEntry> tsEntry = rtWorking.start();
			while (tsEntry.isPresent()) {
				
				List<DatedRtEstimate> priors;
				if (last == null) {
					priors = defaultPriors(tsEntry.get().dateValue(), prof.getId());
				} else {
					priors = priorSelectionStrategy.apply(last);
				}
				List<DatedRtEstimate> results = tsEntry.get().results(priors);
				
				CoriEstimationResultEntry current = new CoriEstimationResultEntry(prof.getId(), tsEntry.get(), results, singleResult);
				singleResult.add(current);
				last = current;
				tsEntry = tsEntry.get().next();
			};
			
			return singleResult;
		
		}).forEach(allResults::addAll);
		System.gc();
		return allResults;
 	}
	
	
	private List<DatedRtEstimate> defaultPriors(LocalDate dateValue, int profileId) {
		return IntStream.range(0,maxTau).mapToObj(i -> defaultPrior(i, dateValue, profileId)).collect(Collectors.toList());
	}
	
	private List<DatedRtEstimate> fixedPriors(double mean, double sd, LocalDate dateValue, int profileId) {
		return IntStream.range(0,maxTau).mapToObj(i -> fixedPrior(mean,sd,i, dateValue, profileId)).collect(Collectors.toList());
	}
	
	protected DatedRtEstimate defaultPrior(int tau, LocalDate dateValue, int profileId) {
		return rZero
				.withDate(tau, dateValue, 0D, 0D)
				.withProfileId(profileId);
	}
	
	protected static DatedRtEstimate fixedPrior(double mean, double sd, int tau, LocalDate dateValue, int profileId) {
		return new GammaMoments(mean,sd).convert()
				.withDate(tau, dateValue, 0D, 0D)
				.withProfileId(profileId);
	}

	protected int getMaxTau() {
		return maxTau;
	}

	
	/**
	 * Get the initial R0 as a GammaParameters object. native Java Api only.
	 * @return R0 gamma distribution
	 */
	public GammaParameters getRZero() {
		return this.rZero;
	}

	/**
	 * Get the maximal window value this estimator supports. native Java Api only.
	 * @return the maximum window size
	 */
	public int getMaxWindow() {
		return this.maxTau+1;
	}
	
}
