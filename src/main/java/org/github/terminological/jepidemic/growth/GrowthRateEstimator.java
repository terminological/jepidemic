package org.github.terminological.jepidemic.growth;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.IntStream;

import org.github.terminological.jepidemic.IncompleteTimeseriesException;
import org.github.terminological.jepidemic.InfectivityProfile;
import org.github.terminological.jepidemic.ParameterOutOfRangeException;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.distributions.BetaPrimeDistribution;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;
import org.github.terminological.jepidemic.distributions.GrowthRateDistribution;
import org.github.terminological.jepidemic.distributions.NegativeBinomialDistribution;
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
 * <br>A. GrowthRate, N. M. Ferguson, C. Fraser, and S. Cauchemez, ‘A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics’, Am J Epidemiol, vol. 178, no. 9, pp. 1505–1512, Nov. 2013, doi: 10.1093/aje/kwt133.<br>
 * 
 * <br> and the implementation of it in the R Library EpiEstim. It assumes that infections (or other type of observation) are a 
 * Poisson distributed process and that the reproduction number is a gamma distributed quantity. The estimates of \(R_t\) may be 
 * aggregated over a specific "window" to provide a more stable estimate.  
 */
@RClass
public class GrowthRateEstimator {

	static Logger log = LoggerFactory.getLogger(GrowthRateEstimator.class);
	
	List<InfectivityProfile> infProf = new ArrayList<>();
	private ExtendedGammaDistribution incidenceZero;
	private int maxTau = 14;
	private int minTau = 7;
	private Map<Integer,Double> weightsForTau;
	private double[] quantiles = new double[] {0.025,0.05,0.25,0.5,0.75,0.95,0.975};
		
	Strategy.PriorSelection priorSelectionStrategy;
	Strategy.PosteriorFiltering posteriorFilteringStrategy;
	Strategy.CombiningStrategy.ForContinuous<ExtendedGammaDistribution> combiningLambdaStrategy;
	Strategy.CombiningStrategy.ForContinuous<GrowthRateDistribution> combiningGrowthRateStrategy;
	Strategy.CombiningStrategy.ForContinuous<BetaPrimeDistribution> combiningRtStrategy;
	Strategy.CombiningStrategy.ForDiscrete<NegativeBinomialDistribution> combiningIncidenceStrategy;

	protected boolean doRtEstimation() {
		boolean doRt = !infProf.isEmpty();
		return doRt;
	}
	
	// ===== Package wide config accessors ====
	
	protected Strategy.PriorSelection getPriorSelectionStrategy() {
		return priorSelectionStrategy;
	}
	
	protected Strategy.PosteriorFiltering getPosteriorFilteringStrategy() {
		return posteriorFilteringStrategy;
	}
	
	protected Strategy.CombiningStrategy.ForContinuous<ExtendedGammaDistribution> getCombiningLambdaStrategy() {
		return combiningLambdaStrategy;
	}

	protected Strategy.CombiningStrategy.ForContinuous<GrowthRateDistribution> getCombiningGrowthRateStrategy() {
		return combiningGrowthRateStrategy;
	}

	protected Strategy.CombiningStrategy.ForDiscrete<NegativeBinomialDistribution> getCombiningIncidenceStrategy() {
		return combiningIncidenceStrategy;
	}
	
	protected Strategy.CombiningStrategy.ForContinuous<BetaPrimeDistribution> getCombiningRtStrategy() {
		return combiningRtStrategy;
	}
	
	protected int getMaxTau() {
		return maxTau;
	}
	
	protected Optional<ExtendedGammaDistribution> getInitialIncidence() {
		return Optional.ofNullable(this.incidenceZero);
	}

	protected int getMaxWindow() {
		return this.maxTau+1;
	}

	protected int getMinTau() {
		return minTau;
	}

	protected double[] getQuantiles() {
		return quantiles;
	}

	protected Optional<Integer> getMaxInfProf() {
		return this.infProf.stream().map(i -> i.length()).reduce(Math::max);
	}

	/** 
	 * Produces a GrowthRateEstimator that is configured to behave in a sensible fashion, with sane defaults. 
	 * In that it estimates an initial poisson rate prior from the first entry and returns a fixed window size and aggregates uncertain SI 
	 * distributions with a random resampling strategy<br><br>
	 * 
	 * 
	 * @param infectivityProfiles - a numeric matrix with the discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)
	 * @param minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.
	 * @param minWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount .
	 * @return a default growth rate estimator
	 */
	@RMethod 
	public GrowthRateEstimator(int minWindow, int maxWindow) {
		
		if (maxWindow <= 1) throw new ParameterOutOfRangeException("The maxWindow parameter must be 2 or more. Probably it should be >=7 unless working with very large number of timeseries");
		this.maxTau = maxWindow;
		this.minTau = minWindow;
		this.withSaneDefaults();
	}
	
	/** 
	 * Produces a GrowthRateEstimator that will estimate growth rate and $R_t$, with sane defaults. The $R_t$ estimation is 
	 * done using bootstrapping, and combining the results using a mixture distribution. This will be relatively slow.
	 *  
	 * @param minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.
	 * @param minWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount.
	 * @param infectivityProfiles - a numeric matrix with each column containing a discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)
	 * @return a strict growth rate and $R_t$ estimator that accounts for uncertainty in the infectivity profile.  
	 */
	@RMethod 
	public static GrowthRateEstimator strictGrowthRateAndRtEstimator(int minWindow, int maxWindow, RNumericArray infectivityProfiles) {
		GrowthRateEstimator out = new GrowthRateEstimator(minWindow, maxWindow)
				.withInfectivityProfileMatrix(infectivityProfiles);
		return out;
	}
	
	/** 
	 * Produces a GrowthRateEstimator that will estimate growth rate and $R_t$, with sane defaults. The $R_t$ estimation is 
	 * done using a single infectivity profile and does not account for infectivity profile uncertainty. This is comparatively fast.<br><br>
	 * 
	 * 
	 * @param minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.
	 * @param minWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount.
	 * @param infectivityProfiles - a numeric matrix with each column containing a discrete probability distribution of the serial interval with \(P_{(t=0)} = 0\)
	 * @return a growth rate and $R_t$ estimator that does not account for uncertainty in the infectivity profile.  
	 */
	@RMethod 
	public static GrowthRateEstimator defaultGrowthRateAndRtEstimator(int minWindow, int maxWindow, double initialIncidence) {
		GrowthRateEstimator out = new GrowthRateEstimator(minWindow, maxWindow);
		out.withInitialIncidence(initialIncidence);
		return out;
	}
	
	/** 
	 * Produces a GrowthRateEstimator that will estimate growth rate only, with sane defaults, unless an infectivity profile is lated added<br><br>
	 * 
	 * @param minWindow - the smallest amount of data considered in estimating the growth rate - 4 is a ssensible default unless you want bang up to date estimate.
	 * @param minWindow - the largest amount of data considered in estimating the growth rate - 14 is a sensible amount.
	 * @return a growth rate and $R_t$ estimator that does not account for uncertainty in the infectivity profile.  
	 */
	@RMethod 
	public static GrowthRateEstimator defaultGrowthRateEstimator(int minWindow, int maxWindow) {
		GrowthRateEstimator out = new GrowthRateEstimator(minWindow, maxWindow);
		return out;
	}
	
	/** 
	 * Produces a GrowthRateEstimator that will estimate growth rate only, with sane defaults, and estimating using windows between 4 and 14 days<br><br>
	 * 
	 * @return a growth rate and $R_t$ estimator that does not account for uncertainty in the infectivity profile.  
	 */
	@RMethod 
	public static GrowthRateEstimator basicGrowthRateEstimator() {
		GrowthRateEstimator out = new GrowthRateEstimator(4, 14);
		return out;
	}
	
	/**
	 * Allows you to change the quantiles returned by the estimator. The defaults are c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
	 * @param quantiles the desired quantiles of the output
	 * @return a fluent method
	 */
	@RMethod 
	public GrowthRateEstimator withOutputQuantiles(RNumericVector quantiles) {
		this.quantiles = quantiles.javaPrimitive(Double.NaN);
		return this;
	}
	
	/**
	 * Provides a sensible baseline configuration for the estimator that involves estimating the poisson rate using the 
	 * previous timestep's posterior with standard deviation scaled by a factor determined by the t-2 value. Multiple posterior
	 * estimates are combined using a weighted average based on a discretised normal distributed weight with mean of 7 and SD of 4.
	 * @return a fluent method
	 */
	@RMethod 
	public GrowthRateEstimator withSaneDefaults() {
		this.withGaussianEstimateWeighting(7, 4);
		this.priorIncidenceFromScaledPreviousPosterior();
		this.useAllPosteriorEstimates();
		this.combineEstimatesWithWeightedMixture();
		return this;
	}
	
	/**
	 * Allows configuration of the initial incidence. If this is omitted the initial incidence is guessed from the first entry in the timeseries
	 * @param incidence
	 * @return a fluent method
	 */
	@RMethod
	public GrowthRateEstimator withInitialIncidence(double incidence) {
		if (incidence < 0) throw new ParameterOutOfRangeException("The mean of the prior incidence must be positive");
		this.incidenceZero = ExtendedGammaDistribution.fromMoments(incidence, incidence);
		return this;
	}
	
	// ===== Combining configuration =========
	
	/**
	 * Allows configuration of the weighting of estimate combination for the weighted mixture combination. 
	 * This uses a discretised normal distribution to generate values for the weighting. The weighting is 
	 * applied to each supported estimate determined by the data window of that estimate.  The default value
	 * is a mean of 7 and SD of 4
	 * @param meanTau the mean of the normal distribution of weights
	 * @param sdTau the sd of the normal distribution of weights
	 * @return a fluent method
	 */
	@RMethod 
	public GrowthRateEstimator withGaussianEstimateWeighting(double meanTau, double sdTau) {
		ExtendedGammaDistribution nd = ExtendedGammaDistribution.fromMoments(meanTau, sdTau); 
		this.weightsForTau = new HashMap<>();
		IntStream.range(minTau,maxTau).forEach(i -> weightsForTau.put(i, nd.probability((double) i-1,(double) i)));
		return this;
	}
	
	/**
	 * Allows configuration of the weighting of estimate combination for the weighted mixture combination. 
	 * This assumes the weights are given as a vector where the first entry is the weight for an estimate using a 
	 * data window of 1. The vector should be the same length as the maxWindow parameter used in configuring the 
	 * GrowthRateEstimator or 14 if the defaults were used.
	 * @param weights - a vector of weights for each data window
	 * @return a fluent method
	 */
	@RMethod 
	public GrowthRateEstimator withEstimateWeighting(RNumericVector weights) {
		this.weightsForTau = new HashMap<>();
		IntStream.range(1,maxTau).forEach(i -> weightsForTau.put(i, weights.get(i-1).javaPrimitive()));
		return this;
	}
	
	/**
	 * Sets the estimator to use a weighted mixture to combine estimates based on the length of different data windows, using weights defined 
	 * by withEstimateWeighting() or withGaussianEstimateWeighting(). If none is set the weighting defaults to a discretised normal
	 * distribution with mean of 7 an SD of 4.
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator combineEstimatesWithWeightedMixture() {
		if (this.weightsForTau == null) this.withGaussianEstimateWeighting(7, 4);
		this.combiningLambdaStrategy = Strategy.CombiningStrategy.ForContinuous.weightedMixture(weightsForTau, 0, 1000000);
		this.combiningGrowthRateStrategy = Strategy.CombiningStrategy.ForContinuous.weightedMixture(weightsForTau,-3, +3);
		this.combiningIncidenceStrategy = Strategy.CombiningStrategy.ForDiscrete.weightedMixture(weightsForTau);
		this.combiningRtStrategy = Strategy.CombiningStrategy.ForContinuous.weightedMixture(weightsForTau, 0, 10);
		return this;
	}
	
	// private method unsed when only using one window
	protected GrowthRateEstimator combineUsingOnlyFirstEstimate() {
		if (this.weightsForTau == null) this.withGaussianEstimateWeighting(7, 4);
		this.combiningLambdaStrategy = Strategy.CombiningStrategy.ForContinuous.collectFirst();
		this.combiningGrowthRateStrategy = Strategy.CombiningStrategy.ForContinuous.collectFirst();
		this.combiningIncidenceStrategy = Strategy.CombiningStrategy.ForDiscrete.collectFirst();
		this.combiningRtStrategy = Strategy.CombiningStrategy.ForContinuous.collectFirst();
		return this;
	}
	
	// ===== Prior selection configuration =========
	
	
	/**
	 * Select a prior estimate for the poisson rate by taking the previous posterior and multiplying 
	 * the SD by a constant factor
	 * @param multiplyingSDby - factor to multiply SD by. Must be greater than or equal to 1.
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator priorIncidenceFromPreviousPosterior(double multiplyingSDby) {
		this.priorSelectionStrategy = Strategy.PriorSelection.posteriorAsPrior(multiplyingSDby);
		return this;
	}
	
	/**
	 * Select a prior estimate for the poisson rate by taking the previous posterior and multiplying 
	 * the SD by a factor determined by the difference between the t-1 and t-2 incidence. This allows for 
	 * rapid rate of change to be detected and less weight given to the prior in this situation 
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator priorIncidenceFromScaledPreviousPosterior() {
		this.priorSelectionStrategy = Strategy.PriorSelection.posteriorAsPrior();
		return this;
	}
	
	/**
	 * Select a prior estimate for the poisson rate by taking the previous posterior and increasing it by  
	 * the exponential previous posterior growth rate and multiplying the resulting disctribution by a 
	 * constant factor. This is a mechanistic prior. 
	 * @param multiplyingSDby - factor to multiply SD by. Must be greater than or equal to 1.
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator priorIncidenceFromPosteriorAndGrowthRate(double multiplyingSDby) {
		this.priorSelectionStrategy = Strategy.PriorSelection.mechanisticPrior(multiplyingSDby);
		return this;
	}
	
	/**
	 * Select a prior estimate for the poisson rate by taking the previous posterior and increasing it by  
	 * the exponential previous posterior growth rate and multiplying the resulting disctribution by a 
	 * factor determined by the change in t-1 and t-2 indidence. This is a mechanistic prior. 
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator priorIncidenceFromScaledPosteriorAndGrowthRate() {
		this.priorSelectionStrategy = Strategy.PriorSelection.combinedMechanisticPrior();
		return this;
	}
	
	// ===== Posterior selection =================
	
	/**
	 * The posterior estimates for each combination of data window and infectivity profile are all retained and combined
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator useAllPosteriorEstimates() {
		this.posteriorFilteringStrategy = Strategy.PosteriorFiltering.allowAll();
		return this;
	}
	
	/**
	 * The posterior estimates for only one data window, but all infectivity profiles, are retained and combined. N.b. the same effect
	 * can be achieved by setting the minWindow and maxWindow in the constructor.
	 * @param tau the window width to retain.
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator usePosteriorEstimatesFromOneWindow(int tau) {
		this.posteriorFilteringStrategy = Strategy.PosteriorFiltering.specificWindow(tau);
		this.combineUsingOnlyFirstEstimate();
		return this;
	}
	
	/**
	 * The posterior estimates for data windows that include a reasonable amount of data are retained and combined. This 
	 * suppresses noisy estimates made on a very few data points when incidence is low.
	 * @param sumIndidence - the minimum number of cases that must be observed in a data window before the estimate is ignored.
	 * @return a fluent method
	 */
	@RMethod public GrowthRateEstimator usePosteriorEstimatesWithEnoughData(int sumIncidence) {
		this.posteriorFilteringStrategy = Strategy.PosteriorFiltering.minIncidence(sumIncidence); //.minIncidence(30);
		return this;
	}
	
	// =========== RT configuration ============
	
	/**
	 * Configure the estimator with an infectivity profile. At least one profile must be added for an estimate of R_t to be made
	 * @param infectivityProfile - A numeric vector of the daily discrete probability of a secondary infection event occurring given an infected individual. It is assumed that this vector starts 
	 * at time zero with a probability of zero. 
	 * @param replace -  replace the existing infectivity profile with this new one? if false then the profile is added to a list.
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod
	public GrowthRateEstimator withInfectivityProfile(RNumericVector infectivityProfile, boolean replace) {
		if(replace) infProf = new ArrayList<>();
		infProf.add(new InfectivityProfile(infectivityProfile, infProf.size()));
		return this;
	}
	
	/**
	 * Configure the estimator with an infectivity profile by specifying it as a matrix. This will replace all previously set profiles.
	 * @param infectivityProfiles - A numeric array of the daily discrete probability of a secondary infection event occurring given an infected individual in each column. It is assumed that this vector starts 
	 * at time zero with a probability of zero. 
	 * @return The estimator itself (a fluent method)
	 */
	@RMethod
	public GrowthRateEstimator withInfectivityProfileMatrix(RNumericArray infectivityProfiles) {
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
	 * Estimates growth rate and, if configured with an infectivity profile, \(R_t\) for a one or more timeseries of infections in an epidemic. If the supplied dataframe is grouped it assumes ach grouping is a timeseries
	 * @param incidence - a dataframe containing at least 2 columns (one containing a date, the other a numeric count of cases). The dataframe may be grouped, in which case grouping is preserved, and each group is treated as a seperate time series 
	 * @param dateColName - the name of the date column
	 * @param incidenceColName - the name of the incidence column
	 * @return a dataframe, containing a summary of the estimates of incidence, growth rate and \(R_t\) with mean, sd, and quantiles for each day.
	 */
	@RMethod
	public RDataframe estimateGrowthRate(RDataframe incidence, String dateColName, String incidenceColName) {
		
		
		RDataframe out = incidence
			.groupModify((d,g) -> {
				try {
					return estimateGrowthRateSingle(d.withColIfAbsent("errors", RCharacter.NA), dateColName, incidenceColName);
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
	 * @param incidence - a dataframe containing at least 2 columns (one containing a date, the other a numeric count of cases).  
	 * @param dateColName - the name of the date column
	 * @param incidenceColName - the name of the incidence column
	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
	 */
	@RMethod
	public RDataframe estimateGrowthRateSingle(RDataframe incidence, String dateColName, String incidenceColName) throws IncompleteTimeseriesException {
		
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
		
		GrowthRateTimeseries rtSeries = new GrowthRateTimeseries(this); //, initialGrowth);
		
		try {
			incidence
				.rename("date", dateColName)
				.rename("I", incidenceColName)
				.stream(TimeseriesEntry.Incidence.class)
				.forEach(x -> rtSeries.add(new GrowthRateTimeseriesEntry(x, rtSeries)));
		} catch (UnconvertableTypeException e) {
			throw new RuntimeException(e);
		}
		
		rtSeries.forEach(s -> {
			s.getSummaryEstimates();
			s.cleanUp();
		});
		
		//Maybe stream here and combine later.
		if (this.doRtEstimation()) {
		
			RDataframe out = rtSeries.stream().collect(Formatter.rtAndGrowthCollector(dateColName, incidenceColName));
			return out;
		
		} else { 
		
			RDataframe out = rtSeries.stream().collect(Formatter.growthCollector(dateColName, incidenceColName));
			return out;
		}
		
	}
	
	
	
//	/**
//	 * This converts the input timeseries into a EstimateResult time series. This is the native Java methd in the api. 
//	 * @param rtSeries - a GrowthRateTimeseries 
//	 * @return a GrowthRateEstimateResult.
//	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
//	 */
//	public GrowthRateEstimationResult estimateRtSingle(GrowthRateTimeseries rtSeries) throws IncompleteTimeseriesException {
//		if (infProf.isEmpty()) throw new RuntimeException("Needs an infectivity profile to be defined before we can estimate Rt");
//		
//		rtSeries.checkComplete();
//		
//		GrowthRateEstimationResult allResults = new GrowthRateEstimationResult(this);
//		infProf.parallelStream().map(prof -> {
//			
//			GrowthRateEstimationResult singleResult = new GrowthRateEstimationResult(this); 
//			GrowthRateTimeseries rtWorking = new GrowthRateTimeseries(rtSeries);
//			GrowthRateUtils conv = new GrowthRateUtils(prof);
//			//TODO: this is a bit of a mess. Probably should embed the InfectionProfile into rtWorking
//			//which woudl trigger the inferStart
//			//and cache the estimation componenets
//			double initialGrowth = conv.RToGrowth(rZero.convert().getMean());
//			rtWorking.inferStart(Optional.ofNullable(atStart ? 0 : null), prof.length(), initialGrowth);
//			rtWorking.estimateForProfile(prof.profile());
//			GrowthRateEstimationResultEntry last = null;
//			Optional<GrowthRateTimeseriesEntry> tsEntry = rtWorking.start();
//			while (tsEntry.isPresent()) {
//				
//				List<DatedGammaEstimate> priors;
//				if (last == null) {
//					priors = defaultPriors(tsEntry.get().dateValue(), prof.getId());
//				} else {
//					priors = priorSelectionStrategy.apply(last);
//				}
//				List<DatedGammaEstimate> results = tsEntry.get().results(priors);
//				
//				GrowthRateEstimationResultEntry current = new GrowthRateEstimationResultEntry(prof.getId(), tsEntry.get(), results, singleResult);
//				singleResult.add(current);
//				last = current;
//				tsEntry = tsEntry.get().next();
//			};
//			
//			return singleResult;
//		
//		}).forEach(allResults::addAll);
//		System.gc();
//		return allResults;
// 	}
//
//	/** 
//	 * Estimates \(R_t\) for a one or more timeseries of infections in an epidemic, where the supplied dataframe is grouped
//	 * @param rates - a dataframe containing at least 2 columns (one containing a date, the other a numeric poission rate of cases). The dataframe may be grouped, in which case grouping is preserved, and each group is treated as a seperate time series 
//	 * @param dateColName - the name of the date column
//	 * @param rateColName - the name of the rate column
//	 * @param samplesPerProfile - the number of times to bootstrap per infectivity profile
//	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
//	 */
//	@RMethod
//	public RDataframe estimateRtFromRates(RDataframe rates, String dateColName, String rateColName, int samplesPerProfile) {
//		
//		
//		RDataframe out = rates
//			.groupModify((d,g) -> {
//				try {
//					return estimateRtFromSingleRate(d.withColIfAbsent("errors", RCharacter.NA), dateColName, rateColName, samplesPerProfile);
//				} catch (IncompleteTimeseriesException e) {
//					return d.mergeWithCol("errors", RConverter.convert("Incomplete timeseries: "+e.getMessage()), (s1,s2) -> RFunctions.paste("; ",s1,s2));
//				} catch (UnexpectedNaValueException e) {
//					return d.mergeWithCol("errors", RConverter.convert("Unexpected NA value"), (s1,s2) -> RFunctions.paste("; ",s1,s2));
//				}
//			});
//		System.gc();
//		return out;
//		
//	}
//	
//	/**
//	 * Estimates \(R_t\) for a single timeseries of infections in an epidemic.
//	 * @param rate - a dataframe containing at least 2 columns (one containing a date, the other a numeric estimate of possion rate of cases). The dataframe may be grouped, in which case grouping is preserved 
//	 * @param dateColName - the name of the date column
//	 * @param rateColName - the name of the estimated poisson rate
//	 * @param samplesPerProfile - the number of times to bootstrap per infectivity profile
//	 * @return a dataframe, containing either a summary of the estimates of \(R_t\) with mean, sd, and quantiles for each day, or a full breakdown of the component \(R_t\) estimates.
//	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
//	 */
//	@RMethod
//	public RDataframe estimateRtFromSingleRate(RDataframe rate, String dateColName, String rateColName, int samplesPerProfile) throws IncompleteTimeseriesException {
//		
//		if (!rate.containsKey(dateColName)) throw new ParameterOutOfRangeException("The dataframe does not contain a column: "+dateColName);
//		if (!rate.getTypeOfColumn(dateColName).equals(RDate.class)) throw new ParameterOutOfRangeException("The date column: "+dateColName+" is not of type RDate");
//		if (!rate.containsKey(rateColName)) throw new ParameterOutOfRangeException("The dataframe does not contain a column: "+rateColName);
//		if (!rate.getTypeOfColumn(rateColName).equals(RNumeric.class)) {
//			if (!rate.getTypeOfColumn(rateColName).equals(RInteger.class)) {
//				rate.mutate(rateColName, RInteger.class, i -> RFunctions.asNumeric(i));
//			} else {
//				throw new ParameterOutOfRangeException("The incidence column: "+rateColName+" is not of type RNumeric (or RInteger)");
//			}
//		}
//		
//		GrowthRateBootstrappedTimeseries rtSeries = new GrowthRateBootstrappedTimeseries(this); //, initialGrowth);
//		
//		try {
//			rate
//				.rename("date", dateColName)
//				.rename("pois", rateColName)
//				.stream(TimeseriesEntry.Rate.class)
//				.forEach(x -> rtSeries.add(new GrowthRateBootstrappedTimeseriesEntry(x, rtSeries)));
//		} catch (UnconvertableTypeException e) {
//			throw new RuntimeException(e);
//		}
//		
//		GrowthRateEstimationResult allResults = estimateRtBootstrapped(rtSeries, samplesPerProfile);
//		GrowthRateEstimationSummary summaryResults = allResults.selectPosterior();
//		
//		//Maybe stream here and combine later. 
//		
//		if (!summarise) {
//			RDataframe out = summaryResults.stream().collect(
//				GrowthRateEstimationSummaryEntry.poissonCollector(dateColName,rateColName,epiEstimMode));
//			return out;
//		
//		} else { 
//		
//			RDataframe out2 = summaryResults.stream().collect(
//				GrowthRateEstimationSummaryEntry.poissonSummaryCollector(dateColName,rateColName,epiEstimMode ));
//			return out2;
//		}
//		
//	}
	
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
	
//	/**
//	 * This converts the input timeseries into a EstimateResult time series. This is the native Java methd in the api. 
//	 * @param rtSeries - a GrowthRateTimeseries
//	 * @param samplesPerProfile - the number of times to bootstrap per infectivity profile 
//	 * @return a GrowthRateEstimateResult.
//	 * @throws IncompleteTimeseriesException if the timeseries does not have contiguous dates
//	 */
//	public GrowthRateEstimationResult estimateRtBootstrapped(GrowthRateBootstrappedTimeseries rtSeries, int samplesPerProfile) throws IncompleteTimeseriesException {
//		if (infProf.isEmpty()) throw new RuntimeException("Needs an infectivity profile to be defined before we can estimate Rt");
//		
//		rtSeries.checkComplete();
//		
//		// Make sure infProf has eough entries to allow for bootstrapping 
//		// of timeseries
////		int infProfLen = lcm(samples, infProf.size());
////		int repeats = infProfLen / (infProf.size());
////		if (infProfLen > 10*samples) infProfLen = 10*samples;
//		List<InfectivityProfile> infProfRep = new ArrayList<>();
//		for (int i=0; i<samplesPerProfile; i++) {
//			infProfRep.addAll(infProf); 
//		}
//				
//		GrowthRateEstimationResult allResults = new GrowthRateEstimationResult(this);
//		infProfRep.parallelStream().map(prof -> {
//			
//			GrowthRateEstimationResult singleResult = new GrowthRateEstimationResult(this);
//			GrowthRateTimeseries rtWorking = rtSeries.sample();
//			GrowthRateUtils conv = new GrowthRateUtils(prof);
//			double initialGrowth = conv.RToGrowth(rZero.convert().getMean());
//			rtWorking.inferStart(Optional.ofNullable(atStart ? 0 : null), prof.length(), initialGrowth);
//			rtWorking.estimateForProfile(prof.profile());
//			GrowthRateEstimationResultEntry last = null;
//			Optional<GrowthRateTimeseriesEntry> tsEntry = rtWorking.start();
//			while (tsEntry.isPresent()) {
//				
//				List<DatedGammaEstimate> priors;
//				if (last == null) {
//					priors = defaultPriors(tsEntry.get().dateValue(), prof.getId());
//				} else {
//					priors = priorSelectionStrategy.apply(last);
//				}
//				List<DatedGammaEstimate> results = tsEntry.get().results(priors);
//				
//				GrowthRateEstimationResultEntry current = new GrowthRateEstimationResultEntry(prof.getId(), tsEntry.get(), results, singleResult);
//				singleResult.add(current);
//				last = current;
//				tsEntry = tsEntry.get().next();
//			};
//			
//			return singleResult;
//		
//		}).forEach(allResults::addAll);
//		System.gc();
//		return allResults;
// 	}
//	
//	
//	private List<DatedGammaEstimate> defaultPriors(LocalDate dateValue) {
//		return IntStream.range(0,maxTau).mapToObj(i -> defaultPrior(i, dateValue, profileId)).collect(Collectors.toList());
//	}
//	
//	private List<DatedGammaEstimate> fixedPriors(double mean, double sd, LocalDate dateValue, int profileId) {
//		return IntStream.range(0,maxTau).mapToObj(i -> fixedPrior(mean,sd,i, dateValue, profileId)).collect(Collectors.toList());
//	}
	
//	protected DatedGammaEstimate defaultPrior(int tau, LocalDate dateValue, int profileId) {
//		return rZero.withDate(tau, dateValue, 0D, 0D, profileId);
//	}
//	
//	protected static DatedGammaEstimate fixedPrior(double mean, double sd, int tau, LocalDate dateValue, int profileId) {
//		return new GammaMoments(mean,sd).convert().withDate(tau, dateValue, 0D, 0D, profileId);
//	}


	
}
