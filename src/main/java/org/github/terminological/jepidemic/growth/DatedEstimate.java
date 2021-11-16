package org.github.terminological.jepidemic.growth;

import java.time.LocalDate;
import java.util.Optional;
import java.util.stream.DoubleStream;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.github.terminological.jepidemic.InfectivityProfile;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.distributions.BetaPrimeDistribution;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;
import org.github.terminological.jepidemic.distributions.GrowthRateDistribution;
import org.github.terminological.jepidemic.distributions.NegativeBinomialDistribution;
import org.github.terminological.jepidemic.distributions.Summary;
import org.github.terminological.jepidemic.growth.EstimateMetadata.GrowthMetadata;
import org.github.terminological.jepidemic.growth.EstimateMetadata.RtFromGrowthMetadata;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class DatedEstimate<X, Y extends TimeseriesEntry, Z extends EstimateMetadata> {

	static Logger log = LoggerFactory.getLogger(DatedEstimate.class);
	
	Optional<X> distribution = Optional.empty();
	Y timeseriesEntry;
	Z metadata;
	
	public DatedEstimate(X distribution, Y timeseriesEntry, Z metadata) {
		this.distribution = Optional.ofNullable(distribution);
		this.timeseriesEntry = timeseriesEntry;
		this.metadata = metadata;
	}
	
	public DatedEstimate(Y timeseriesEntry, Z metadata) {
		this.timeseriesEntry = timeseriesEntry;
		this.metadata = metadata;
		// this.doEstimate(timeseriesEntry, metadata);
	}
	
	// public abstract void doEstimate(Y timeseriesEntry, Z metadata);
	
	public void setDistribution(X distribution) {
		this.distribution = Optional.ofNullable(distribution);
	}
	
	Optional<X> getDist() {return distribution;}
	Y getTsEntry() {return timeseriesEntry;}
	Z getMeta() {return metadata;}
	LocalDate getDate() {return timeseriesEntry.getDate();}
	LocalDate getStartDate() {return getDate().minusDays(metadata.getTau());}
	LocalDate getEndDate() {return metadata.isTwoSided() ? getDate() : getDate().plusDays(metadata.getTau());}
	int getWindowSize() {return metadata.isTwoSided() ? 2*metadata.getTau()+1 : metadata.getTau()+1;}
	abstract Summary getSummary();
	
	
	public static abstract class Continuous<A extends AbstractRealDistribution,B extends TimeseriesEntry, C extends EstimateMetadata> extends DatedEstimate<A,B,C> {

		public Continuous(A distribution, B timeseriesEntry, C metadata) {
			super(distribution, timeseriesEntry, metadata);
		}
		
		public Continuous(B timeseriesEntry, C metadata) {
			super(timeseriesEntry, metadata);
		}
		
		public Summary getSummary() {return getDist().map(Summary::of).orElse(Summary.nan());}
		
		
	}
	
	public static abstract class Discrete<A extends AbstractIntegerDistribution,B extends TimeseriesEntry, C extends EstimateMetadata> extends DatedEstimate<A,B,C>  {

		public Discrete(A distribution, B timeseriesEntry, C metadata) {
			super(distribution, timeseriesEntry, metadata);
		}
		
		public Discrete(B timeseriesEntry, C metadata) {
			super(timeseriesEntry, metadata);
		}
		
		public Summary getSummary() {return getDist().map(Summary::of).orElse(Summary.nan());}
		
		
	}
	
	
// 	
//	public static class GammaRt extends Continuous<ExtendedGammaDistribution, CoriTimeseriesEntry, RtMetadata> {
//
//		public GammaRt(GammaRt prior, CoriTimeseriesEntry forEntry, RtMetadata metadata) {
//			// start with distribution unset
//			super(null, forEntry, metadata);
//			int tau = metadata.getTau();
//			// TODO: infectivityProfile must be defined elsewhere and is not used here.
//			// this is because the profile is fixed for a CoriTimeseries as a lot of the calculation can be done up front.
//			// this needs to be rethought.
//			ExtendedGammaDistribution p = prior.getDist().orElse(metadata.defaultR0Prior());
//			double posteriorShape = p.getShape()+forEntry.sumI_s(tau);
//			double posteriorRate = p.getScale()+forEntry.sumLambda_s(tau);
//			
//			// assumes NaN in parameters return null in static methods
//			this.setDistribution(ExtendedGammaDistribution.fromShapeAndRate(posteriorShape, posteriorRate));
//			
//			
//		}
//		
//	}
	
//	public static interface DatedMetadataProvider<A extends EstimateMetadata, B extends TimeseriesEntry> {
//		A getMeta();
//		B getTsEntry();
//	}
	
	public static class PoissonRate extends Continuous<ExtendedGammaDistribution, GrowthRateTimeseriesEntry, GrowthMetadata> {

		public PoissonRate(GrowthRateTimeseriesEntry forEntry, GrowthMetadata metadata) {
			// start with distribution unset
			super(forEntry, metadata);
			
			int tau = metadata.getTau();
			
			ExtendedGammaDistribution p = forEntry
					.getLambdaPrior(metadata)
					.orElse(metadata.defaultLambdaPrior());
		
			double posteriorShape = p.getShape()+forEntry.sumI_s(tau);
			double posteriorRate = p.getRate()+2*tau+1;
			
			// assumes NaN in parameters return null in static methods
			this.setDistribution(ExtendedGammaDistribution.fromShapeAndRate(posteriorShape, posteriorRate));
			
		}
		
		public PoissonRate(ExtendedGammaDistribution givenDistribution, GrowthRateTimeseriesEntry forEntry, GrowthMetadata metadata) {
			super(givenDistribution, forEntry, metadata);
		}
		
	}
	
	public static class GrowthRate extends Continuous<GrowthRateDistribution, GrowthRateTimeseriesEntry, GrowthMetadata> {

		Optional<PoissonRate> prior;
		
		public GrowthRate(GrowthRateTimeseriesEntry forEntry, GrowthMetadata metadata) {
			// start with distribution unset
			super(forEntry, metadata);
			
			int tau = metadata.getTau();
			ExtendedGammaDistribution p = forEntry
					.getLambdaPrior(metadata)
					.orElse(metadata.defaultLambdaPrior());
			
			double alphaPlus = forEntry.forwardSumI_s(2*tau)+p.getShape();
			double alphaMinus = forEntry.backwardSumI_s(2*tau)+p.getShape();
			
			// assumes NaN in parameters return null in static methods
			this.setDistribution(GrowthRateDistribution.withForwardEstimateAndBackwardEstimate(alphaPlus, alphaMinus, tau));
			
		}
		
		public GrowthRate(GrowthRateDistribution givenDistribution, GrowthRateTimeseriesEntry forEntry, GrowthMetadata metadata) {
			super(givenDistribution, forEntry, metadata);
			this.prior = Optional.empty();
		}
		
	}
	
	public static class Incidence extends Discrete<NegativeBinomialDistribution, GrowthRateTimeseriesEntry, GrowthMetadata> {

		public Incidence(GrowthRateTimeseriesEntry forEntry, GrowthMetadata metadata) {
			// start with distribution unset
			super(forEntry, metadata);
			
			Optional<ExtendedGammaDistribution> posterior = forEntry.getLambdaEstimate(metadata).getDist();
			
			double alphaPrime = posterior.map(d -> d.getShape()).orElse(Double.NaN);
			double betaPrime = posterior.map(d -> d.getRate()).orElse(Double.NaN);
			
			// assumes NaN in parameters return null in static methods
			this.setDistribution(NegativeBinomialDistribution.fromSizeAndProbability(
				alphaPrime,
				1 / (betaPrime+1)
			));
			
		}
		
	}
	
	public static class LambdaBasedRt extends Continuous<BetaPrimeDistribution, GrowthRateTimeseriesEntry, RtFromGrowthMetadata> {

		public LambdaBasedRt(GrowthRateTimeseriesEntry forEntry, RtFromGrowthMetadata metadata) {
			// start with distribution unset
			super(forEntry, metadata);
			
			// TODO: growth rate based Rt estimate using Wallinga (2007) and sampling for growth rate uncertainty.
			// Also possible to make method an item in the metadata although distributional form of result here is stated
			// Probably can make it an abstract real distribution
			// Third method exists which is the Poisson ratio, again would have to be sampled
			
			InfectivityProfile ip = metadata.getInfectivityProfile();
			
			// This relies on inheritance of RtFromGrowthMetadata to generate key to find relevant lambda estimate
			GrowthMetadata metaForTau = new GrowthMetadata(metadata); 
			Optional<ExtendedGammaDistribution> t0posterior = forEntry.getLambdaEstimate(metaForTau).getDist();
			
			if (!t0posterior.isPresent()) return;
			double alphaPrime = t0posterior.get().getShape();
			double betaPrime = t0posterior.get().getRate();
			
			double sumAlphaPrimeOmegaOnBetaPrime = 0;
			double sumAlphaPrimeOmega2OnBetaPrime2 = 0;
			
			// The 0 index case is included but by infecivity profile definition 
			// omega will be zero.
			double[] ipProfile = ip.profile();
			if (ipProfile[0] != 0) log.warn("Infectivity profile was not zero at time zero. Ignoring first value.");
			
			double total = DoubleStream.of(ipProfile).sum();
			
			for(int i = 1; i< ipProfile.length; i++) {
				double omega = ipProfile[i]/total;
				
				// Find the posterior lamda estimate for the i'th previous day
				Optional<ExtendedGammaDistribution> tNposterior = forEntry.lag(i)
						.flatMap(ts -> ts.getLambdaEstimate(metaForTau).getDist());
						
				// If any of the posterior values for previous days are empty we can't calculate R_t
				// This exits early leaving distribution unset - which is handled as no estimate possible.
				if (!tNposterior.isPresent()) {
					this.setDistribution(null);
					return;		
				}
				
				double tN_alphaPrime = tNposterior.get().getShape();
				double tN_betaPrime = tNposterior.get().getRate();
				
				sumAlphaPrimeOmegaOnBetaPrime += tN_alphaPrime * omega / tN_betaPrime;
				sumAlphaPrimeOmega2OnBetaPrime2 += tN_alphaPrime * Math.pow(omega,2) / Math.pow(tN_betaPrime,2);
			}
			
			double alphaPrimePrime = Math.pow(sumAlphaPrimeOmegaOnBetaPrime, 2) / sumAlphaPrimeOmega2OnBetaPrime2;
			double betaPrimePrime = alphaPrimePrime / sumAlphaPrimeOmegaOnBetaPrime;
			
			BetaPrimeDistribution est = BetaPrimeDistribution
					.extraDistrParameters(alphaPrime, alphaPrimePrime, betaPrimePrime/betaPrime);
			
			this.setDistribution(est);
			
			
			return;
		}
		
		
		
	}
}
