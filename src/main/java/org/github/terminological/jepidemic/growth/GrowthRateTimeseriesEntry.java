package org.github.terminological.jepidemic.growth;

import java.time.LocalDate;
import java.util.Arrays;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.github.terminological.jepidemic.Timeseries;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;
import org.github.terminological.jepidemic.distributions.Summary;
import org.github.terminological.jepidemic.distributions.Summary.Impl;
import org.github.terminological.jepidemic.growth.DatedEstimate.LambdaBasedRt;
import org.github.terminological.jepidemic.growth.EstimateMetadata.GrowthMetadata;
import org.github.terminological.jepidemic.growth.EstimateMetadata.RtFromGrowthMetadata;

import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RNumeric;


/**
 * A timeseries entry in a ordered structure with accessors for the previous and next items in the time series. <br>
 * 
 * This is where the component calculations for the GrowthRate method are performed. The calculations for the cori method are not
 * performed until the recache() method is called, which is done by the GrowthRateEstimator::estimateForProfile method when considering a new infectivity profile.
 *  
 *
 */
public class GrowthRateTimeseriesEntry implements TimeseriesEntry.Incidence, Comparable<GrowthRateTimeseriesEntry> {

	private LocalDate date;
	private GrowthRateTimeseries series;
	private Optional<GrowthRateTimeseriesEntry> next;
	private Optional<GrowthRateTimeseriesEntry> prev;
	private boolean inferred;
	private double[] backwardSums;
	private double[] forwardSums;
	
	private double value;
	
	public GrowthRateTimeseriesEntry(TimeseriesEntry.Incidence e, GrowthRateTimeseries series) {
		this(e.date().get(), e.incidence().javaPrimitive(),  series, false);	
	}
	
	public GrowthRateTimeseriesEntry(GrowthRateTimeseriesEntry e, GrowthRateTimeseries series) {
		this(e.dateValue(), e.value(), series, e.isInferred());
		this.inferred = e.inferred;
	}
	
	public GrowthRateTimeseriesEntry(LocalDate date, double value, GrowthRateTimeseries series, boolean inferred) {
		this.value = value;
		this.date = date;
		this.setTimeseries(series);
		this.inferred = inferred;
		this.backwardSums = new double[series.getEstimator().getMaxTau()*2+1];
		Arrays.fill(backwardSums, -1);
		this.forwardSums = new double[series.getEstimator().getMaxTau()*2+1];
		Arrays.fill(forwardSums, -1);
	}

	@Override
	public void setTimeseries(Timeseries<?> ts) {
		this.series = (GrowthRateTimeseries) ts;
	}
	
	public Optional<GrowthRateTimeseriesEntry> next() {
		if (next == null) next = Optional.ofNullable(series.higher(this));
		return next;
	}
	
	public Optional<GrowthRateTimeseriesEntry> prev() {
		if (prev == null) prev = Optional.ofNullable(series.lower(this));
		return prev;
	}
	
	public Optional<GrowthRateTimeseriesEntry> lead(int i) {
		if (i == 0) return Optional.of(this);
		return next().flatMap(n -> n.lead(i-1));
	}
	
	public Optional<GrowthRateTimeseriesEntry> lag(int i) {
		if (i == 0) return Optional.of(this);
		return prev().flatMap(p -> p.lag(i-1));
	}
	
	public double value() { return value;}
	public LocalDate dateValue() { return date;}
	
	@Override
	public RNumeric incidence() {return RConverter.convert(this.value);}

	@Override
	public RDate date() {return RConverter.convert(this.date);}

	@Override
	public int compareTo(GrowthRateTimeseriesEntry o) {
		return this.dateValue().compareTo(o.dateValue());
	}

	public String toString() {
		return dateValue()+"\t"+value();
	}
	
	// I_s is the sum of incidence over the window.
	double sumI_s(int tau) {
		if (this.inferred) return Double.NaN;
		if (tau < 0) return 0;
		double tmp = forwardSumI_s(tau)+backwardSumI_s(tau)-value;
		return tmp;
	}
	
	double backwardSumI_s(int tau) {
		if (tau < 0) return 0;
		if (this.backwardSums[tau] == -1) {
			this.backwardSums[tau]= 
				(this.inferred ? Double.NaN : this.value) + 
				(prev().isPresent() ? prev().get().backwardSumI_s(tau-1) : Double.NaN);
			
		}
		return this.backwardSums[tau];
	}
	
	double forwardSumI_s(int tau) {
		if (tau < 0) return 0;
		if (this.forwardSums[tau] == -1) {
			this.forwardSums[tau] = 
				(this.inferred ? Double.NaN : this.value) +
				(next().isPresent() ? next().get().forwardSumI_s(tau-1) : Double.NaN);
		}
		return this.forwardSums[tau];
	}

	public void setInferred() {
		this.inferred = true;
	}

	public boolean isInferred() {
		return inferred;
	}

	@Override
	public Optional<Timeseries<?>> getTimeseries() {
		return Optional.ofNullable(series);
	}
	
	public GrowthRateTimeseries getActualTimeseries() {
		return series;
	}
	
	
//	Optional<ExtendedGammaDistribution> getRtPrior(RtMetadata meta) {}
//	
//	Collection<DatedEstimate.PoissonRate> getAllLambdaEstimates() {
//		return this.getPrimaryEstimates().getLambda().values();
//	}
	DatedEstimate.PoissonRate getLambdaEstimate(GrowthMetadata meta) {
		return this.getPrimaryEstimates().getLambda().get(meta);
	}
	
	@SuppressWarnings("unchecked")
	Summary.Impl<AbstractRealDistribution> getLambdaSummary() {
		return (Impl<AbstractRealDistribution>) this.getSummaryEstimates().getLambda();
	}
	
//	public Optional<ExtendedGammaDistribution> getLambdaPosterior(GrowthMetadata metadata) {
//		return this.getPrimaryEstimates().lambdaEstimates.optional(metadata).flatMap(l -> l.getDist());
//	}
//	
//	Collection<DatedEstimate.GrowthRate> getAllGrowthRateEstimates() {
//		return this.getPrimaryEstimates().getGrowth().values();
//	}
	DatedEstimate.GrowthRate getGrowthRate(GrowthMetadata meta) {
		return this.getPrimaryEstimates().getGrowth().get(meta);
	}
	
	@SuppressWarnings("unchecked")
	Summary.Impl<AbstractRealDistribution> getGrowthRateSummary() {
		return (Impl<AbstractRealDistribution>) this.getSummaryEstimates().getGrowth();
	}
//	
//	List<DatedEstimate.Incidence> getAllIncidenceEstimates() {}
//	DatedEstimate.Incidence getIncidenceEstimate(GrowthMetadata meta) {}
//	Summary.Impl<? extends AbstractIntegerDistribution> getIncidenceSummary() {}
//	
//	List<DatedEstimate.LambdaBasedRt> getAllRtFromLambdaEstimates() {}
//	DatedEstimate.Incidence getRtFromLambdaEstimate(GrowthMetadata meta) {}
//	Summary.Impl<? extends AbstractRealDistribution> getRtFromLambdaSummary() {}
//	
//	List<DatedEstimate.GammaRt> getAllCoriRtEstimates() {}
//	DatedEstimate.GammaRt getCoriRtEstimate(GrowthMetadata meta) {}
//	Summary.Impl<? extends AbstractRealDistribution> getCoriRtSummary() {}
//	
//	
	
//	public static class Estimation<X extends DatedEstimate<X,GrowthRateTimeseriesEntry,Y>, Y extends EstimateMetadata> {
//		
//		OptionalMap<Y,X> map = new OptionalMap<Y,X>();
//		List<Y> support;
//		X defaultValue;
//		GrowthRateTimeseriesEntry ts;
//		Strategy.PosteriorFiltering<X> posteriorFilter;
//		Strategy.CombiningStrategy<X, ?> combiner;
//		BiFunction<GrowthRateTimeseriesEntry,Y,X> posteriorConstructor;
//		
//		Collection<X> getAllEstimates() {
//			return map.values();
//		}
//		Collection<X> getFilteredEstimates() {
//			return map.values().stream().filter(posteriorFilter).collect(Collectors.toList());
//		}
//		
//		X getEstimate(Y meta) {
//			if (!support.contains(meta)) return defaultValue;
//			if (map.containsKey(meta)) return map.get(meta); 
//			X tmp = posteriorConstructor.apply(ts, meta);
//			map.put(meta, tmp);
//			return tmp;
//		}
//		
//		Summary getSummary() {
//			return combiner.apply(this.getFilteredEstimates()).map(x -> (Summary) x).orElse(Summary.nan());
//		}
//		
//		
//	}
	
	// TODO: SMMR paper
	// TODO: Science paper
	// TODO: jepidemic current CORI implementation testing
	// DONE: deal with memory issues with clean up.
	// TODO: synthetic series with growth rate baseline ? as java class
	// TODO: write up implementation - welch scatterthwaithe approximation test & impact.
	// TODO: Cori method in new framework
	// TODO: Wallinga method of Rt
	// DONE: Write R api 
	
	
	

	
	// Caches
	
	OptionalMap<GrowthMetadata, ExtendedGammaDistribution> priorCache;
	PrimaryEstimates primaryEstimateCache;
	DerivedEstimates derivedEstimateCache;
	CombinedSummary summaryCombinationCache;
	
	
	Optional<ExtendedGammaDistribution> getLambdaPrior(GrowthMetadata meta) {
		if (priorCache == null) {
			// Get previous days results for prior selection.
			priorCache = new OptionalMap<GrowthMetadata, ExtendedGammaDistribution>();
			
			// THIS IS WHERE WE APPLY THE PRIOR SELECTION STRATEGY
			// This will typically be called from teh Dated
			
			Strategy.PriorSelection strat = this.getActualTimeseries()
					.getEstimator()
					.getPriorSelectionStrategy();
			
			this.getActualTimeseries().getGrowthSupport().stream()
					.map(m -> strat.apply(this,m))
					.filter(pr -> pr.getDist().isPresent())
					.forEach(pr -> priorCache.put(pr.getMeta(), pr.getDist().get()));
		}
		return priorCache.optional(meta);
	}
	
	PrimaryEstimates getPrimaryEstimates() {
		
		if (primaryEstimateCache == null) {
		
			// This is where me manage workflow.
			// Find a set of priors for the estimation using a prior selection strategy
			
			primaryEstimateCache = new PrimaryEstimates(this);
			
			for (GrowthMetadata metadata: this.getActualTimeseries().getGrowthSupport()) {
				// Define metadata for posteriors (which may be metadata from priors)
				// Construct the posteriors
				DatedEstimate.PoissonRate posteriorLamda = new DatedEstimate.PoissonRate(this, metadata);
				DatedEstimate.GrowthRate posteriorGrowth = new DatedEstimate.GrowthRate(this, metadata);
				primaryEstimateCache.getLambda().put(metadata, posteriorLamda);
				primaryEstimateCache.getGrowth().put(metadata, posteriorGrowth);
			}
			
			
			
		}
		return primaryEstimateCache;
		
	}	
	
	DerivedEstimates getDerivedEstimates() {
		
		if (derivedEstimateCache == null) {
			
			derivedEstimateCache = new DerivedEstimates(this);
			// This is where me manage workflow.
			// Find a set of priors for the estimation using a prior selection strategy
		
			for (GrowthMetadata metadata: this.getActualTimeseries().getGrowthSupport()) {
				
				DatedEstimate.Incidence posteriorIncidence = new DatedEstimate.Incidence(this, metadata);
				derivedEstimateCache.getIncidence().put(metadata, posteriorIncidence);
				
			}
			
			if (this.getActualTimeseries().getEstimator().doRtEstimation()) {
				for (RtFromGrowthMetadata meta: this.getActualTimeseries().getRtFromGrowthSupport()) {
					DatedEstimate.LambdaBasedRt rtEstimate = new DatedEstimate.LambdaBasedRt(this, meta); 
					derivedEstimateCache.getRt().put(meta, rtEstimate);
				}
			}
			
			
		}
		return derivedEstimateCache;
	}
	
//	DerivedEstimates filteredEstimates = null;
//	PrimaryEstimates filteredLambdaEstimates = null;
	
//	public PrimaryEstimates filteredLambdaEstimates() {
//		if (filteredLambdaEstimates == null) {
//			PrimaryEstimates lambdaEstimates = getPrimaryEstimates();
//			PrimaryEstimates out = new PrimaryEstimates(this);
//			Strategy.PosteriorFiltering<DatedEstimate.PoissonRate> filter = this.getActualTimeseries().getEstimator().getPosteriorLambdaFilteringStrategy();
//			Strategy.PosteriorFiltering<DatedEstimate.GrowthRate> filter3 = this.getActualTimeseries().getEstimator().getPosteriorGrowthFilteringStrategy();
//			lambdaEstimates.lambdaEstimates.values().stream()
//				.filter(filter)
//				.forEach(l -> out.lambdaEstimates.put(l.getMeta(), l));
//			lambdaEstimates.growthEstimates.values().stream()
//				.filter(filter3)
//				.forEach(l -> out.growthEstimates.put(l.getMeta(), l));
//			filteredLambdaEstimates = out;
//		}
//		return filteredLambdaEstimates;
//	}
//	
//	public DerivedEstimates filteredGrowthEstimates() {
//		if (filteredEstimates == null) {
//			DerivedEstimates estimates = getDerivedEstimates();
//			DerivedEstimates out = new DerivedEstimates(this);
//			Strategy.PosteriorFiltering<DatedEstimate.Incidence> filter2 = this.getActualTimeseries().getEstimator().getPosteriorIncidenceFilteringStrategy();
//			Strategy.PosteriorFiltering<LambdaBasedRt> filter4 = this.getActualTimeseries().getEstimator().getPosteriorRtFilteringStrategy();
//			estimates.incidenceEstimates.values().stream()
//				.filter(filter2)
//				.forEach(l -> out.incidenceEstimates.put(l.getMeta(), l));
//			if (this.getActualTimeseries().getEstimator().doRtEstimation()) {
//				estimates.rtEstimates.values().stream()
//					.filter(filter4)
//					.forEach(l -> out.rtEstimates.put(l.getMeta(), l));	
//			}
//			filteredEstimates = out;
//		}
//		return filteredEstimates;
//	}
	
	
	public CombinedSummary getSummaryEstimates() {
		
		if (this.summaryCombinationCache == null) {
			CombinedSummary out = new CombinedSummary(this);
			
			// THIS APPLIES POSTERIOR FILTERING AND COMBINATION STRATEGIES
			// TODO: this has been cut down to just be a filter.
			// unfortunately this measn we can't select outside of estimates for the 
			// day and can;t 
			
			GrowthRateEstimator es = this.getActualTimeseries().getEstimator();
			Strategy.PosteriorFiltering filt = es.getPosteriorFilteringStrategy();
			
			
			PrimaryEstimates lambdaEstimates = getPrimaryEstimates();
			out.setGrowth( 
					es.getCombiningGrowthRateStrategy().apply(
						lambdaEstimates.getGrowth().values().stream().filter(filt)
							.collect(Collectors.toList())));
			
			out.setLambda(
					es.getCombiningLambdaStrategy().apply(
							lambdaEstimates.lambdaEstimates.values().stream().filter(filt)
								.collect(Collectors.toList())));
			
			DerivedEstimates estimates = getDerivedEstimates();
			out.setIncidence( 
					es.getCombiningIncidenceStrategy().apply(
						estimates.incidenceEstimates.values().stream().filter(filt)
							.collect(Collectors.toList())));
			
			if (this.getActualTimeseries().getEstimator().doRtEstimation()) {
				out.setRt( 
					es.getCombiningRtStrategy().apply(
							estimates.rtEstimates.values().stream().filter(filt)
								.collect(Collectors.toList())));
			}
			summaryCombinationCache = out;
		}
		return summaryCombinationCache;
	}
	
	public static class PrimaryEstimates {
		
		PrimaryEstimates(GrowthRateTimeseriesEntry ts) {
			this.ts = ts;
		}
				
		GrowthRateTimeseriesEntry ts;
		OptionalMap<GrowthMetadata,DatedEstimate.PoissonRate> lambdaEstimates = new OptionalMap<>();
		OptionalMap<GrowthMetadata,DatedEstimate.GrowthRate> growthEstimates = new OptionalMap<>();
		
		GrowthRateTimeseriesEntry getTsEntry() {return ts;}
		OptionalMap<GrowthMetadata,DatedEstimate.PoissonRate> getLambda() {return lambdaEstimates;}
		OptionalMap<GrowthMetadata,DatedEstimate.GrowthRate> getGrowth() {return growthEstimates;}

		
	}
	
	public static class DerivedEstimates {
		
		DerivedEstimates(GrowthRateTimeseriesEntry ts) {
			this.ts = ts;
		}
		
		GrowthRateTimeseriesEntry ts;
		OptionalMap<GrowthMetadata,DatedEstimate.Incidence> incidenceEstimates = new OptionalMap<>();
		OptionalMap<RtFromGrowthMetadata,DatedEstimate.LambdaBasedRt> rtEstimates = new OptionalMap<>();
		
		GrowthRateTimeseriesEntry getTsEntry() {return ts;}
		OptionalMap<GrowthMetadata,DatedEstimate.Incidence> getIncidence() {return incidenceEstimates;}
		OptionalMap<RtFromGrowthMetadata, LambdaBasedRt> getRt() {return rtEstimates;}
		
	}
	
	public static class CombinedSummary {
		
		CombinedSummary(GrowthRateTimeseriesEntry ts) {
			this.ts = ts;
		}
		
		GrowthRateTimeseriesEntry ts;
		//GrowthMetadata metadata;
		Optional<Summary> lambdaEstimate = Optional.empty();
		Optional<Summary> growthEstimate = Optional.empty();
		Optional<Summary> incidenceEstimate = Optional.empty();
		Optional<Summary> rtEstimate = Optional.empty();
		
		public void setRt(Optional<AbstractRealDistribution> rtEstimate) {this.rtEstimate = rtEstimate.map(Summary::of);}
		public void setGrowth(Optional<AbstractRealDistribution> growthEstimate) {this.growthEstimate = growthEstimate.map(Summary::of);}
		public void setIncidence(Optional<AbstractIntegerDistribution> incidenceEstimate) {this.incidenceEstimate = incidenceEstimate.map(Summary::of);}
		public void setLambda(Optional<AbstractRealDistribution> lambdaEstimate) {this.lambdaEstimate = lambdaEstimate.map(Summary::of);}
		
		public Summary getRt() {return rtEstimate.orElse(Summary.nan());}
		public Summary getGrowth() {return growthEstimate.orElse(Summary.nan());}
		public Summary getIncidence() {return incidenceEstimate.orElse(Summary.nan());}
		public Summary getLambda() {return lambdaEstimate.orElse(Summary.nan());}
		
		public void collapseAll(double[] quantiles) {
			lambdaEstimate = lambdaEstimate.map(c -> c.collapse(quantiles)); 
			growthEstimate = growthEstimate.map(c -> c.collapse(quantiles));
			incidenceEstimate = incidenceEstimate.map(c -> c.collapse(quantiles));
			rtEstimate = rtEstimate.map(c -> c.collapse(quantiles));
			
		}
		
	}

	public Stream<GrowthRateTimeseriesEntry> streamForward(int maxDaysAhead) {
		if (maxDaysAhead < 0) return Stream.empty();
		return Stream.concat(
				Stream.of(this), 
				this.next().map(o -> Stream.of(o)).orElse(Stream.empty()).flatMap(n -> n.streamForward(maxDaysAhead-1)));
		
	}
	
	public Stream<GrowthRateTimeseriesEntry> streamBackward(int maxDaysBefore) {
		if (maxDaysBefore < 0) return Stream.empty();
		return Stream.concat(
				Stream.of(this), 
				this.prev().map(o -> Stream.of(o)).orElse(Stream.empty()).flatMap(n -> n.streamForward(maxDaysBefore-1)));
		
	}

	public void cleanUp() {
		int maxTau = this.getActualTimeseries().getEstimator().getMaxTau();
		int postLength = this.getActualTimeseries().getEstimator().getMaxInfProf().orElse(0);
		int retain = Math.max(maxTau*2+1, postLength+5);
		if (this.next().isPresent()) {
			
			this.lag(retain).ifPresent(c->c.freeSpace());
		} else {
			this.streamBackward(retain).forEach(c-> c.freeSpace());		}
		
	}
	
	public void freeSpace() {
		double[] quantiles = this.getActualTimeseries().getEstimator().getQuantiles();
		this.getSummaryEstimates().collapseAll(quantiles);
		this.primaryEstimateCache = null;
		this.derivedEstimateCache = null;
		this.priorCache = null;
	}
	
//	// calculate posteriors basd on priors
//	public DatedEstimate.PoissonRate lambdaPosterior(DatedEstimate.PoissonRate lambdaPrior, int tau) {
//		return new DatedEstimate.PoissonRate(tau, lambdaPrior, this);
//	}
//	
//	// This constructs all the posteriors for a list of priors 
//	public List<DatedEstimate.PoissonRate> lambdaPosteriors() {
//		List<DatedEstimate.PoissonRate> lambdaPriors = this.series.getEstimator().priorSelectionStrategy.apply(this);
//		return lambdaPriors.stream().map(
//				p -> new DatedEstimate.PoissonRate(p.tau, p, this)).collect(Collectors.toList());
//	}
//	
////	// This constructs all the posteriors for a list of priors 
////	public List<DatedEstimate.PoissonRate> lambdaPosteriors(DatedEstimate.PoissonRate lambdaPrior) {
////		return IntStream.range(1,series.getMaxTau()).mapToObj(
////				i -> new DatedEstimate.PoissonRate(i, lambdaPrior, this)).collect(Collectors.toList());
////	}
//	
//	public DatedEstimate.GrowthRate growthPosterior(DatedEstimate.PoissonRate lambdaPrior) {
//		return new DatedEstimate.GrowthRate(lambdaPrior, this);
//	}
//	
//	public DatedEstimate.Incidence incidencePosterior(DatedEstimate.PoissonRate lambdaPrior) {
//		return new DatedEstimate.Incidence(lambdaPrior, this);
//	}
	
	
	
	
	
}
