package org.github.terminological.jepidemic.growth;

import java.util.Collection;
import java.util.Map;
import java.util.Optional;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.github.terminological.jepidemic.distributions.EmpiricalDistribution;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;
import org.github.terminological.jepidemic.distributions.GrowthRateDistribution;
import org.github.terminological.jepidemic.distributions.Summary;
import org.github.terminological.jepidemic.distributions.WeightedMixtureIntegerDistribution;
import org.github.terminological.jepidemic.distributions.WeightedMixtureRealDistribution;
import org.github.terminological.jepidemic.growth.DatedEstimate.Continuous;
import org.github.terminological.jepidemic.growth.DatedEstimate.Discrete;
import org.github.terminological.jepidemic.growth.DatedEstimate.PoissonRate;
import org.github.terminological.jepidemic.growth.EstimateMetadata.GrowthMetadata;


public class Strategy {

	public static interface PriorSelection extends BiFunction<GrowthRateTimeseriesEntry, GrowthMetadata, DatedEstimate.PoissonRate> {
		
//		//Make a list of priors associated with a time series entry
//		public static PriorSelection fixedPrior(int minTau, int maxTau, ExtendedGammaDistribution prior) {
//			return ts -> IntStream
//					.range(minTau,maxTau)
//					.mapToObj(tau -> new PoissonRate(prior,tau,ts))
//					.collect(Collectors.toList());
//		}
		
		
		//Make a list of priors associated with a time series entry
		public static PriorSelection mechanisticPrior(double factor) {
			return (ts,meta) -> {
				ExtendedGammaDistribution prior = ts.prev().map(prev -> {
					Optional<ExtendedGammaDistribution> lambda = prev.getLambdaEstimate(meta).getDist();
					Optional<GrowthRateDistribution> growth = prev.getGrowthRate(meta).getDist();
				 	if (!lambda.isPresent() || !growth.isPresent()) return meta.defaultLambdaPrior();
					EmpiricalDistribution ePrior = EmpiricalDistribution
						.bySampling(
							lambda.get(), 
							growth.get(), 
							(x,y) -> x*Math.exp(y), 50);
					Summary moments = Summary.of(ePrior);
					ExtendedGammaDistribution gPrior = ExtendedGammaDistribution.fromMoments(moments.getMean(), moments.getSD()*factor);
					return gPrior;
				}).orElse(meta.defaultLambdaPrior());
				return new DatedEstimate.PoissonRate(prior,ts,meta);
			};
		}
		
		public static PriorSelection combinedMechanisticPrior() {
			return (ts,meta) -> {
				Optional<Double> gr2 = ts.lag(2).flatMap(p -> Optional.ofNullable(p.getLambdaSummary().getMean()));
				Optional<Double> gr1 = ts.lag(1).flatMap(p -> Optional.ofNullable(p.getLambdaSummary().getMean()));
				Double factor = gr1.flatMap(g1 -> gr2.map(
						g2 -> Math.exp(Math.abs(g1-g2)/g2))
						).orElse(2.0);
				return combinedMechanisticPrior(factor).apply(ts,meta);
			};
		}
		
		//Make a list of priors associated with a time series entry
		public static PriorSelection combinedMechanisticPrior(double factor) {
			if(factor < 1) throw new RuntimeException("Factor must be greater than 1");
			return ((ts,meta) -> { 
				Optional<AbstractRealDistribution> mixLambda = ts.prev().flatMap(p -> p.getLambdaSummary().getDistribution());
				Optional<AbstractRealDistribution> mixGrowth = ts.prev().flatMap(p -> p.getGrowthRateSummary().getDistribution());
				EmpiricalDistribution ePrior = EmpiricalDistribution
						.bySampling(
							mixLambda.orElse(meta.defaultLambdaPrior()), 
							mixGrowth.orElse(new NormalDistribution(0,0.2)), 
							(x,y) -> x*Math.exp(y), 50);
				Summary moments = Summary.of(ePrior);
				ExtendedGammaDistribution prior = ExtendedGammaDistribution.fromMoments(moments.getMean(), moments.getSD()*factor);
				
				return new PoissonRate(prior, ts, meta);
				
			});
		}
		
		public static PriorSelection posteriorAsPrior(double factor) {
			if(factor < 1) throw new RuntimeException("Factor must be greater than 1");
			return (ts,meta) -> { 
				Optional<AbstractRealDistribution> ePost1 = ts.prev().flatMap(t -> t.getLambdaSummary().getDistribution());
				Optional<ExtendedGammaDistribution> prior = ePost1
						.map(s-> ExtendedGammaDistribution.fromMoments(s.getNumericalMean(), Math.sqrt(s.getNumericalVariance())*factor));
				
				return new PoissonRate(prior.orElse(meta.defaultLambdaPrior()), ts, meta);
				
			};
		}
		
		public static PriorSelection posteriorAsPrior() {
			return (ts,meta) -> { 
				Optional<Summary> ePost1 = ts.lag(1).map(t -> t.getLambdaSummary());
				Optional<Double> ePost2 = ts.lag(2).flatMap(t -> Optional.ofNullable(t.getLambdaSummary().getMean()));
				double factor =  ePost1.map(s -> s.getMean())
						.flatMap(
							p1 -> ePost2.map(
									p2 -> Math.exp(Math.abs(p1-p2)/p1)
							)
						).orElse(2.0);
				
				Optional<ExtendedGammaDistribution> prior = ePost1
						.map(s-> ExtendedGammaDistribution.fromMoments(s.getMean(), s.getSD()*factor));
				
				return new PoissonRate(prior.orElse(meta.defaultLambdaPrior()), ts, meta);
				
			};
		}
		
	}
	
	//TODO: the use of this this is not implemented as it is for the Cori method when
	// it has been recreated in this framework.
	// it will be implemented in the GrowthRateTimeseriesEntry.getSummaryEstimates function
	// but will also need some additional configuration.
//	public static interface PosteriorSelection<X extends DatedEstimate<?,? extends GrowthRateTimeseriesEntry,?>> extends Function<GrowthRateTimeseriesEntry, Collection<? extends X>> {
//		
//		public static PosteriorSelection<LambdaBasedRt> rtWithinWindowSelector(int maxDaysAhead) {
//			return (ts -> 
//				ts.streamForward(maxDaysAhead)
//					.flatMap(tsN -> tsN.getDerivedEstimates().getRt().values().stream())
//					.filter(rt -> !rt.getStartDate().isAfter(ts.getDate()) && rt.getEndDate().isBefore(ts.getDate()))
//					.collect(Collectors.toList())
//			);
//		}
//		
//		public static PosteriorSelection<LambdaBasedRt> rtFromTodaySelector() {
//			return (ts -> ts.getDerivedEstimates().getRt().values().stream()
//					.collect(Collectors.toList()));
//		}
//		
//	}
	
	public static interface PosteriorFiltering extends Predicate<DatedEstimate<?,? extends GrowthRateTimeseriesEntry,?>> {
		
		// allows an estimate through if its window > minTau and the timeseries entry has a certain minimum amount of incidence associated with it
		public static PosteriorFiltering 
			minIncidence(double totalIncidence) {
				return pr -> {
					int tau = pr.getMeta().getTau();
					double incidenceSum = pr.getTsEntry().sumI_s(tau); 
					return incidenceSum > totalIncidence;
				};
			}
		
		// allows an estimate through if its window > minTau and the timeseries entry have a certain minimum amount of incidence associated with it
		public static PosteriorFiltering 
				specificWindow(int tau) {
					return pr -> tau == pr.getMeta().getTau();
				}
		
		// allows an estimate through if its window > minTau and the timeseries entry have a certain minimum amount of incidence associated with it
		public static PosteriorFiltering 
				allowAll() {
					return pr -> true;
				}
		
	}
	
	public static interface CombiningStrategy<X,Y> extends Function<Collection<? extends X>,Optional<Y>> {
		
		// TODO: a resampling strategy?
		
		public static interface ForContinuous<X extends AbstractRealDistribution> extends CombiningStrategy<Continuous<X,?,?>,AbstractRealDistribution> {
			// can be combined with a posterior selection
			public static <A extends AbstractRealDistribution> ForContinuous<A> collectFirst() {
				return list -> {
					if (list.size() == 0) return Optional.empty();
					// if (list.size() > 1) throw new RuntimeException("Only one result expected");
					return list.iterator().next().getDist().map(a -> (AbstractRealDistribution) a);
				};
			}
			
			// can be combined with a posterior selection
			public static <A extends AbstractRealDistribution> ForContinuous<A> weightedMixture(Map<Integer,Double> weightsForTau, double hardLowerLimit, double hardUpperLimit) {
				return list -> {
					
					WeightedMixtureRealDistribution out = WeightedMixtureRealDistribution.empty(hardLowerLimit, hardUpperLimit);
					for (Continuous<A,?,?> item : list) {
						int tau = item.getMeta().getTau();
						double weight = weightsForTau.getOrDefault(tau,0D);
						if( weight>0 && item.getDist().isPresent()) {
							out.add(item.getDist().get(), weight);
						}
					}
					return out.ifNotEmpty();
				};
			}
		}
		
		public static interface ForDiscrete<X extends AbstractIntegerDistribution> extends CombiningStrategy<Discrete<X,?,?>,AbstractIntegerDistribution> {
			// can be combined with a posterior selection
			public static <A extends AbstractIntegerDistribution> ForDiscrete<A> collectFirst() {
				return list -> {
					if (list.size() == 0) return Optional.empty();
					// if (list.size() > 1) throw new RuntimeException("Only one result expected");
					return list.iterator().next().getDist().map(a -> (AbstractIntegerDistribution) a);
				};
			}
			
			// can be combined with a posterior selection
			public static <A extends AbstractIntegerDistribution> ForDiscrete<A> weightedMixture(Map<Integer,Double> weightsForTau) {
				return list -> {
					
					WeightedMixtureIntegerDistribution out = WeightedMixtureIntegerDistribution.empty();
					for (Discrete<A,?,?> item : list) {
						int tau = item.getMeta().getTau();
						double weight = weightsForTau.getOrDefault(tau,0D);
						if( weight>0 && item.getDist().isPresent()) {
							out.add(item.getDist().get(), weight);
						}
					}
					return out.ifNotEmpty();
				};
			}
			
			
		}
	}
}
