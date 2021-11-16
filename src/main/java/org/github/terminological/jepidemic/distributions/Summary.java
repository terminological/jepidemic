package org.github.terminological.jepidemic.distributions;

import java.util.Arrays;
import java.util.Map;
import java.util.Optional;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.distribution.RealDistribution;

public interface Summary {

	public double getMean();
	public double getSD();
	public double quantile(double q);
	
	default public Summary collapse(double[] quantiles) {
		return Summary.collapse(this, quantiles);
	}
	
	default public Map<Double,Double> quantiles(double[] quantiles) {
		SortedMap<Double,Double> out = new TreeMap<Double,Double>();
		Arrays.stream(quantiles).boxed().forEach(p -> out.put(p, (double) this.quantile(p)));
		return out;
	}
	
	public static abstract class Impl<X> implements Summary {
		public abstract Optional<X> getDistribution();
		public String toString() {
				return String.format("%1.3f +- %1.3f [%1.3f - %1.3f]", this.getMean(), this.getSD(), this.quantile(0.025), this.quantile(0.975));
		}
	}
	
//	public default CoriEstimationSummaryEntry withDate(LocalDate date, int window, double incidence) {
//		return new CoriEstimationSummaryEntry(this, date, window, incidence);
//	}
	
	public static <X extends RealDistribution> Summary.Impl<X> of(X dist) {
		return new Summary.Impl<X>() {
			@Override
			public double getMean() {
				return dist.getNumericalMean();
			}
			@Override
			public double getSD() {
				return Math.sqrt(dist.getNumericalVariance());
			}
			@Override
			public double quantile(double q) {
				try {
					return dist.inverseCumulativeProbability(q);
				} catch (Exception e) {
					return Double.NaN;
				}
			}
			@Override
			public Optional<X> getDistribution() {
				return Optional.of(dist);
			}
		};
	}
	
	public static <X> Summary.Impl<X> nan() {
		return new Summary.Impl<X>() {

			@Override
			public double getMean() {
				return Double.NaN;
			}

			@Override
			public double getSD() {
				return Double.NaN;
			}

			@Override
			public double quantile(double q) {
				return Double.NaN;
			}

			@Override
			public Optional<X> getDistribution() {
				return Optional.empty();
			}
			
		};
	}
	
	public static <X extends IntegerDistribution> Summary.Impl<X> of(X dist) {
		return new Summary.Impl<X>() {
			@Override
			public double getMean() {
				return dist.getNumericalMean();
			}
			@Override
			public double getSD() {
				return Math.sqrt(dist.getNumericalVariance());
			}
			@Override
			public double quantile(double q) {
				return dist.inverseCumulativeProbability(q);
			}
			@Override
			public Optional<X> getDistribution() {
				return Optional.of(dist);
			}
		};
	}
	
	public static Summary collapse(Summary summary, double[] quantiles) {
		return new Summary() {

			double mean = summary.getMean();
			double sd = summary.getSD();
			Map<Double,Double> quantileMap = summary.quantiles(quantiles);
			
			@Override
			public double getMean() {
				return mean;
			}

			@Override
			public double getSD() {
				return sd;
			}

			@Override
			public double quantile(double q) {
				return quantileMap.getOrDefault(q, Double.NaN);
			}
			
		};
	}
	
	
}
