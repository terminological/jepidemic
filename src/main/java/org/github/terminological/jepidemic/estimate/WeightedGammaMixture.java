package org.github.terminological.jepidemic.estimate;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.util.Pair;
import org.github.terminological.jepidemic.distributions.Summary;

public class WeightedGammaMixture extends AbstractRealDistribution implements Summary {

	public static final double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;
    
	private List<Pair<GammaParameters, Double>> weightedMap = new ArrayList<>();
	private double totalWeight = 0;
	private double mean = -1;
	private double var = -1;
	private Map<Double,Double> quantiles = new HashMap<>();
	private double acc_max = Double.NEGATIVE_INFINITY;
//	private double acc_min = Double.POSITIVE_INFINITY;
//	private UnivariateIntegrator intr;
	
	/**
     * Creates a Mixture of Gamma distributions.
     *
     * @param <X> - Dated or undated gamma parameters
     * @param gammas - a list of the gamma distributions of the mixture
     */
    public <X extends GammaParameters> WeightedGammaMixture(List<X> gammas) {
    	this(gammas, g -> 1.0);
    }
    
    public <X extends GammaParameters> WeightedGammaMixture(List<X> gammas, Function<X,Double> mapper) {
    	super(null);
    	gammas.stream().filter(g -> g.isDefined()).forEach(g -> {
    		double wt = mapper.apply(g);
    		weightedMap.add(new Pair<>(g, wt));
    		totalWeight += wt;
    		double ul = g.convert().getMean()+10*g.convert().getSD();
    		if (acc_max<ul) acc_max=ul;
//    		double ll = g.convert().getMean()-10*g.convert().getSD();
//    		if (acc_min>ll) acc_min=ll;
    	});
 //   	this.intr = new TrapezoidIntegrator(0.001,0.001,TrapezoidIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,TrapezoidIntegrator.TRAPEZOID_MAX_ITERATIONS_COUNT); 
    }

    
	@Override
	public double density(double x) {
		if (totalWeight == 0) return Double.NaN;
		if (x <= 0.0) return 0;
		if (x == Double.POSITIVE_INFINITY) return 0;
		double out = weightedMap.stream().mapToDouble(kv -> kv.getKey().dist().density(x)*kv.getValue()).sum()/totalWeight;
		return out;
	}

	@Override
	public double cumulativeProbability(double x) {
		if (totalWeight == 0) return Double.NaN;
		if (x <= 0.0) return 0;
		if (x == Double.POSITIVE_INFINITY) return 1;
		double out = weightedMap.stream().mapToDouble(kv -> kv.getKey().dist().cumulativeProbability(x)*kv.getValue()).sum()/totalWeight;
		return out;
	}

	@Override
	// https://en.wikipedia.org/wiki/Mixture_distribution
	public double getNumericalMean() {
		if (totalWeight == 0) return Double.NaN;
		if (mean < 0) {
			Double sum = this.weightedMap.stream()
					.map(x -> x.getSecond()/totalWeight*x.getFirst().dist().getNumericalMean())
					.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
					.reduce(Double::sum).orElse(Double.NaN);
			mean = sum;
//			UnivariateFunction f = new UnivariateFunction() {
//				@Override
//				public double value(double x) {
//					return x * WeightedGammaMixture.this.density(x);
//				}};
//			mean = this.intr.integrate(100000, f, this.getSupportLowerBound(), this.acc_max);
		}
		return mean;
	}

	@Override
	// https://en.wikipedia.org/wiki/Mixture_distribution
	public double getNumericalVariance() {
		if (totalWeight == 0) return Double.NaN;
		if (var < 0) {
			Double sum = this.weightedMap.stream()
					.map(x -> x.getSecond()/totalWeight*(x.getFirst().dist().getNumericalVariance()+Math.pow(x.getFirst().dist().getNumericalMean(),2)))
					.reduce(Double::sum).orElse(Double.NaN);
			var = sum-Math.pow(getNumericalMean(),2);
//			UnivariateFunction f = new UnivariateFunction() {
//				@Override
//				public double value(double x) {
//					return Math.pow(x,2) * WeightedGammaMixture.this.density(x);
//				}};
//			var = this.intr.integrate(100000, f, this.getSupportLowerBound(), this.acc_max) - 
//					Math.pow(this.getNumericalMean(),2);
//			// https://en.wikipedia.org/wiki/Variance#Absolutely_continuous_random_variable
		}
		return var;
	}

	@Override
	public double getSupportLowerBound() {
		return 0;
	}

	@Override
	public double getSupportUpperBound() {
		return Double.POSITIVE_INFINITY;
	}

	@Override
	public boolean isSupportLowerBoundInclusive() {
		return true;
	}

	@Override
	public boolean isSupportUpperBoundInclusive() {
		return false;
	}

	@Override
	public boolean isSupportConnected() {
		return true;
	}

	@Override
	public double getMean() {
		return this.getNumericalMean();
	}

	@Override
	public double getSD() {
		return Math.sqrt(this.getNumericalVariance());
	}

	@Override
	public double quantile(double q) {
		if (totalWeight == 0) return Double.NaN;
		if (!quantiles.containsKey(q)) {
			quantiles.put(q, this.inverseCumulativeProbability(q));
		}
		return quantiles.get(q);
	}
    
	public GammaParameters gammaApproximation() {
		return new GammaMoments(this.getMean(),this.getSD()).convert();
	}

}
