package org.github.terminological.jepidemic.distributions;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.exception.MultiDimensionMismatchException;
import org.apache.commons.math3.util.Pair;

public class WeightedMixtureRealDistribution extends AbstractRealDistribution {

	public static final double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;
    
	private List<Pair<? extends TruncatedRealDistribution, Double>> weightedMap = new ArrayList<>();
	private double totalWeight = 0;
	private List<Double> breaks = null; 
	private Double mean = null;
	private Double var = null;
	private double max;
	private double min;
	//private UnivariateIntegrator ti;
	double hardLowerLimit;
	double hardUpperLimit;
	
    
    public WeightedMixtureRealDistribution(List<Pair<? extends AbstractRealDistribution,Double>> distributions, double hardLowerLimit, double hardUpperLimit) {
    	super(Random.RNG);
    	this.min = Double.MAX_VALUE;
    	this.max = Double.MIN_VALUE;
    	this.hardLowerLimit = hardLowerLimit;
    	this.hardUpperLimit = hardUpperLimit;
    	//this.ti = new RombergIntegrator(0.001,0.001,RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
    	breaks = new ArrayList<>();
    	totalWeight = 0;
    	distributions.stream().forEach(li -> this.add(li.getFirst(), li.getSecond()));
    }

    
    public WeightedMixtureRealDistribution add(AbstractRealDistribution distribution, Double weight) {
    	TruncatedRealDistribution trunc;
    	if(distribution instanceof TruncatedRealDistribution) {
    		trunc = (TruncatedRealDistribution) distribution;
    	} else {
    		trunc = TruncatedRealDistribution.of(distribution);
    	}
    	trunc.enforceHardLimits(hardLowerLimit, hardUpperLimit);
    	
    	this.weightedMap.add(Pair.create(trunc, weight));
    	min = Math.min(min, trunc.getMinValue());
    	max = Math.max(max, trunc.getMaxValue());
//    	min = Math.min(min, distribution.getSupportLowerBound());
//    	max = Math.max(max, distribution.getSupportUpperBound());
    	mean = null;
    	var = null;
    	totalWeight = totalWeight+weight;
    	breaks.add(totalWeight);
    	return this;
    }
    
	@Override
	public double density(double x) {
		if (totalWeight == 0) return Double.NaN;
		if (x<=min || x >=max) return 0;
		double sum =  weightedMap.stream()
				.mapToDouble(kv -> kv.getKey().density(x)*kv.getValue()/totalWeight)
				.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
				.sum();
		return sum;
		
	}

	@Override
	public double cumulativeProbability(double x) {
		if (totalWeight == 0) return Double.NaN;
		if (x<=min) return 0;
		if (x>=max) return 1;
		double out = weightedMap
				.stream()
				.mapToDouble(
					kv -> kv.getKey().cumulativeProbability(x)*kv.getValue()/totalWeight
				)
				.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
				.sum();
		return out;
	}

	@Override
	public double sample() {
		double selector = this.random.nextDouble()*totalWeight;
		int selection = 0;
		while (selector > breaks.get(selection) ) selection +=1;
		return weightedMap.get(selection).getFirst().sample();
	}

	@Override
	public double[] sample(int sampleSize) {
		// TODO: better way to do this exists by determining size of component samples first rather than repeated sampling.
		double[] out = new double[sampleSize];
		for (int i=0; i<sampleSize; i++) {
			out[i] = sample();
		}
		return out;
	}
	
	@Override
	// https://en.wikipedia.org/wiki/Mixture_distribution
	public double getNumericalMean() {
		if (totalWeight == 0) return Double.NaN;
		if (mean == null) {
			Double sum = this.weightedMap.stream()
					.map(x -> x.getSecond()/totalWeight*x.getFirst().getNumericalMean())
					.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
					.reduce(Double::sum).orElse(Double.NaN);
			mean = sum;
//			UnivariateFunction f = new UnivariateFunction() {
//				@Override
//				public double value(double x) {
//					return x * WeightedMixtureRealDistribution.this.density(x);
//				}};
//			mean = this.ti.integrate(100000, f, min, max);
		}
		return mean;
	}

	@Override
	// https://en.wikipedia.org/wiki/Mixture_distribution
	public double getNumericalVariance() {
		if (totalWeight == 0) return Double.NaN;
		if (var == null) {
			Double sum = this.weightedMap.stream()
					.map(x -> x.getSecond()/totalWeight*(x.getFirst().getNumericalVariance()+Math.pow(x.getFirst().getNumericalMean(),2)))
					.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
					.reduce(Double::sum).orElse(Double.NaN);
			var = sum-Math.pow(getNumericalMean(),2);
//			UnivariateFunction f = new UnivariateFunction() {
//				@Override
//				public double value(double x) {
//					return Math.pow(x,2) * WeightedMixtureRealDistribution.this.density(x);
//				}};
//			var = this.ti.integrate(100000, f, min, max) - Math.pow(this.getNumericalMean(),2);
//			// https://en.wikipedia.org/wiki/Variance#Absolutely_continuous_random_variable
		}
		return var;
	}

	@Override
	public double getSupportLowerBound() {
		return min;
	}

	@Override
	public double getSupportUpperBound() {
		return max;
	}

	@Override
	public boolean isSupportLowerBoundInclusive() {
		return false;
	}

	@Override
	public boolean isSupportUpperBoundInclusive() {
		return false;
	}

	@Override
	public boolean isSupportConnected() {
		return true;
	}
	
	public static WeightedMixtureRealDistribution empty(double hardLowerLimit, double hardUpperLimit) {
		return new WeightedMixtureRealDistribution(new ArrayList<>(), hardLowerLimit, hardUpperLimit);
	}
	
	public static WeightedMixtureRealDistribution empty() {
		return empty(-Double.MAX_VALUE,Double.MAX_VALUE);
	}


	public Optional<AbstractRealDistribution> ifNotEmpty() {
		if(weightedMap.isEmpty()) return Optional.empty();
		return Optional.of(this);
	}
	
	/**
     * Creates a Mixture of distributions.
     *
     * @param <X> - Dated or undated gamma parameters
     * @param gammas - a list of the gamma distributions of the mixture
     */
    public static WeightedMixtureRealDistribution unweighted(Collection<? extends AbstractRealDistribution> distributions, double hardLowerLimit, double hardUpperLimit) {
    	return(
			new WeightedMixtureRealDistribution(
					distributions.stream().map(x -> (Pair<? extends AbstractRealDistribution, Double>) Pair.create(x, 1.0)).collect(Collectors.toList()),
					hardLowerLimit, hardUpperLimit
			)
    	);
    }
    
    /**
     * Creates a Mixture of distributions.
     *
     * @param <X> - Dated or undated gamma parameters
     * @param gammas - a list of the gamma distributions of the mixture
     */
    public static WeightedMixtureRealDistribution weighted(List<? extends AbstractRealDistribution> distributions, List<Double> weights, double hardLowerLimit, double hardUpperLimit) {
    	if (weights.size() != distributions.size()) throw new MultiDimensionMismatchException(new Integer[] {distributions.size()}, new Integer[] {weights.size()});
    	List<Pair<? extends AbstractRealDistribution,Double>> out = new ArrayList<>();
    	int i =0;
    	for (AbstractRealDistribution distribution: distributions) {
    		out.add(Pair.create(distribution, weights.get(i)));
    		i++;
    	}
    	return(new WeightedMixtureRealDistribution(out, hardLowerLimit, hardUpperLimit));
    }


	
}
