package org.github.terminological.jepidemic.distributions;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.exception.MultiDimensionMismatchException;
import org.apache.commons.math3.util.Pair;

public class WeightedMixtureIntegerDistribution extends AbstractIntegerDistribution {

	public static final double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;
    
	private List<Pair<TruncatedIntegerDistribution, Double>> weightedMap = new ArrayList<>();
	private double totalWeight = 0;
	private List<Double> breaks = null; 
	private Double mean = null;
	private Double var = null;
	private int max;
	private int min;
	
	
    
    public WeightedMixtureIntegerDistribution(List<Pair<? extends AbstractIntegerDistribution,Double>> distributions) {
    	super(Random.RNG);
    	this.min = Integer.MAX_VALUE;
    	this.max = Integer.MIN_VALUE;
    	breaks = new ArrayList<>();
    	totalWeight = 0;
    	distributions.stream().forEach(li -> this.add(li.getFirst(), li.getSecond()));
    }

    
    public WeightedMixtureIntegerDistribution add(AbstractIntegerDistribution distribution, Double weight) {
    	TruncatedIntegerDistribution trunc;
    	if(distribution instanceof TruncatedIntegerDistribution) {
    		trunc = (TruncatedIntegerDistribution) distribution;
    	} else {
    		trunc = TruncatedIntegerDistribution.of(distribution);
    	}
    	this.weightedMap.add(Pair.create(trunc, weight));
    	min = Math.min(min, trunc.getMinValue());
    	max = Math.max(max, trunc.getMaxValue());
    	mean = null;
    	var = null;
    	totalWeight = totalWeight+weight;
    	breaks.add(totalWeight);
    	return this;
    }
    
	@Override
	public double probability(int x) {
		if (totalWeight == 0) return Double.NaN;
		if (x<=min || x >=max) return 0;
		return weightedMap.stream()
				.mapToDouble(kv -> kv.getKey().probability(x)*kv.getValue()/totalWeight)
				.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
				.sum();
	}

	@Override
	public double cumulativeProbability(int x) {
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
	public int sample() {
		double selector = this.random.nextDouble()*totalWeight;
		int selection = 0;
		while (selector > breaks.get(selection) ) selection +=1;
		return weightedMap.get(selection).getFirst().sample();
	}

	@Override
	public int[] sample(int sampleSize) {
		// TODO: better way to do this exists by determining size of component samples first rather than repeated sampling.
		int[] out = new int[sampleSize];
		for (int i=0; i<sampleSize; i++) {
			out[i] = sample();
		}
		return out;
	}
	
	@Override
	public double getNumericalMean() {
		if (mean == null) {
			Double sum = this.weightedMap.stream()
					.map(x -> x.getSecond()/totalWeight*x.getFirst().getNumericalMean())
					.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
					.reduce(Double::sum).orElse(Double.NaN);
			mean = sum;
//			double sum = 0;
//			double sum2 = 0;
//			for (int i = min; i<=max; i++) {
//				sum = sum+i*probability(i);
//				sum2 = sum2+i*i*probability(i);
//			}
//			mean = sum;
//			var = sum2 - mean*mean;
		}
		return mean;
	}

	@Override
	public double getNumericalVariance() {
		if (var == null) {
			// getNumericalMean();
			Double sum = this.weightedMap.stream()
					.map(x -> x.getSecond()/totalWeight*(x.getFirst().getNumericalVariance()+Math.pow(x.getFirst().getNumericalMean(),2)))
					.filter(v -> !Double.isNaN(v) && Double.isFinite(v))
					.reduce(Double::sum).orElse(Double.NaN);
			var = sum-Math.pow(getNumericalMean(),2);
		}
		return var;
	}

	@Override
	public int getSupportLowerBound() {
		return min;
	}

	@Override
	public int getSupportUpperBound() {
		return max;
	}

	@Override
	public boolean isSupportConnected() {
		return true;
	}
	
	public static WeightedMixtureIntegerDistribution empty() {
		return new WeightedMixtureIntegerDistribution(new ArrayList<>());
	}

	public Optional<AbstractIntegerDistribution> ifNotEmpty() {
		if(weightedMap.isEmpty()) return Optional.empty();
		return Optional.of(this);
	}
	
	/**
     * Creates a Mixture of distributions.
     *
     * @param distributions - a collection of the distributions of the mixture
     * @return a weighted mixture distribtion
     */
    public static WeightedMixtureIntegerDistribution unweighted(Collection<? extends AbstractIntegerDistribution> distributions) {
    	return(
			new WeightedMixtureIntegerDistribution(
					distributions.stream().map(x -> (Pair<? extends AbstractIntegerDistribution, Double>) Pair.create(x, 1.0)).collect(Collectors.toList())
			)
    	);
    }
    
    /**
     * Creates a Weighted mixture of distributions.
     *
     * @param distributions - a list of the distributions of the mixture
     * @param weights - a list of the distribution weights of the mixture
     * @return a weighted mixture distribtion
     */
    public static WeightedMixtureIntegerDistribution weighted(List<? extends AbstractIntegerDistribution> distributions, List<Double> weights) {
    	if (weights.size() != distributions.size()) throw new MultiDimensionMismatchException(new Integer[] {distributions.size()}, new Integer[] {weights.size()});
    	List<Pair<? extends AbstractIntegerDistribution,Double>> out = new ArrayList<>();
    	int i =0;
    	for (AbstractIntegerDistribution distribution: distributions) {
    		out.add(Pair.create(distribution, weights.get(i)));
    		i++;
    	}
    	return(new WeightedMixtureIntegerDistribution(out));
    }

}
