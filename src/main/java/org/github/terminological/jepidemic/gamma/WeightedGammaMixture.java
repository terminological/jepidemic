package org.github.terminological.jepidemic.gamma;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.util.Pair;

public class WeightedGammaMixture extends AbstractRealDistribution implements StatSummary {

	public static final double DEFAULT_INVERSE_ABSOLUTE_ACCURACY = 1e-9;
    
	private List<Pair<GammaParameters, Double>> weightedMap = new ArrayList<>();
	private double totalWeight = 0;
	private double mean = -1;
	private double var = -1;
	private Map<Double,Double> quantiles = new HashMap<>();
	private double acc = Double.NEGATIVE_INFINITY;
	
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
    		double ul = g.convert().getMean()+5*g.convert().getSD();
    		if (acc<ul) acc=ul;
    	});
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
	public double getNumericalMean() {
		if (totalWeight == 0) return Double.NaN;
		if (mean < 0) {
			mean = this.weightedMap.stream().mapToDouble(d -> d.getKey().getMean()*d.getValue()).sum()/totalWeight;
		}
		return mean;
	}

	@Override
	public double getNumericalVariance() {
		if (totalWeight == 0) return Double.NaN;
		if (var < 0) {
			var = this.weightedMap.stream().mapToDouble(d ->
				d.getValue()*( //weight * (
					d.getKey().getVar()+ // var_i   
					d.getKey().getMean()*d.getKey().getMean()) //+mu_i^2
				).sum()/totalWeight - this.getNumericalMean()*this.getNumericalMean();
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
