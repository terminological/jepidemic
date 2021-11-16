package org.github.terminological.jepidemic.distributions;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;

public class TruncatedIntegerDistribution extends AbstractIntegerDistribution {

	AbstractIntegerDistribution d;
	int ul;
	int ll;
	double mean = Double.NaN;
	double var = Double.NaN;
	boolean useApproximateMean;
	double up;
	double lp;
	
	
protected TruncatedIntegerDistribution(AbstractIntegerDistribution d, int ll, int ul, boolean useOriginalMean) {
	super(Random.RNG);
		this.d = d;
		this.useApproximateMean = useOriginalMean;
		this.changeLimits(ll, ul);
	}
	
	public TruncatedIntegerDistribution changeLimits(int newMin, int newMax) {
		if (newMin != ll || newMax != ul) {
			this.ll = Math.max(d.getSupportLowerBound(), newMin);
			this.ul = Math.min(d.getSupportUpperBound(), newMax);
			this.lp = d.cumulativeProbability(ll-1);
			this.up = d.cumulativeProbability(ul);
			
			this.mean = Double.NaN;
			this.var = Double.NaN;
		}
		return this;
	}
	
	public TruncatedIntegerDistribution enforceHardLimits(int hardLowerLimit, int hardUpperLimit) {
		if ( Double.isNaN(this.getMaxValue()) || this.getMaxValue() > hardUpperLimit)  
    		this.changeLimits(this.getMinValue(), hardUpperLimit);
    	
    	if ( Double.isNaN(this.getMinValue()) || this.getMinValue() < hardLowerLimit)
    		this.changeLimits(hardLowerLimit, this.getMaxValue());
    	
    	return this;
	}
	
	@Override
	public double probability(int x) {
		if (x<ll || x >ul) return 0;
		return d.probability(x)*1/(up-lp);
	}

	@Override
	public double cumulativeProbability(int x) {
		if (x<ll) return 0;
		if (x>ul) return 1;
		return (d.cumulativeProbability(x) - lp) / (up-lp); 
	}
	
	

	@Override
	public double getNumericalMean() {
		if (useApproximateMean) return d.getNumericalMean();
		if (Double.isNaN(mean)) {
			double sum = 0;
			double sum2 = 0;
			for (int i = ll; i<=ul; i++) {
				sum = sum+i*probability(i);
				sum2 = sum2+i*i*probability(i);
			}
			mean = sum;
			var = sum2 - mean*mean;
		}
		return mean;
	}

	@Override
	public double getNumericalVariance() {
		if (useApproximateMean) return d.getNumericalVariance();
		if (Double.isNaN(var)) {
			getNumericalMean();
		}
		return var;
	}

	@Override
	public int getSupportLowerBound() {
		return ll;
	}

	@Override
	public int getSupportUpperBound() {
		return ul;
	}

	@Override
	public boolean isSupportConnected() {
		return d.isSupportConnected();
	}
	
	public int getMinValue() {return ll;}
	public int getMaxValue() {return ul;}

	public static TruncatedIntegerDistribution of(TruncatedIntegerDistribution distribution, int low, int high) {
		return distribution.changeLimits(low, high);
	}

	public static TruncatedIntegerDistribution of(AbstractIntegerDistribution d) {
		int ll = (int) Math.floor(d.getNumericalMean() - 10* Math.sqrt(d.getNumericalVariance()));
		int ul = (int) Math.ceil(d.getNumericalMean() + 10* Math.sqrt(d.getNumericalVariance()));
		return of(d,ll,ul, true);
	}

	public static TruncatedIntegerDistribution of(AbstractIntegerDistribution distribution, int low, int high, boolean useOriginalMean) {
		return new TruncatedIntegerDistribution(distribution, low, high, useOriginalMean);
	}

	public static TruncatedIntegerDistribution of(TruncatedIntegerDistribution d) {
		return d;
	}

}
