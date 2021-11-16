package org.github.terminological.jepidemic.distributions;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.exception.OutOfRangeException;

public class TruncatedRealDistribution extends AbstractRealDistribution {

	AbstractRealDistribution d;
	double ul;
	double ll;
//	double LOWER_TRUNC = 0.00001;
	BracketingNthOrderBrentSolver s;
	double mean = Double.NaN;
	double var = Double.NaN;
//	double uSupport = Double.NaN;
//	double lSupport = Double.NaN;
	boolean useApproximateMean;
	RombergIntegrator ti;
	double up;
	double lp;
	
	
//	double TOL = 10;
	
//	public double getActualSupportUpperBound() {
//		if (Double.isNaN(uSupport)) {
//			double hi = Double.isInfinite(d.getSupportUpperBound()) ? Double.MAX_VALUE : Math.min(Double.MAX_VALUE, d.getSupportUpperBound());
//			double lo = Double.isInfinite(d.getSupportLowerBound()) ? -Double.MAX_VALUE : Math.max(-Double.MAX_VALUE, d.getSupportLowerBound());
//			while (Math.abs(hi-lo) > TOL) {
//				double test = (hi+lo)/2;
//				try {
//					double value = d.cumulativeProbability(test);
//					if (Double.isNaN(value)) {
//						hi = test;
//					} else {
//						lo = test;
//					}
//				} catch (Exception e) {
//					hi = test;
//				}
//			}
//			uSupport= lo;
//		}
//		return uSupport;
//	}
//	
//	public double getActualSupportLowerBound() {
//		if (Double.isNaN(lSupport)) {
//			double hi = Double.isInfinite(d.getSupportUpperBound()) ? Double.MAX_VALUE : Math.min(Double.MAX_VALUE, d.getSupportUpperBound());
//			double lo = Double.isInfinite(d.getSupportLowerBound()) ? -Double.MAX_VALUE : Math.max(-Double.MAX_VALUE, d.getSupportLowerBound());
//			while (Math.abs(hi-lo) > TOL) {
//				double test = (hi+lo)/2;
//				try {
//					if (Double.isNaN(d.cumulativeProbability(test))) {
//						lo = test;
//					} else {
//						hi = test;
//					}
//				} catch (Exception e) {
//					lo = test;
//				}
//			}
//			lSupport= hi;
//		}
//		return lSupport;
//	}
	
	protected TruncatedRealDistribution(AbstractRealDistribution d, double ll, double ul, boolean useOriginalMean) {
		super(Random.RNG);
		this.d = d;
		this.s = new BracketingNthOrderBrentSolver();
		this.ti = new RombergIntegrator(0.001,0.001,RombergIntegrator.DEFAULT_MIN_ITERATIONS_COUNT,RombergIntegrator.ROMBERG_MAX_ITERATIONS_COUNT);
		this.useApproximateMean = useOriginalMean;
		this.changeLimits(ll, ul);
		this.enforceHardLimits(-Double.MAX_VALUE, Double.MAX_VALUE);
		if( Double.isNaN(this.ll) || Double.isInfinite(this.ll) || Double.isNaN(this.ul) || Double.isInfinite(this.ul) ) {
			throw new RuntimeException("Unspecified or non finite limits in truncated distribution.");
		}
	}
	
	public TruncatedRealDistribution changeLimits(double newMin, double newMax) {
		if (newMin != ll || newMax != ul) {
			this.ll = Math.max(d.getSupportLowerBound(), newMin);
			this.ul = Math.min(d.getSupportUpperBound(), newMax);
			this.lp = d.cumulativeProbability(ll);
			this.up = d.cumulativeProbability(ul);
			
			this.mean = Double.NaN;
			this.var = Double.NaN;
		}
		return this;
	}
	
	public TruncatedRealDistribution enforceHardLimits(double hardLowerLimit, double hardUpperLimit) {
		if ( Double.isNaN(this.getMaxValue()) || this.getMaxValue() > hardUpperLimit)  
    		this.changeLimits(this.getMinValue(), hardUpperLimit);
    	
    	if ( Double.isNaN(this.getMinValue()) || this.getMinValue() < hardLowerLimit)
    		this.changeLimits(hardLowerLimit, this.getMaxValue());
    	
    	return this;
	}
	
	@Override
	public double inverseCumulativeProbability(double p) throws OutOfRangeException {
		return s.solve(Integer.MAX_VALUE, x -> cumulativeProbability(x)-p, ll,ul,(ll+ul)/2);
	}
	
	@Override
	public double density(double x) {
		if (x<=ll || x >=ul) return 0;
		return d.density(x)*1/(up-lp);
	}

	@Override
	public double cumulativeProbability(double x) {
		if (x<=ll) return 0;
		if (x>=ul) return 1;
		return (d.cumulativeProbability(x) - lp) / (up-lp); 
	}

	@Override
	public double getNumericalMean() {
		if (useApproximateMean) return d.getNumericalMean();
		if (Double.isNaN(mean)) {
			UnivariateFunction f = new UnivariateFunction() {
				@Override
				public double value(double x) {
					return x * d.density(x);
				}};
			mean = this.ti.integrate(100000, f, ll, ul);
		}
		return mean;
	}

	@Override
	public double getNumericalVariance() {
		if (useApproximateMean) return d.getNumericalVariance();
		if (Double.isNaN(var)) {
			UnivariateFunction f = new UnivariateFunction() {
				@Override
				public double value(double x) {
					return x * x * d.density(x);
				}};
			var = this.ti.integrate(100000, f, ll, ul) - Math.pow(this.getNumericalMean(),2);
		}
		return var;
	}

	@Override
	public double getSupportLowerBound() {
		return ll;
	}

	@Override
	public double getSupportUpperBound() {
		return ul;
	}

	@SuppressWarnings("deprecation")
	@Override
	public boolean isSupportLowerBoundInclusive() {
		return d.isSupportLowerBoundInclusive();
	}

	@SuppressWarnings("deprecation")
	@Override
	public boolean isSupportUpperBoundInclusive() {
		return d.isSupportUpperBoundInclusive();
	}

	@Override
	public boolean isSupportConnected() {
		return d.isSupportConnected();
	}
	
	public double getMinValue() {return ll;}
	public double getMaxValue() {return ul;}

	public static TruncatedRealDistribution of(TruncatedRealDistribution distribution, double low, double high) {
		return distribution.changeLimits(low, high);
	}

	public static TruncatedRealDistribution of(AbstractRealDistribution d) {
		// double ll = d.getNumericalMean() - 10* Math.sqrt(d.getNumericalVariance());
		// double ul = d.getNumericalMean() + 10* Math.sqrt(d.getNumericalVariance());
		
		double ll = Double.NaN;
		try {
			ll = d.inverseCumulativeProbability(0.0001);
		} catch (Exception e) {
			// Not a well behaved distribution - which is why we are truncating it anyway
		}			
		if (!Double.isFinite(ll) || Double.isNaN(ll)) ll = d.getNumericalMean() - 10* Math.sqrt(d.getNumericalVariance());
		
		double ul = Double.NaN;
		try {
			ul= d.inverseCumulativeProbability(0.9999);
		} catch (Exception e) {
			// Not a well behaved distribution - which is why we are truncating it anyway
		}
		if (!Double.isFinite(ul) || Double.isNaN(ul)) ul = d.getNumericalMean() + 10* Math.sqrt(d.getNumericalVariance());
		
		return of(d,ll,ul, true);
	}

	public static TruncatedRealDistribution of(AbstractRealDistribution distribution, double low, double high, boolean useOriginalMean) {
		return new TruncatedRealDistribution(distribution, low, high, useOriginalMean);
	}

	public static TruncatedRealDistribution of(TruncatedRealDistribution d) {
		return d;
	}

	public RealDistribution getUnderlyingDistribution() {
		return d;
	}

}
