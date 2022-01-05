package org.github.terminological.jepidemic.distributions;

public class GrowthRateDistribution extends TransformerDistribution<BetaPrimeDistribution> {

	double alphaMinus; 
	double alphaPlus;
	int tau;
	
	public GrowthRateDistribution(double alphaPlus, double alphaMinus, int tau) {
		super(
				BetaPrimeDistribution.extraDistrParameters(alphaPlus, alphaMinus, 1),
				/*y=g(x)*/ x -> Math.log(x) / (2.0*tau),
				/*x=g^-1(y)*/ y -> Math.exp(2.0*tau*y),
				/*dg-1/dy =*/ y -> 2.0*tau*Math.exp(2.0*tau*y), 
				TransformerDistribution.MEAN_PRECISION, TransformerDistribution.SD_PRECISION);
		this.alphaPlus = alphaPlus;
		this.alphaMinus = alphaMinus;
		this.tau = tau;
	}

	@Override
	public double getNumericalMean() {
		
		double betap_mean = this.getBareDistribution().getNumericalMean();
		double betap_variance = this.getBareDistribution().getNumericalVariance();
		// https://stats.stackexchange.com/questions/57715/expected-value-and-variance-of-loga
		// taylor exansion for log
		double mean = 1/(2.0*tau)*(Math.log(betap_mean) - betap_variance/(2.0*Math.pow(betap_mean,2)));
		return mean;
	}



	@Override
	public double getNumericalVariance() {
		
		double betap_mean = this.getBareDistribution().getNumericalMean();
		double betap_variance = this.getBareDistribution().getNumericalVariance();
		// https://stats.stackexchange.com/questions/57715/expected-value-and-variance-of-loga
		// taylor exansion for log
		 double variance = 1/(4.0*Math.pow(tau,2))*(betap_variance/Math.pow(betap_mean,2));
		return variance;
	}
	
//	//public class GrowthRateDistribution extends TruncatedRealDistribution {
//
//
//	//	private GrowthRateDistribution(GrowthRateDistribution.ExponentialBetaPrime dist, double lower, double upper) {
//	//		super(dist,lower,upper);
//	//	}
//	//	
//	//	public static GrowthRateDistribution withShapeAndShape(double alphaPlus, double alphaMinus, int tau) {
//	//		if(Double.isNaN(alphaPlus) || Double.isNaN(alphaMinus)) return null;
//	//		ExponentialBetaPrime dist = new GrowthRateDistribution.ExponentialBetaPrime(alphaPlus, alphaMinus, tau);
//	//		double lower = dist.inverseCumulativeProbability(0.0001);
//	//		double upper = dist.inverseCumulativeProbability(0.9999);
//	//		return new GrowthRateDistribution(dist,lower,upper);
//	//	}
//	//	
//	//	public static class GrowthRateDistribution extends AbstractRealDistribution {
//
//	double shape1; // on wolfram this is p, on wikipedia alpha
//	double shape2; // on wolfram this is q, on wikipedia beta
//	double tau;
//	UnivariateIntegrator ti;
//	Double mean = null;
//	Double var = null;
//	//		BrentSolver s;
////	private double ll;
////	private double ul;
//	private double lp;
//	private double up;
//	//		double LOWER_TRUNC = 0.00001;
//	//		double UPPER_TRUNC = 0.99999;
//	//		
//	// private ExponentialBetaPrime(double shape1, double shape2, int tau) {
//	private GrowthRateDistribution(double shape1, double shape2, int tau) {
//		super(new Well19937c());
//		if(shape1 <= 0) throw new NotStrictlyPositiveException(shape1);
//		if(shape2 <= 0) throw new NotStrictlyPositiveException(shape2);
//		this.shape1 = shape1;
//		this.shape2 = shape2;
//		this.tau = (double) tau;
//		this.up = rawCumulativeProbability(this.getSupportUpperBound());
//		this.lp = rawCumulativeProbability(this.getSupportLowerBound());
////		ul = inverseCumulativeProbability(up);
////		//			s = new BrentSolver();
//		//			ll = s.solve(10000, x -> actualCumulativeProbability(x)-LOWER_TRUNC, -Double.MAX_VALUE, Double.MAX_VALUE, 0);
//		//			ul = s.solve(10000, x -> actualCumulativeProbability(x)-UPPER_TRUNC, -Double.MAX_VALUE, Double.MAX_VALUE, 0);
//		this.ti = new RombergIntegrator();
//	}
//
//	@Override
//	public double density(double x) {
//		// https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
//		// the generalised form but parameterised differently
//		// https://reference.wolfram.com/language/ref/BetaPrimeDistribution.html
//		if (x<=this.getSupportLowerBound() || x >=this.getSupportUpperBound()) return 0;
//		return rawDensity(x) / (up-lp);
//	}
//	
//	public double rawDensity(double x) {
//		// https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
//		// the generalised form but parameterised differently
//		// https://reference.wolfram.com/language/ref/BetaPrimeDistribution.html
//		double p = shape1;
//		double q = shape2;
//		return Math.pow((1+Math.exp(x)), -p-q) * Math.exp(x*(p+2*tau-1)) * 
//				2*tau /	Math.exp(Beta.logBeta(p, q)); //Math.exp(Beta.logBeta(p, q));
//	}
//
//	@Override
//	public double cumulativeProbability(double x) {
//		if (x<=this.getSupportLowerBound()) return 0;
//		if (x>=this.getSupportUpperBound()) return 1;
//		return (rawCumulativeProbability(x)-lp)/(up-lp);
//	}
//	
//	public double rawCumulativeProbability(double x) {
//		double p = shape1;
//		double q = shape2;
//		double x2 = Math.exp(2*tau*x) / (1+Math.exp(2*tau*x));
//		if (Double.isNaN(x2)) x2=1;
//		double out = Beta.regularizedBeta(x2, p, q);
//		return out;
//	}
//
//	@Override
//	public double getNumericalMean() {
//		if (mean == null) {
//			UnivariateFunction f = new UnivariateFunction() {
//				@Override
//				public double value(double x) {
//					return x * density(x);
//				}};
//			mean = this.ti.integrate(100000, f, getSupportLowerBound(), getSupportUpperBound());
//		}
//		return mean;
//	}
//
//	@Override
//	public double getNumericalVariance() {
//		if (var == null) {
//			UnivariateFunction f = new UnivariateFunction() {
//				@Override
//				public double value(double x) {
//					return x * x * density(x);
//				}};
//				var = this.ti.integrate(100000, f, getSupportLowerBound(), getSupportUpperBound()) - Math.pow(this.getNumericalMean(),2);
//		}
//		return var;
//	}
//
//	@Override
//	public double getSupportLowerBound() {
//		return -5;
//	}
//
//	@Override
//	public double getSupportUpperBound() {
//		return +5;
//	}
//
//	@Override
//	public boolean isSupportLowerBoundInclusive() {
//		return false;
//	}
//
//	@Override
//	public boolean isSupportUpperBoundInclusive() {
//		return false;
//	}
//
//	@Override
//	public boolean isSupportConnected() {
//		return true;
//	}
//
//	@Override
//	public int hashCode() {
//		final int prime = 31;
//		int result = 1;
//		long temp;
//		temp = Double.doubleToLongBits(shape1);
//		result = prime * result + (int) (temp ^ (temp >>> 32));
//		temp = Double.doubleToLongBits(shape2);
//		result = prime * result + (int) (temp ^ (temp >>> 32));
//		temp = Double.doubleToLongBits(tau);
//		result = prime * result + (int) (temp ^ (temp >>> 32));
//		return result;
//	}
//
//	@Override
//	public boolean equals(Object obj) {
//		if (this == obj)
//			return true;
//		if (obj == null)
//			return false;
//		if (!(obj instanceof GrowthRateDistribution))
//			return false;
//		GrowthRateDistribution other = (GrowthRateDistribution) obj;
//		if (Double.doubleToLongBits(shape1) != Double.doubleToLongBits(other.shape1))
//			return false;
//		if (Double.doubleToLongBits(shape2) != Double.doubleToLongBits(other.shape2))
//			return false;
//		if (Double.doubleToLongBits(tau) != Double.doubleToLongBits(other.tau))
//			return false;
//		return true;
//	}

	public static GrowthRateDistribution withForwardEstimateAndBackwardEstimate(double alphaPlus, double alphaMinus, int tau) {
		if(Double.isNaN(alphaPlus) || Double.isNaN(alphaMinus)) return null;
		return new GrowthRateDistribution(alphaPlus, alphaMinus, tau);
	}



	
	

}
