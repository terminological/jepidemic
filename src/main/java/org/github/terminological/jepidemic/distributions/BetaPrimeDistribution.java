package org.github.terminological.jepidemic.distributions;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.Precision;

public class BetaPrimeDistribution extends AbstractRealDistribution {

	private static final double DEFAULT_EPSILON = 1E-14;
	
	double shape1; // on wolfram this is p, on wikipedia alpha
	double shape2; // on wolfram this is q, on wikipedia beta
	double scale; // on wolfram this is beta (parameter 3), on wikipedia this is q (parameter 4)
	double shape3; // on wolfram this is alpha (parameter 4) on wikidepia p (parameter3)
	
	
	private BetaPrimeDistribution(double shape1, double shape2) {
		this(shape1,shape2,1D,1D);
	}
	
	private BetaPrimeDistribution(double shape1, double shape2, double scale) {
		this(shape1,shape2,1D,scale);
	}
	
	private BetaPrimeDistribution(double shape1, double shape2, double shape3, double scale) {
		super(Random.RNG);
		if(shape1 <= 0) throw new NotStrictlyPositiveException(shape1);
		if(shape2 <= 0) throw new NotStrictlyPositiveException(shape2);
		if(scale <=0) throw new NotStrictlyPositiveException(scale);
		if(shape3 <=0) throw new NotStrictlyPositiveException(shape3);
		this.shape1 = shape1;
		this.shape2 = shape2;
		this.shape3 = shape3;
		this.scale = scale;
		
	}

	@Override
	public double density(double x) {
		// https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization
		// the generalised form but parameterised differently
		// https://reference.wolfram.com/language/ref/BetaPrimeDistribution.html
		if (x < getSupportLowerBound()) return 0;
		if (x > getSupportUpperBound()) return 0;
		double p = shape1;
		double q = shape2;
		//double alpha = scale;
		//double beta = shape3;
		double beta = scale;
		double alpha = shape3;
				
		return alpha * Math.pow((1+Math.pow(x/beta, alpha)), -p-q) * Math.pow(x/beta, p*alpha-1)  /
				(beta*Math.exp(Beta.logBeta(p, q)));
	}

	@Override
	public double cumulativeProbability(double x) {
		if (x < getSupportLowerBound()) return 0;
		if (x > getSupportUpperBound()) return 1;
		double p = shape1;
		double q = shape2;
		//double alpha = scale;
		//double beta = shape3;
		double beta = scale;
		double alpha = shape3;
		
		double xPrime = Math.pow(x, alpha) / (Math.pow(x, alpha)+Math.pow(beta, alpha));
		double test = Math.abs(p-q)/p;
		
		
		try {
			if(xPrime==0) return 0;
			if(xPrime==1) return 1;
			if(Precision.equals(test,0,DEFAULT_EPSILON) && Precision.equals(xPrime,0.5,DEFAULT_EPSILON)) {
				return 0.5;
			}
			return Beta.regularizedBeta(
					xPrime, 
					p, q, 10000);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public double getNumericalMean() {
		double p = shape1; // alpha on wikipedia
		double q = shape2; // beta on wikipedia
		//double alpha = scale;
		//double beta = shape3;
		double beta = scale; // q on wikipedia
		double alpha = shape3; // p on wikipedia
		
		
		if (q < 1/alpha) return Double.POSITIVE_INFINITY;
		
		double g1 = Gamma.logGamma(q - 1/alpha);
		double g2 = Gamma.logGamma(p + 1/alpha);
		double g3 = Gamma.logGamma(p);
		double g4 = Gamma.logGamma(q);
		
		return beta * Math.exp(g1 + g2 - g3 - g4);
	}

	@Override
	public double getNumericalVariance() {
		double p = shape1;
		double q = shape2;
		//double alpha = scale;
		//double beta = shape3;
		double beta = scale;
		double alpha = shape3;
		
		if (q < 2/alpha) return Double.NaN;
		
		double g1 = 2*Gamma.logGamma(q-1/alpha)+2*Gamma.logGamma(p+1/alpha);
		double divisor = 2*Gamma.logGamma(p)+2*Gamma.logGamma(q);
		double g2 = Gamma.logGamma(p)+Gamma.logGamma(q)+Gamma.logGamma(q-2/alpha)+Gamma.logGamma(p+2/alpha);
		
		return 
				Math.pow(beta,2) * 
				(
					-Math.exp(
							g1
							-divisor
					)
					+Math.exp(
							g2
							-divisor
					)
				);
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
		return false;
	}

	@Override
	public boolean isSupportUpperBoundInclusive() {
		return true;
	}

	@Override
	public boolean isSupportConnected() {
		return true;
	}

	
	// wikipedia gives this in the order shape1 (alpha), shape2 (beta), shape3 (p), scale (q)
	public static BetaPrimeDistribution wikipediaParameters(double alpha, double beta, double p, double q) {
		if(Double.isNaN(alpha) || Double.isNaN(beta) || Double.isNaN(p) || Double.isNaN(q)) return null;
		return new BetaPrimeDistribution(alpha,beta,q,p);
	}
	
	public static BetaPrimeDistribution wikipediaParameters(double alpha, double beta) {
		if(Double.isNaN(alpha) || Double.isNaN(beta)) return null;
		return new BetaPrimeDistribution(alpha,beta,1,1);
	}
	
	
	// alpha is shape3; beta is scale
	public static BetaPrimeDistribution wolframParameters(double p, double q, double alpha, double beta) {
		if(Double.isNaN(alpha) || Double.isNaN(beta) || Double.isNaN(p) || Double.isNaN(q)) return null;
		return new BetaPrimeDistribution(p,q,beta,alpha);
	}
	
	// alpha is shape3 = 1; beta is scale
	public static BetaPrimeDistribution wolframParameters(double p, double q, double beta) {
		if(Double.isNaN(beta) || Double.isNaN(p) || Double.isNaN(q)) return null;
		return new BetaPrimeDistribution(p,q,beta,1);
	}
	
	public static BetaPrimeDistribution wolframParameters(double p, double q) {
		if(Double.isNaN(p) || Double.isNaN(q)) return null;
		return new BetaPrimeDistribution(p,q,1,1);
	}
	
	public static BetaPrimeDistribution extraDistrParameters(double shape1, double shape2, double scale) {
		if(Double.isNaN(shape1) || Double.isNaN(shape2) || Double.isNaN(scale)) return null;
		return new BetaPrimeDistribution(shape1,shape2,scale);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(scale);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(shape1);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(shape2);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(shape3);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof BetaPrimeDistribution))
			return false;
		BetaPrimeDistribution other = (BetaPrimeDistribution) obj;
		if (Double.doubleToLongBits(scale) != Double.doubleToLongBits(other.scale))
			return false;
		if (Double.doubleToLongBits(shape1) != Double.doubleToLongBits(other.shape1))
			return false;
		if (Double.doubleToLongBits(shape2) != Double.doubleToLongBits(other.shape2))
			return false;
		if (Double.doubleToLongBits(shape3) != Double.doubleToLongBits(other.shape3))
			return false;
		return true;
	}
	
	
	
}
