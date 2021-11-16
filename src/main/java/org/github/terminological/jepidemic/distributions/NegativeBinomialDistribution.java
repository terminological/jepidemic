package org.github.terminological.jepidemic.distributions;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

public class NegativeBinomialDistribution extends AbstractIntegerDistribution {

	private NegativeBinomialDistribution(double size, double probability) {
		super(Random.RNG);
		this.size = size;
		this.probability = probability;
	}

	double size;
	double probability;
	
	@Override
	public double getNumericalMean() {
		return probability*size / (1- probability);
	}

	@Override
	public double getNumericalVariance() {
		return probability*size / Math.pow(1- probability,2);
	}

	@Override
	public boolean isSupportConnected() {
		return true;
	}

	@Override
	// https://en.wikipedia.org/wiki/Negative_binomial_distribution#Gamma%E2%80%93Poisson_mixture
	public double probability(int x) {
		if (x < getSupportLowerBound() || x > getSupportUpperBound()) return 0;
		try {
			return Gamma.gamma(size+x) / 
				(CombinatoricsUtils.factorial(x) * Gamma.gamma(size)) *
				Math.pow(probability, x) *
				Math.pow(1-probability,size);
		} catch (MathArithmeticException e) {
			return 0;
		}
	}

	@Override
	// https://en.wikipedia.org/wiki/Negative_binomial_distribution
	public double cumulativeProbability(int x) {
		if (x < getSupportLowerBound()) return 0;
		if (x > getSupportUpperBound()) return 1;
		double out = 1-Beta.regularizedBeta(probability, x+1, size);
		if (Double.isNaN(out)) {
			return 1;
		}
		return out;
	}

	@Override
	public int getSupportLowerBound() {
		return 0;
	}

	@Override
	public int getSupportUpperBound() {
		return Integer.MAX_VALUE;
	}
	

	public NegativeBinomialDistribution wider(double factor) {
		Summary summ = Summary.of(this);
		return NegativeBinomialDistribution.fromMoments(summ.getMean(), summ.getSD()*factor);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(probability);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(size);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof NegativeBinomialDistribution))
			return false;
		NegativeBinomialDistribution other = (NegativeBinomialDistribution) obj;
		if (Double.doubleToLongBits(probability) != Double.doubleToLongBits(other.probability))
			return false;
		if (Double.doubleToLongBits(size) != Double.doubleToLongBits(other.size))
			return false;
		return true;
	}
	

	public static NegativeBinomialDistribution fromMoments(double mean, double sd) {
		if(Double.isNaN(mean) || Double.isNaN(sd)) return null;
		return new NegativeBinomialDistribution(
			FastMath.pow(mean,2) / (FastMath.pow(sd,2)-mean), //rate 
			1-mean/FastMath.pow(sd,2) //probability
		);
	}
	
	public static NegativeBinomialDistribution fromSizeAndProbability(double size, double probability) {
		if(Double.isNaN(size) || Double.isNaN(probability)) return null;
		return new NegativeBinomialDistribution(size, probability);
	}
	
}
