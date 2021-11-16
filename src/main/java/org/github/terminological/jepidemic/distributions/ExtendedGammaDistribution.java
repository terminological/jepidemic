package org.github.terminological.jepidemic.distributions;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.util.FastMath;

public class ExtendedGammaDistribution extends GammaDistribution {

	private ExtendedGammaDistribution(double shape, double scale) throws NotStrictlyPositiveException {
		super(shape, scale);
	}
	
	public double getRate() {
		return 1/getScale();
	}
	
	
	public ExtendedGammaDistribution wider(double factor) {
		Summary summ = Summary.of(this);
		return ExtendedGammaDistribution.fromMoments(summ.getMean(), summ.getSD()*factor);
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(getScale());
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(getShape());
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof ExtendedGammaDistribution))
			return false;
		ExtendedGammaDistribution other = (ExtendedGammaDistribution) obj;
		if (Double.doubleToLongBits(getScale()) != Double.doubleToLongBits(other.getScale()))
			return false;
		if (Double.doubleToLongBits(getShape()) != Double.doubleToLongBits(other.getShape()))
			return false;
		return true;
	}
	
	public static ExtendedGammaDistribution fromMoments(double mean, double sd) {
		if(Double.isNaN(mean) || Double.isNaN(sd)) return null;
		return ExtendedGammaDistribution.fromShapeAndScale(
			FastMath.pow(mean/sd,2), //shape 
			FastMath.pow(sd,2)/mean //scale
		);
	}
	
	public static ExtendedGammaDistribution fromShapeAndScale(double shape, double scale) {
		if(Double.isNaN(shape) || Double.isNaN(scale)) return null;
		return new ExtendedGammaDistribution(shape,scale);
	}
	
	public static ExtendedGammaDistribution fromShapeAndRate(double shape, double rate) {
		if(Double.isNaN(shape) || Double.isNaN(rate)) return null;
		return new ExtendedGammaDistribution(shape,1.0/rate);
	}
	
}
