package org.github.terminological.jepidemic.distributions;

import java.util.function.Function;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.RombergIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.exception.OutOfRangeException;

public class TransformerDistribution<X extends AbstractRealDistribution> extends AbstractRealDistribution {

	Function<Double,Double> g_x;
	// Function<Double,Double> dg_x_dx;
	Function<Double,Double> inv_g_y;
	Function<Double,Double> dinv_g_y_dy;
	UnivariateIntegrator ti;
	X big_x;
	Double mean;
	Double var;
	Double lb;
	Double ub;
	double meanPrecision;
	double sdPrecision;
	static double MEAN_PRECISION = 0.00001;
	static double SD_PRECISION = 0.0000001;
	
	// TODO: this only set up for strictly increasing functions of real distributions
	// Needs adjustment for decreaing functions.
	// https://www.statlect.com/fundamentals-of-probability/functions-of-random-variables-and-their-distribution
	/**
	 * @param distribution - the pre transformation distribution X. when we want the distribution of Y
	 * @param transform - Y = g(X) what is the transformation of X that gives us Y? e.g. Y=e^X
	 * @param inverseTransform - Y = g^{-1}(Y) what is the inverse? X=log(Y)
	 * @param inverseTransformDifferential - dg^{-1}/dy the differential of the inverse with respect to Y? e.g. DY = 1/Y
	 * @param meanPrecision - the quantile limit of the integration to estimate the mean
	 * @param sdPrecision - the quantile limit of the integration to estimate the sd
	 */
	public TransformerDistribution(
			X distribution, 
			Function<Double, Double> transform, 
			Function<Double, Double> inverseTransform,
			Function<Double, Double> inverseTransformDifferential,
			double meanPrecision,
			double sdPrecision
		) {
		super(Random.RNG);
		this.g_x = transform;
		this.inv_g_y = inverseTransform;
		this.dinv_g_y_dy = inverseTransformDifferential;
		this.big_x = distribution;
		this.ub = g_x.apply(big_x.getSupportUpperBound());
		this.lb = g_x.apply(big_x.getSupportLowerBound());
		this.ti = new RombergIntegrator();
		this.meanPrecision = meanPrecision;
		this.sdPrecision = sdPrecision;
	}
	@Override
	public double density(double y) {
		if (y < getSupportLowerBound()) return 0;
		if (y > getSupportUpperBound()) return 0;
		double tmp1 = big_x.density(inv_g_y.apply(y)); 
		double tmp2 = dinv_g_y_dy.apply(y);
		return tmp1 * tmp2;
	}
	@Override
	public double cumulativeProbability(double y) {
		if (y < getSupportLowerBound()) return 0;
		if (y > getSupportUpperBound()) return 1;
		return big_x.cumulativeProbability(inv_g_y.apply(y));
	}
	
	public X getBareDistribution() {
		return big_x;
	}
	
	@Override
	public double inverseCumulativeProbability(double p) throws OutOfRangeException {
		double x = big_x.inverseCumulativeProbability(p);
		double y = g_x.apply(x);
		return y;
	}
	@Override
	public double sample() {
		return g_x.apply(big_x.sample());
	}
	@Override
	public double getNumericalMean() {
		if (mean == null) {
			UnivariateFunction f = new UnivariateFunction() {
				@Override
				public double value(double x) {
					return x * TransformerDistribution.this.density(x) / (1-2*meanPrecision);
				}};
			double llimit = inverseCumulativeProbability(meanPrecision);
			double ulimit = inverseCumulativeProbability(1-meanPrecision); 
			mean = this.ti.integrate(Integer.MAX_VALUE, f, llimit, ulimit);
		}
		return(mean);
	}
	@Override
	public double getNumericalVariance() {
		if (var == null) {
			UnivariateFunction f = new UnivariateFunction() {
				@Override
				public double value(double x) {
					return Math.pow(x,2) * TransformerDistribution.this.density(x) / (1-2*sdPrecision);
				}};
			double llimit = inverseCumulativeProbability(sdPrecision);
			double ulimit = inverseCumulativeProbability(1-sdPrecision);
			var = this.ti.integrate(Integer.MAX_VALUE, f, llimit, ulimit) - Math.pow(getNumericalMean(),2);
		}
		return(var);
	}
	@Override
	public double getSupportLowerBound() {
		return lb;
	}
	@Override
	public double getSupportUpperBound() {
		return ub;
	}
	@SuppressWarnings("deprecation")
	@Override
	public boolean isSupportLowerBoundInclusive() {
		return big_x.isSupportLowerBoundInclusive();
	}
	@SuppressWarnings("deprecation")
	@Override
	public boolean isSupportUpperBoundInclusive() {
		return big_x.isSupportUpperBoundInclusive();
	}
	@Override
	public boolean isSupportConnected() {
		return big_x.isSupportConnected();
	}
	
	public static <Y extends AbstractRealDistribution> TransformerDistribution<Y> of(
			Y distribution, 
			Function<Double, Double> transform, 
			Function<Double, Double> inverseTransform,
			Function<Double, Double> inverseTransformDifferential 
		) {
		return new TransformerDistribution<Y>(
				distribution, 
				transform, 
				inverseTransform,
				inverseTransformDifferential,
				MEAN_PRECISION,
				SD_PRECISION
			);
	}
	
	public static <Y extends AbstractRealDistribution> TransformerDistribution<Y> of(
			Y distribution, 
			Function<Double, Double> transform, 
			Function<Double, Double> inverseTransform,
			Function<Double, Double> inverseTransformDifferential,
			double meanPrecision,
			double sdPrecision
		) {
		return new TransformerDistribution<Y>(
				distribution, 
				transform, 
				inverseTransform,
				inverseTransformDifferential,
				meanPrecision,
				sdPrecision
			);
	}
	
}
