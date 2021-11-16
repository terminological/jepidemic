package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.github.terminological.jepidemic.distributions.Summary;

public class GammaParameters implements Summary {
	
	private double scale;
	private double shape;
	GammaDistribution g;
	
	public GammaParameters(double shape, double scale) {
		this.scale = scale; 
		this.shape = shape;
	}
	
	public GammaParameters(RandomGenerator rng,double shape, double scale) {
		this.scale = scale; 
		this.shape = shape;
	}
	
	public GammaMoments convert() {
		return new GammaMoments(
				getMean(),
				Math.sqrt(getVar()));
	}
	
	public GammaDistribution dist() {
		if (g==null) g=new GammaDistribution(getShape(), getScale());
		return g;
	}
	
//	Map<Double,Double> quantiles(double[] quantiles) {
//		Map<Double,Double> out = new HashMap<>();
//		Arrays.stream(quantiles).boxed().forEach(p -> out.put(p, dist().inverseCumulativeProbability(p)));
//		return out;
//	}
	
	
	
	public boolean isDefined() {
		return Double.isFinite(getScale()) && Double.isFinite(getShape()) &&
				getScale() >= 0 && getShape() >= 0;
	}
	
	public String toString() {return this.convert().toString();}
	
	public DatedGammaEstimate withDate(int tau, LocalDate date, double incidence,double incidenceInWindow) {
		return new DatedGammaEstimate(this, tau, date, incidence, incidenceInWindow);
	}
	
//	@Deprecated
//	public static GammaParameters welchSatterthwaiteCombination(List<? extends GammaParameters> distributions) {
//		//https://stats.stackexchange.com/questions/72479/generic-sum-of-gamma-random-variables
//		double sumShapeScale = 0;
//		double sumShapeScale2 = 0;
//		double sumShape = 0;
//		double N=0;
//		for (GammaParameters dist: distributions) {
//			if (dist.isDefined()) {
//				sumShape = sumShape+dist.getShape();
//				sumShapeScale = sumShapeScale+dist.getShape()*dist.getScale();
//				sumShapeScale2 = sumShapeScale2+dist.getShape()*dist.getScale()*dist.getScale();
//				N+=1;
//			}
//		}
//		
//		return new GammaParameters(
//			sumShapeScale*sumShapeScale/sumShapeScale2/N,
//			sumShapeScale2/sumShapeScale
//		); //.convert().wider(2).convert();
//	}

	public double getShape() {
		return shape;
	}

	public double getScale() {
		return scale;
	}
	
	public double getRate() {
		return 1.0/scale;
	}
	
//	public static GammaParameters normalAssumptionCombination(List<? extends GammaParameters> distributions) {
//		return GammaMoments.normalAssumptionEstimate(
//				distributions.stream().map(d -> d.convert()).collect(Collectors.toList())
//				).convert();
//	}
	
	public static Summary resamplingCombination(List<? extends GammaParameters> distributions, int sampleSize) {
		DescriptiveStatistics ds = new DescriptiveStatistics(); 
		distributions.parallelStream().flatMapToDouble(d -> Arrays.stream(d.dist().sample(sampleSize))).forEach(ds::addValue);
		return new Summary() {

			@Override
			public double getMean() {
				return ds.getMean();
			}

			@Override
			public double getSD() {
				return Math.sqrt(ds.getVariance());
			}

			@Override
			public double quantile(double q) {
				return ds.getPercentile(q*100);
			}
		};
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(scale);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(shape);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof GammaParameters))
			return false;
		GammaParameters other = (GammaParameters) obj;
		if (Double.doubleToLongBits(scale) != Double.doubleToLongBits(other.scale))
			return false;
		if (Double.doubleToLongBits(shape) != Double.doubleToLongBits(other.shape))
			return false;
		return true;
	}

	public static <X extends DatedGammaEstimate> Summary mixtureDistribution(List<X> cere, Function<X,Double> weights) {
		return new WeightedGammaMixture(cere, weights);
	}

	public static Summary mixtureDistribution(List<? extends DatedGammaEstimate> cere) {
		return new WeightedGammaMixture(cere);
	}

	public double getMean() {
		return getShape()*getScale();
	}

	public double getVar() {
		return getShape()*getScale()*getScale();
	}

	public static Summary mixtureApproximation(List<? extends DatedGammaEstimate> cere) {
		return new WeightedGammaMixture(cere).gammaApproximation().convert();
	}

	@Override
	public double getSD() {
		return Math.sqrt(getVar());
	}

	@Override
	public double quantile(double q) {
		if (!this.isDefined()) return Double.NaN;
		return dist().inverseCumulativeProbability(q);
	}
	
	public static GammaParameters fromShapeAndScale(double shape, double scale) {
		return new GammaParameters(shape,scale);
	}
	
	public static GammaParameters fromShapeAndRate(double shape, double rate) {
		return new GammaParameters(shape,1.0/rate);
	}
}