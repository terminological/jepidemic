package org.github.terminological.jepidemic.estimate;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.util.FastMath;
import org.github.terminological.jepidemic.distributions.Summary;

public class GammaMoments implements Summary {
	private double mean;
	private double sd;
	
	public GammaMoments(double mean, double sd) {this.mean = mean; this.sd = sd;}
	public GammaParameters convert() {
		return new GammaParameters(
				FastMath.pow(mean/sd,2), //shape 
				FastMath.pow(sd,2)/mean); //scale
	}
	public GammaMoments wider(double factor) {
		return new GammaMoments(mean, sd*factor);
	}
	public double getMean() {return mean;}
	public double getSD() {return sd;}
	
	public String toString() {return mean+"\u00B1"+sd;}
	
	GammaDistribution g;
	
	public GammaDistribution dist() {
		if (g==null) g=new GammaDistribution(this.convert().getShape(), this.convert().getScale());
		return g;
	}
	
	public double quantile(double p) {
		if (!this.convert().isDefined()) return Double.NaN;
		return dist().inverseCumulativeProbability(p);
	}
	
//	public static GammaMoments normalAssumptionEstimate(List<? extends GammaMoments> distributions) {
//		double meanOfMean = distributions.stream().mapToDouble(d -> d.getMean()).summaryStatistics().getAverage();
//		double meanOfSD = distributions.stream().mapToDouble(d -> d.getSD()).summaryStatistics().getAverage();
//		return new GammaMoments(meanOfMean, meanOfSD);	
//	}
	
}