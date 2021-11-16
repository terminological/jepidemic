package org.github.terminological.jepidemic.estimate;

import org.github.terminological.jepidemic.InfectivityProfile;

public class GrowthRateUtils {

	GammaMoments gammaSI;
	InfectivityProfile profile;
	
	public GrowthRateUtils(GammaMoments serialInterval) {
		this.gammaSI = serialInterval;
	}
	
	public GrowthRateUtils(InfectivityProfile profile) {
		this.gammaSI = gammaSerialInterval(profile);
		this.profile = profile;
	}
	
//	private static GammaMoments gammaSerialInterval(List<InfectivityProfile> profiles) {
//		double meanOfMean = profiles.stream().map(p -> p.getNumericalMean()).collect(Collectors.averagingDouble(d -> d));
//		double meanOfSd = profiles.stream().map(p -> Math.sqrt(p.getNumericalVariance())).collect(Collectors.averagingDouble(d -> d));
//		return new GammaMoments(meanOfMean, meanOfSd);
//	}
	
	private static GammaMoments gammaSerialInterval(InfectivityProfile profile) {
		double meanOfMean = profile.getNumericalMean();		
		double meanOfSd = Math.sqrt(profile.getNumericalVariance());
		return new GammaMoments(meanOfMean, meanOfSd);
	}
	
	public double growthRateToR(double r) {
		return doGrowthToR(r,gammaSI);
	}
	
	
	public double RToGrowth(double R) {
		return doRToGrowth(R,gammaSI);
	}
	
	// Convert Growth Rates to Reproduction numbers.
	//
	// @description See [here](https://www.medrxiv.org/content/10.1101/2020.01.30.20019877v3.full.pdf) 
	// for justification.
	// @param r Numeric, rate of growth estimates
	// @param gamma_mean Numeric, mean of the gamma distribution
	// @param gamma_sd Numeric, standard deviation of the gamma distribution
	//
	// @return Numeric vector of reproduction number estimates
	// @export
	//
	// @examples
	// 
	// growth_to_R(0.2, 4, 1)
	public static double doGrowthToR(double r, GammaMoments gammaSI) {
	  double k = Math.pow((gammaSI.getSD() / gammaSI.getMean()),2);
	  double R = Math.pow(1 + k * r * gammaSI.getMean(),(1/k));
	  return R;
	}
	  
	// Convert Reproduction Numbers to Growth Rates
	//
	// @description See [here](https://www.medrxiv.org/content/10.1101/2020.01.30.20019877v3.full.pdf) 
	// for justification.
	// @param R Numeric, Reproduction number estimates
	// @inheritParams growth_to_R
	// @return Numeric vector of reproduction number estimates
	// @export
	//
	// @examples
	// 
	// R_to_growth(2.18, 4, 1)  
	public static double doRToGrowth(double R, GammaMoments gammaSI) {
	  double k = Math.pow((gammaSI.getSD()/ gammaSI.getMean()),2);
	  double r = (Math.pow(R,k) - 1) / (k * gammaSI.getMean());
	  return r;
	}
	
}
