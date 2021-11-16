package org.github.terminological.jepidemic.estimate;

import java.util.Arrays;

import org.junit.jupiter.api.Test;

public class TestMixture {

	@Test
	public void smokeTest() {
		
		WeightedGammaMixture mix = new WeightedGammaMixture(Arrays.asList(
//				new GammaParameters(2,3),
//				new GammaParameters(5,2),
//				new GammaParameters(7,4)
				new GammaMoments(4, 2).convert(),
				new GammaMoments(3, 5).convert(),
				new GammaMoments(7, 1).convert()
				));
		
		System.out.println("=PDF=");
		Arrays.stream(new double[] {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0}).map(m -> mix.density(m)).forEach(System.out::println);
		System.out.println("=CDF=");
		Arrays.stream(new double[] {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0}).map(m -> mix.cumulativeProbability(m)).forEach(System.out::println);
		
		System.out.println("=MOMENTS=");
		System.out.println(mix.getNumericalMean());
		System.out.println(mix.getNumericalVariance());
		
		System.out.println("=QUANTILES=");
		System.out.println(mix.quantile(0.025));
		System.out.println(mix.quantile(0.25));
		System.out.println(mix.quantile(0.5));
		System.out.println(mix.quantile(0.75));
		System.out.println(mix.quantile(0.975));
		
	}
	
//	@Test
//	public void testWelchScatterbachEstimate() {
//		GammaParameters tmp = GammaParameters.welchSatterthwaiteCombination(Arrays.asList(
//				new GammaParameters(3,1),
//				new GammaParameters(4,2),
//				new GammaParameters(5,1)
//		));
//		System.out.println(tmp.getShape()+" "+tmp.getScale());
//		System.out.println(tmp.toString());		
//	}
}
