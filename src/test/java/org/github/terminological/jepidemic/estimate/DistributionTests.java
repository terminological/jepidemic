package org.github.terminological.jepidemic.estimate;

import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.Stream;

import org.apache.commons.math3.distribution.AbstractIntegerDistribution;
import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.github.terminological.jepidemic.distributions.BetaPrimeDistribution;
import org.github.terminological.jepidemic.distributions.EmpiricalDistribution;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;
import org.github.terminological.jepidemic.distributions.GrowthRateDistribution;
import org.github.terminological.jepidemic.distributions.NegativeBinomialDistribution;
import org.github.terminological.jepidemic.distributions.Summary;
import org.github.terminological.jepidemic.distributions.TransformerDistribution;
import org.github.terminological.jepidemic.distributions.TruncatedIntegerDistribution;
import org.github.terminological.jepidemic.distributions.TruncatedRealDistribution;
import org.github.terminological.jepidemic.distributions.WeightedMixtureIntegerDistribution;
import org.github.terminological.jepidemic.distributions.WeightedMixtureRealDistribution;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;


class DistributionTests {

	@BeforeEach
	void setUp() throws Exception {
	}

	@Test
	final void testDistributionsAreNullIfParametersNA() {
		
		assert(ExtendedGammaDistribution.fromShapeAndRate(Double.NaN, Double.NaN) == null);
		assert(ExtendedGammaDistribution.fromShapeAndScale(Double.NaN, Double.NaN) == null);
		assert(ExtendedGammaDistribution.fromMoments(Double.NaN, Double.NaN) == null);
		assert(GrowthRateDistribution.withForwardEstimateAndBackwardEstimate(Double.NaN, Double.NaN, 1) == null);
		assert(NegativeBinomialDistribution.fromMoments(Double.NaN, Double.NaN) == null);
		assert(NegativeBinomialDistribution.fromSizeAndProbability(Double.NaN, Double.NaN) == null);
		assert(BetaPrimeDistribution.wikipediaParameters(Double.NaN, Double.NaN) == null);
	}

	
	@Test
	final void testWeightedMixture() {
		WeightedMixtureRealDistribution mix = WeightedMixtureRealDistribution.empty();
		mix.add(new NormalDistribution(-1,3), 2.0);
		mix.add(new NormalDistribution(1,5), 3.0);
		EmpiricalDistribution emp = new EmpiricalDistribution(mix.sample(100000));
		
		System.out.println("mean: "+mix.getNumericalMean());
		System.out.println("mean: "+emp.getNumericalMean());
		
		System.out.println("var: "+mix.getNumericalVariance());
		System.out.println("var: "+emp.getNumericalVariance());
		
		System.out.println("q.0.95: "+mix.inverseCumulativeProbability(0.95));
		System.out.println("q.0.95: "+emp.inverseCumulativeProbability(0.95));
		
		System.out.println("cdf @ 5: "+mix.cumulativeProbability(5));
		System.out.println("cdf @ 5: "+emp.cumulativeProbability(5));
		
		System.out.println("pdf @ 5: "+mix.density(5));
		System.out.println("pdf @ 5: "+emp.density(5));
	}
	
	private abstract class PrintHelper implements BiFunction<String,Function<AbstractRealDistribution,Object>,String> {}
	private abstract class PrintHelper2 implements BiFunction<String,Function<AbstractIntegerDistribution,Object>,String> {}
	
	PrintHelper printHelper(AbstractRealDistribution norm,AbstractRealDistribution trunc) {
		return new PrintHelper() {
			@Override
			public String apply(String t, Function<AbstractRealDistribution, Object> u) {
				Object out1 = u.apply(norm);
				if (out1 instanceof String) return String.format("%s, Norm: %s, Trunc: %s",t, u.apply(norm),u.apply(trunc));
				return String.format("%s, Norm: %1.6f, Trunc: %1.6f",t, u.apply(norm),u.apply(trunc));
			}
		};
	}
	
	PrintHelper2 printHelper(AbstractIntegerDistribution norm,AbstractIntegerDistribution trunc) {
		return new PrintHelper2() {
			@Override
			public String apply(String t, Function<AbstractIntegerDistribution, Object> u) {
				Object out1 = u.apply(norm);
				if (out1 instanceof String) return String.format("%s, Norm: %s, Trunc: %s",t, u.apply(norm),u.apply(trunc));
				return String.format("%s, Norm: %1.6f, Trunc: %1.6f",t, u.apply(norm),u.apply(trunc));
			}
		};
	}
	
	
	@Test
	final void testTruncation() {
		
		Stream.of(
				new NormalDistribution(0,3),
				ExtendedGammaDistribution.fromMoments(3, 2)
				).forEach(norm -> {
					
					//TruncatedRealDistribution<NormalDistribution> trunc = TruncatedRealDistribution.of(norm, -4, 4);
					TruncatedRealDistribution trunc = TruncatedRealDistribution.of(norm);
					
					// assert(norm.getNumericalMean() == trunc.getNumericalMean());
					PrintHelper p = printHelper(norm,trunc);
					
					
					Stream.of(
							p.apply("Mean", d -> d.getNumericalMean()),
							p.apply("Variance", d -> d.getNumericalVariance()),
							p.apply("q.0.05", d -> d.inverseCumulativeProbability(0.05)),
							p.apply("q.0.95", d -> d.inverseCumulativeProbability(0.95))
					).forEach(System.out::println);;
		
				});
	}
	
	@Test
	final void testIntegerTruncation() {
		
		Stream.of(
				new PoissonDistribution(4),
				NegativeBinomialDistribution.fromMoments(3, 2)
				).forEach(norm -> {
					
					//TruncatedRealDistribution<NormalDistribution> trunc = TruncatedRealDistribution.of(norm, -4, 4);
					TruncatedIntegerDistribution trunc = TruncatedIntegerDistribution.of(norm);
					
					// assert(norm.getNumericalMean() == trunc.getNumericalMean());
					PrintHelper2 p = printHelper(norm,trunc);
					
					
					Stream.of(
							p.apply("Mean", d -> d.getNumericalMean()),
							p.apply("Variance", d -> d.getNumericalVariance()),
							p.apply("q.0.05", d -> (double) d.inverseCumulativeProbability(0.05)),
							p.apply("q.0.95", d -> (double) d.inverseCumulativeProbability(0.95))
					).forEach(System.out::println);;
		
				});
	}
	
	
	@Test
	final void testWeightedIntMixture() {
		WeightedMixtureIntegerDistribution mix = WeightedMixtureIntegerDistribution.empty();
		mix.add(new PoissonDistribution(4), 2.0);
		mix.add(new PoissonDistribution(5), 3.0);
		EnumeratedIntegerDistribution emp = new EnumeratedIntegerDistribution(mix.sample(100000));
		System.out.println("mean: "+mix.getNumericalMean());
		System.out.println("mean: "+emp.getNumericalMean());
		
		System.out.println("var: "+mix.getNumericalVariance());
		System.out.println("var: "+emp.getNumericalVariance());
		
		System.out.println("q.0.95: "+mix.inverseCumulativeProbability(0.95));
		System.out.println("q.0.95: "+emp.inverseCumulativeProbability(0.95));
		
		System.out.println("cdf @ 5: "+mix.cumulativeProbability(5));
		System.out.println("cdf @ 5: "+emp.cumulativeProbability(5));
		
		System.out.println("pdf @ 5: "+mix.probability(5));
		System.out.println("pdf @ 5: "+emp.probability(5));
	}
	
	@Test void testTransformerDistribution() {
		
//		
		
		TransformerDistribution<NormalDistribution> trans = TransformerDistribution.of(
				new NormalDistribution(0.5,0.25),
				x -> Math.exp(x),
				y -> Math.log(y),
				dy -> 1/dy
				);
		
		LogNormalDistribution lnorm = new LogNormalDistribution(0.5,0.25);
		
		
		PrintHelper p = printHelper(lnorm,trans);
		
		
		Stream.of(
				p.apply("Summary", d -> Summary.of(d).toString()),
				p.apply("Mean", d -> d.getNumericalMean()),
				p.apply("Variance", d -> d.getNumericalVariance()),
				p.apply("q.0.05", d -> d.inverseCumulativeProbability(0.05)),
				p.apply("q.0.95", d -> d.inverseCumulativeProbability(0.95)),
				p.apply("d.1", d -> d.density(1D))
		).forEach(System.out::println);
		
	}
	
	@Test void testBetaPrime() {
//		extraDistr::qbetapr(c(0.025,0.5,0.975), shape1 = 40, shape2 = 40)
//		[1] 0.6431356 1.0000000 1.5548820
		System.out.println(Summary.of(BetaPrimeDistribution.extraDistrParameters(3, 5, 1)));
		
		System.out.println("====== 3,5,1");
		BetaPrimeDistribution bp = BetaPrimeDistribution.extraDistrParameters(3, 5, 1);
		DoubleStream.of(0,0.1,0.2,0.3,0.4,0.5).map(x-> bp.density(x)).forEach(System.out::println);
		System.out.println("====== 3,5,1: CDF");
		DoubleStream.of(0,0.1,0.2,0.3,0.4,0.5).map(x-> bp.cumulativeProbability(x)).forEach(System.out::println);
		
		System.out.println("====== 3,5,1: MEAN: "+bp.getNumericalMean());
		System.out.println("====== 3,5,1: MEAN: "+bp.getNumericalVariance());
//		extraDistr::dbetapr(c(0,0.1,0.2,0.3,0.4,0.5), shape1 = 3, shape2 = 5)
		// [1] 0.0000000 0.4898327 0.9767858 1.1584705 1.1383741 1.0242341
//		> integrate(f = function(x) x*extraDistr::dbetapr(x, shape1 = 3, shape2 = 5, scale = 1),lower = -Inf, upper = Inf)
//		0.75 with absolute error < 8.3e-15
//		> tmp = integrate(f = function(x) x^2*extraDistr::dbetapr(x, shape1 = 3, shape2 = 5, scale = 1),lower = -Inf, upper = Inf)
//		> tmp$value
//		[1] 1
//		> tmp$value - 0.75^2
//		[1] 0.4375
//		> 
		
		System.out.println("====== 3,5,0.5");
		BetaPrimeDistribution bp2 = BetaPrimeDistribution.extraDistrParameters(3, 5, 0.5);
		DoubleStream.of(0,0.1,0.2,0.3,0.4,0.5).map(x-> bp2.density(x)).forEach(System.out::println);
		System.out.println("====== 3,5,0.5: CDF");
		DoubleStream.of(0,0.1,0.2,0.3,0.4,0.5).map(x-> bp2.cumulativeProbability(x)).forEach(System.out::println);
//		> extraDistr::pbetapr(c(0,0.1,0.2,0.3,0.4,0.5), shape1 = 3, shape2 = 5, scale = 0.5)
//		[1] 0.00000000 0.09577546 0.32077014 0.52465296 0.67266671 0.77343750
		System.out.println("====== 3,5,0.5: MEAN: "+bp2.getNumericalMean());
		System.out.println("====== 3,5,0.5: MEAN: "+bp2.getNumericalVariance());
//		> extraDistr::dbetapr(c(0,0.1,0.2,0.3,0.4,0.5), shape1 = 3, shape2 = 5, scale = 0.5)
//		[1] 0.0000000 1.9535715 2.2767481 1.7601997 1.2196051 0.8203125
		
	}
	
	@Test void testGrowthRateDistribution() {
		
		double Is_minus = 40;
		double Is_plus = 90;
		ExtendedGammaDistribution prior = ExtendedGammaDistribution.fromMoments(40, 5);
		
		System.out.println(prior.getShape());
		
		double tau = 5; 
		double alphaPlus = Is_plus+prior.getShape();
		double alphaMinus = Is_minus+prior.getShape();
		double betaP = prior.getRate() + 2*tau + 1;
		
		ExtendedGammaDistribution num = ExtendedGammaDistribution.fromShapeAndRate(alphaPlus, betaP);
		ExtendedGammaDistribution denom = ExtendedGammaDistribution.fromShapeAndRate(alphaMinus, betaP);
		
		System.out.println("Gammas numerator");
		System.out.println(Summary.of(num));
		System.out.println("Gammas denominator");
		System.out.println(Summary.of(denom));
		
		EmpiricalDistribution bpComp = EmpiricalDistribution.bySampling(num, denom, (x,y) -> x/y, 10000);
		BetaPrimeDistribution betaPrime = BetaPrimeDistribution.extraDistrParameters(alphaPlus, alphaMinus, 1);
		
		System.out.println("Ratio of gammas by sampling");
		System.out.println(Summary.of(bpComp));
		
		System.out.println("Ratio of gammas as beta prime");
		System.out.println(Summary.of(betaPrime));
		
		
		System.out.println("GR by sampling");
		EmpiricalDistribution grComp = EmpiricalDistribution.bySampling(num, denom, (x,y) -> 1/(2*tau)*Math.log(x/y), 10000);
		System.out.println(Summary.of(grComp));
		
		//System.out.println(Summary.of(BetaPrimeDistribution.extraDistrParameters(alphaPlus, alphaMinus, 1)).quantiles(new double[] {0.025,0.5,0.975}));
		
		System.out.println("GR by growth rate distribution");
		GrowthRateDistribution dist1 = GrowthRateDistribution.withForwardEstimateAndBackwardEstimate(alphaPlus, alphaMinus, (int) tau);
		System.out.println(Summary.of(dist1));
		
		// System.out.println(alphaPlus);
		// System.out.println(alphaMinus);
		
		System.out.println("GR by beta transformer");
		TransformerDistribution<BetaPrimeDistribution> grEstimate = TransformerDistribution.of(
				betaPrime,
				/*y=g(x)*/ x -> Math.log(x) / (2.0*tau),
				/*x=g^-1(y)*/ y -> Math.exp(2*tau*y),
				/*dg-1/dy =*/ y -> 2.0*tau*Math.exp(2*tau*y)
				);
		System.out.println(Summary.of(grEstimate));
		
		
		PrintHelper p = printHelper(dist1,grEstimate);
		
		
		Stream.of(
				p.apply("Summary", d -> Summary.of(d).toString()),
				p.apply("Mean", d -> d.getNumericalMean()),
				p.apply("Variance", d -> d.getNumericalVariance()),
				p.apply("q.0.05", d -> d.inverseCumulativeProbability(0.05)),
				p.apply("q.0.95", d -> d.inverseCumulativeProbability(0.95)),
				p.apply("d.1", d -> d.density(1D))
		).forEach(System.out::println);
	}
}
