package org.github.terminological.jepidemic.estimate;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;

import java.io.IOException;
import java.io.InputStream;

import org.github.terminological.jepidemic.IncompleteTimeseriesException;
import org.github.terminological.jepidemic.growth.GrowthRateEstimator;
import org.junit.jupiter.api.Test;

import uk.co.terminological.rjava.IncompatibleTypeException;
import uk.co.terminological.rjava.ZeroDimensionalArrayException;
import uk.co.terminological.rjava.types.RDataframe;
import uk.co.terminological.rjava.types.RNamedList;
import uk.co.terminological.rjava.types.RNumericArray;
import uk.co.terminological.rjava.types.RNumericVector;
import uk.co.terminological.rjava.types.RObject;


public class TestGrowthRateEstimator {

	static RDataframe getFlu2009() throws IOException {
		InputStream is = TestGrowthRateEstimator.class.getResourceAsStream("/flu2009.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RDataframe getFlu2009SerialInterval() throws IOException {
		InputStream is = TestGrowthRateEstimator.class.getResourceAsStream("/flu2009SerialInterval.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RDataframe getFlu2009Result() throws IOException {
		InputStream is = TestGrowthRateEstimator.class.getResourceAsStream("/flu2009EpiEstim.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RDataframe getPheApiOutput() throws IOException {
		InputStream is = TestGrowthRateEstimator.class.getResourceAsStream("/pheApi.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RNamedList getCovidSerialInterval() throws IOException {
		InputStream is = TestGrowthRateEstimator.class.getResourceAsStream("/epiestimCovidConfig.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RNamedList.class, is);
	}
	
	@Test
	public void testTestData() throws IOException, IncompatibleTypeException, ZeroDimensionalArrayException {
		assertDoesNotThrow(() -> {
			System.out.println("===== FLU2009 incidence =====");
			System.out.println(getFlu2009());
			System.out.println("\n===== FLU2009 serial interval =====");
			System.out.println(getFlu2009SerialInterval());
			System.out.println("\n===== FLU2009 Epiestim output =====");
			System.out.println(getFlu2009Result());
			System.out.println("\n===== PHE API standard output =====");
			System.out.println(getPheApiOutput());
			System.out.println("\n===== EpiEstim complex config =====");
			//System.out.println(getCovidSerialInterval().getAs("si_sample",RNumericArray.class));
			System.out.println(getCovidSerialInterval());//.getAs("si_sample",RNumericArray.class));
			//System.out.println(Arrays.toString(getCovidSerialInterval().getAs("si_sample",RNumericArray.class).rPrimitive()));
		});
	}
	
	private static GrowthRateEstimator basicEst() throws IOException {
		GrowthRateEstimator est = GrowthRateEstimator.defaultGrowthRateEstimator(1, 14);
		
		return est;
	}
	
	@Test
	public void testGrowthRateEstimator() throws IOException, IncompleteTimeseriesException {
		GrowthRateEstimator est = basicEst();
		
		for (int i=1;i<=4;i++) {
			String iLab = null;
			switch (i) {
			case 1: 
				est.priorIncidenceFromScaledPreviousPosterior();
				iLab = "auto scaled posterior";
				break;
			case 2: 
				est.priorIncidenceFromPreviousPosterior(1.5);
				iLab = "fixed scaled posterior";
				break;
			case 3: 
				est.priorIncidenceFromScaledPosteriorAndGrowthRate();
				iLab = "auto scaled mechanistic";
				break;
			case 4: 
				est.priorIncidenceFromPosteriorAndGrowthRate(1.5);
				iLab = "fixed scaled mechanistic";
				break;
			}
		
			for (int j=1; j<=1; j++) {
				String jLab = null;
				switch(j) {
				case 1:
					est.combineEstimatesWithWeightedMixture();
					jLab = "weighted mixture";
					break;
				}
				
				
				for (int k=1; k<=3; k++) {
					String kLab = null;
					switch(k) {
					case 1:
						est.useAllPosteriorEstimates();
						kLab = "all estimates";
						break;
					case 2:
						est.usePosteriorEstimatesFromOneWindow(7);
						kLab = "fixed window";
						jLab = "collect first";
						break;
					case 3:
						est.usePosteriorEstimatesWithEnoughData(1);
						kLab = "minimum data";
						break;
					}
					
				
					RDataframe tmp = getFlu2009();
					RDataframe out = est.estimateGrowthRateSingle(tmp, "dates", "I");
					RNumericVector res = out.pull("Growth.Quantile.0.5.value",RNumericVector.class);
					System.out.println(iLab+" "+jLab+" "+kLab+" "+res.rCode());
				}
			}
		}
		
		
	}
	
	@Test
	public void testSingleCovidTs() throws IOException, IncompleteTimeseriesException, IncompatibleTypeException, ZeroDimensionalArrayException {
		GrowthRateEstimator est = GrowthRateEstimator.defaultGrowthRateEstimator(5, 14);
		est.withInfectivityProfile(
				getCovidSerialInterval().getAs("si_sample", RNumericArray.class).get(0).getVector(),true);
//		est.withInfectivityProfileMatrix(
//				getCovidSerialInterval().getAs("si_sample", RNumericArray.class));
		RDataframe out = est.estimateGrowthRateSingle(
				getPheApiOutput()
					.filter("name", n->n.get().equals("England"))
					.filter("statistic", n->n.get().equals("case"))
						, "date", "Imputed.value");
		System.out.println(out.asCsv());
	}
	
//	@Test
//	public void testMultipleSI() throws IOException, ZeroDimensionalArrayException {
//		RNamedList nl = getCovidSerialInterval();
//		RNumericArray tmp = nl.getAs("si_sample", RNumericArray.class);
//		CoriEstimator est = new CoriEstimator(5,5,28).legacySupport(true);
//		est.withInfectivityProfileMatrix(tmp);
//		est.withAdaptivePrior(1.1);
//		est.selectAdaptiveWindow(100, 3);
//		est.collectMixtureApproximation();
//		RDataframe df = getPheApiOutput();
//		df.mutate("value", RInteger.class, v -> RFunctions.asNumeric(v));
//		RDataframe out = est.estimateRt(df, "date", "value");
//		System.out.println(out.asCsv());
//	}
//	
//	@Test
//	public void testDetailedSI() throws IOException, ZeroDimensionalArrayException {
//		RNamedList nl = getCovidSerialInterval();
//		RNumericArray tmp = nl.getAs("si_sample", RNumericArray.class);
//		CoriEstimator est = new CoriEstimator(5,5,28).legacySupport(true);
//		est.withInfectivityProfileMatrix(tmp);
//		est.withAdaptivePrior(1.1);
//		est.selectAdaptiveWindow(100, 3);
//		est.detailedOutput();
//		RDataframe df = getFlu2009();
//		RDataframe out = est.estimateRt(df, "dates", "I");
//		System.out.println(out.asCsv());
//	}
//	
//	@Test
//	public void testBootsSI() throws IOException, ZeroDimensionalArrayException {
//		RNamedList nl = getCovidSerialInterval();
//		RNumericArray tmp = nl.getAs("si_sample", RNumericArray.class);
//		CoriEstimator est = new CoriEstimator(5,5,28).legacySupport(true);
//		est.withInfectivityProfileMatrix(tmp);
//		est.withAdaptivePrior(1.1);
//		est.selectAdaptiveWindow(100, 3);
//		est.collectMixtureApproximation();
//		RDataframe df = getPheApiOutput();
//		df.mutate("value", RInteger.class, v -> RFunctions.asNumeric(v));
//		RDataframe out = est.estimateRtFromRates(df, "date", "RollMean.value", 10);
//		System.out.println(out.asCsv());
//	}
//	
//	@Test
//	public void testFluResult() throws IOException, IncompleteTimeseriesException {
//		
//		CoriEstimator est = CoriEstimator.defaultEpiEstim(getFlu2009SerialInterval().pull("prob",RNumericVector.class),5,5,7);
//		RDataframe out = est.estimateRtSingle(getFlu2009(), "dates", "I");
//		RDataframe out2 = out.filter("Mean(R)", n -> Double.isFinite(n.get()));
//		
//		RNumericVector ref = out2.pull("Mean(R)", RNumericVector.class);
//		RNumericVector ref2 = getFlu2009Result().pull("Mean(R)", RNumericVector.class);
//		
//		assertTrue(all((v1,v2) -> precisionEquals(v1,v2,0.00001), ref, ref2));
//		
//		System.out.println("Matches R output to 5 dp:");
//		System.out.println(ref);
//		System.out.println(ref2);
//	}
	
	
}
