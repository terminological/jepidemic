package org.github.terminological.jepidemic.estimate;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertTrue;

import static uk.co.terminological.rjava.RFunctions.*;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import org.github.terminological.jepidemic.IncompleteTimeseriesException;
import org.github.terminological.jepidemic.gamma.GammaParameters;
import org.junit.jupiter.api.Test;

import uk.co.terminological.rjava.IncompatibleTypeException;
import uk.co.terminological.rjava.RFunctions;
import uk.co.terminological.rjava.ZeroDimensionalArrayException;
import uk.co.terminological.rjava.types.RDataframe;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNamedList;
import uk.co.terminological.rjava.types.RNumeric;
import uk.co.terminological.rjava.types.RNumericArray;
import uk.co.terminological.rjava.types.RNumericVector;
import uk.co.terminological.rjava.types.RObject;


public class TestEstimator {

	static RDataframe getFlu2009() throws IOException {
		InputStream is = TestEstimator.class.getResourceAsStream("/flu2009.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RDataframe getFlu2009SerialInterval() throws IOException {
		InputStream is = TestEstimator.class.getResourceAsStream("/flu2009SerialInterval.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RDataframe getFlu2009Result() throws IOException {
		InputStream is = TestEstimator.class.getResourceAsStream("/flu2009EpiEstim.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RDataframe getPheApiOutput() throws IOException {
		InputStream is = TestEstimator.class.getResourceAsStream("/pheApi.ser");
		if(is==null) throw new IOException("Could not locate input file");
		return RObject.readRDS(RDataframe.class, is);
	}
	
	static RNamedList getCovidSerialInterval() throws IOException {
		InputStream is = TestEstimator.class.getResourceAsStream("/epiestimCovidConfig.ser");
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
	
	private static CoriEstimator basicEst() throws IOException {
		CoriEstimator est = new CoriEstimator(5,5,14);
		est.atStartOfTimeseries();
		est.withInfectivityProfile(getFlu2009SerialInterval().pull("prob",RNumericVector.class));
		return est;
	}
	
	@Test
	public void testCoriEstimator1() throws IOException, IncompleteTimeseriesException {
		CoriEstimator est = basicEst();
		
		for (int i=1;i<=3;i++) {
			String iLab = null;
			switch (i) {
			case 1: 
				est.withAdaptivePrior(1.25);
				iLab = "adaptive";
				break;
			case 2: 
				est.withDefaultPrior();
				iLab = "default";
				break;
			case 3: 
				est.withFixedPrior(2, 2);
				iLab = "fixed";
				break;
			}
		
			for (int j=1; j<=4; j++) {
				String jLab = null;
				switch(j) {
				case 1:
					est.selectSpecificWindow(7);
					jLab = "specific";
					break;
				case 2:
					est.selectAdaptiveWindow(24,3);
					jLab = "adaptive";
					break;
				case 3:
					est.selectMinimumUncertainty(0.5,3);
					jLab = "minimum";
					break;
				case 4:
					est.selectMixtureCombination();
					jLab = "mixture";
					break;
				}
			
				
				RDataframe out = est.estimateRtSingle(getFlu2009(), "dates", "I");
				RNumericVector res = out.pull("Rt.Mean",RNumericVector.class);
				System.out.println(iLab+" "+jLab+" "+res.rCode());
			}
		}
		
		
	}
	
//	@Test
//	public void testCoriEstimator() throws IOException, IncompleteTimeseriesException {
//		CoriEstimator est = basicEst();
//		//est.selectSpecificWindow(7);
//		est.selectAdaptiveWindow(24,3);
//		//est.selectMinimumUncertainty(0.5,3);
//		est.selectMixtureCombination();
//		est.atStartOfTimeseries();
//		est.collectMixtureQuantiles();
//		//est.collectResampledQuantiles(1000);
//		//est.collectNormalApproximateQuantiles();
//		System.out.println(getFlu2009());
//		RDataframe out = est.estimateRtSingle(getFlu2009(), "dates", "I");
//		System.out.println(out.asCsv());
//	}
	
	@Test
	public void testMultipleSI() throws IOException, ZeroDimensionalArrayException {
		RNamedList nl = getCovidSerialInterval();
		RNumericArray tmp = nl.getAs("si_sample", RNumericArray.class);
		CoriEstimator est = new CoriEstimator(5,5,28).legacySupport(true);
		tmp.get().forEach(v -> est.withInfectivityProfile(v.getVector()));
		est.withAdaptivePrior(1.1);
		est.selectAdaptiveWindow(100, 3);
		est.collectMixtureApproximation();
		RDataframe df = getPheApiOutput();
		df.mutate("value", RInteger.class, v -> RFunctions.asNumeric(v));
		RDataframe out = est.estimateRt(df, "date", "value");
		System.out.println(out.asCsv());
	}
	
	@Test
	public void testFluResult() throws IOException, IncompleteTimeseriesException {
		
		CoriEstimator est = CoriEstimator.defaultEpiEstim(getFlu2009SerialInterval().pull("prob",RNumericVector.class),5,5,7);
		RDataframe out = est.estimateRtSingle(getFlu2009(), "dates", "I");
		RDataframe out2 = out.filter("Mean(R)", n -> Double.isFinite(n.get()));
		
		RNumericVector ref = out2.pull("Mean(R)", RNumericVector.class);
		RNumericVector ref2 = getFlu2009Result().pull("Mean(R)", RNumericVector.class);
		
		assertTrue(all((v1,v2) -> precisionEquals(v1,v2,0.00001), ref, ref2));
		
		System.out.println("Matches R output to 5 dp:");
		System.out.println(ref);
		System.out.println(ref2);
	}
}
