package org.github.terminological.jepidemic.simulation;

import java.time.LocalDate;
import java.util.function.Function;

import org.github.terminological.jepidemic.RtTimeseriesEntry;
import org.github.terminological.jepidemic.estimate.CoriEstimationSummary;
import org.github.terminological.jepidemic.estimate.CoriEstimator;
import org.github.terminological.jepidemic.gamma.GammaParameters;

import uk.co.terminological.rjava.UnconvertableTypeException;
import uk.co.terminological.rjava.types.RDataframe;

// @RClass
public class Projection {

	CoriEstimationSummary series;
	Function<LocalDate,GammaParameters> rtGenerator;
	
	// @RMethod
	public Projection(RDataframe dataframe) throws UnconvertableTypeException {
		this.series = new CoriEstimationSummary((CoriEstimator) null);
		
		dataframe.attach(RtTimeseriesEntry.class)
			.streamCoerce()
			.map(s -> s.assumeGamma(series,0))
			.forEach(series::add);
		
	}
	
	public Projection(CoriEstimationSummary series) {
		this.series = series;
	}
	
	
//	private Double nextIncidence(CoriEstimationSummaryEntry entry, InfectivityProfile infectivityProfile) {
//		int size = infectivityProfile.length();
//		Stream<CoriEstimationSummaryEntry> s = entry.laggingWindowInclusive(size).map(x -> (CoriEstimationSummaryEntry) x);
//		
//	}
	
}
