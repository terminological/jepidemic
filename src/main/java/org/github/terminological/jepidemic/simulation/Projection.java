package org.github.terminological.jepidemic.simulation;

import java.time.LocalDate;
import java.util.function.Function;
import java.util.stream.Stream;

import org.github.terminological.jepidemic.InfectivityProfile;
import org.github.terminological.jepidemic.RtTimeseriesEntry;
import org.github.terminological.jepidemic.estimate.CoriEstimationSummary;
import org.github.terminological.jepidemic.estimate.CoriEstimationSummaryEntry;
import org.github.terminological.jepidemic.gamma.GammaParameters;

import uk.co.terminological.rjava.RClass;
import uk.co.terminological.rjava.RMethod;
import uk.co.terminological.rjava.UnconvertableTypeException;
import uk.co.terminological.rjava.types.RDataframe;

// @RClass
public class Projection {

	CoriEstimationSummary series;
	Function<LocalDate,GammaParameters> rtGenerator;
	
	// @RMethod
	public Projection(RDataframe dataframe) throws UnconvertableTypeException {
		this.series = new CoriEstimationSummary();
		
		dataframe.attach(RtTimeseriesEntry.class)
			.streamCoerce()
			.map(s -> s.assumeGamma())
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
