package org.github.terminological.jepidemic.estimate;

import org.github.terminological.jepidemic.Timeseries;

public class CoriEstimationSummary extends Timeseries<CoriEstimationSummaryEntry> {

	public CoriEstimationSummary() {
		super((d1,d2) -> d1.date.compareTo(d2.date));
	}

}
