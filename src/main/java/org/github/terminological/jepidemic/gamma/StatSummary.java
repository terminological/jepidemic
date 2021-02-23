package org.github.terminological.jepidemic.gamma;

import java.time.LocalDate;

import org.github.terminological.jepidemic.estimate.CoriEstimationSummaryEntry;

public interface StatSummary {

	public double getMean();
	public double getSD();
	public double quantile(double q);
	
	public default CoriEstimationSummaryEntry withDate(LocalDate date, int window, double incidence) {
		return new CoriEstimationSummaryEntry(this, date, window, incidence);
	}
}
