package org.github.terminological.jepidemic.gamma;

public interface StatSummary {

	public double getMean();
	public double getSD();
	public double quantile(double q);
	
//	public default CoriEstimationSummaryEntry withDate(LocalDate date, int window, double incidence) {
//		return new CoriEstimationSummaryEntry(this, date, window, incidence);
//	}
}
