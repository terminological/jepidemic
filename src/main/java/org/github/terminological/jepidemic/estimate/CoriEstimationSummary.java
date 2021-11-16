package org.github.terminological.jepidemic.estimate;

import java.util.ArrayList;

import org.github.terminological.jepidemic.Timeseries;

public class CoriEstimationSummary extends Timeseries<CoriEstimationSummaryEntry> {

	private CoriEstimator est;

	public CoriEstimationSummary(CoriEstimator est) {
		super((d1,d2) -> d1.date.compareTo(d2.date));
		this.est = est;
	}
	
	public CoriEstimationSummary(CoriEstimationResult unfiltered) {
		super((d1,d2) -> d1.date.compareTo(d2.date));
		est = unfiltered.estimator;
		unfiltered.forEach(cere -> 
			this.add(new CoriEstimationSummaryEntry(
					new ArrayList<>(cere.posteriors()), this))
				);
	}

	public CoriEstimator getEstimator() {
		return est;
	}

	public void merge(CoriEstimationSummaryEntry t) {
		if (this.contains(t)) {
			this.forEach(t2 -> {
				if (t2.getDate().equals(t.getDate())) {
					t2.merge(t);
				}
			});
		} else {
			this.add(t);
		}
	}

	
	
}
