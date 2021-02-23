package org.github.terminological.jepidemic;

import org.github.terminological.jepidemic.estimate.CoriEstimationSummaryEntry;
import org.github.terminological.jepidemic.gamma.GammaMoments;

import uk.co.terminological.rjava.RName;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;

public interface RtTimeseriesEntry extends TimeseriesEntry {

	@RName("window") RInteger window();
	@RName("endDate") RDate endDate();
	@RName("startDate") RDate startDate();
	@RName("Rt.Mean") RNumeric mean();
	@RName("Rt.SD") RNumeric sd();
	
	
	
	public default CoriEstimationSummaryEntry assumeGamma() {
		if (this instanceof CoriEstimationSummaryEntry) return (CoriEstimationSummaryEntry) this;
		return new GammaMoments(this.mean().get(), this.sd().get()).convert()
		.withDate(-1, date().get(), incidence().get()).toStatSummary();
	}
	
	
}
