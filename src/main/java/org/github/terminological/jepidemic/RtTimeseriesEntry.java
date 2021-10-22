package org.github.terminological.jepidemic;

import java.util.List;

import org.github.terminological.jepidemic.estimate.CoriEstimationSummary;
import org.github.terminological.jepidemic.estimate.CoriEstimationSummaryEntry;
import org.github.terminological.jepidemic.gamma.GammaMoments;

import uk.co.terminological.rjava.RName;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;

public interface RtTimeseriesEntry extends TimeseriesEntry.Incidence {

	@RName("window") RInteger window();
	@RName("endDate") RDate endDate();
	@RName("startDate") RDate startDate();
	@RName("Rt.Mean") RNumeric mean();
	@RName("Rt.SD") RNumeric sd();
	
	
	
	public default CoriEstimationSummaryEntry assumeGamma(CoriEstimationSummary summ, int profileId) {
		if (this instanceof CoriEstimationSummaryEntry) return (CoriEstimationSummaryEntry) this;
		return 
			new	CoriEstimationSummaryEntry(
				List.of(
						new GammaMoments(this.mean().get(), this.sd().get()).convert()
						.withDate(-1, date().get(), incidence().get(),profileId)
						),
				summ);
	}
	
	
}
