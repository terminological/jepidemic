package org.github.terminological.jepidemic.estimate;

public class CoriBootstrappedTimeseries extends CoriTimeseries implements Cloneable {

	public CoriBootstrappedTimeseries(CoriEstimator coriEstimator) {
		super(coriEstimator);
	}

	public CoriBootstrappedTimeseries sample() {
		CoriBootstrappedTimeseries ts = new CoriBootstrappedTimeseries(this.getEstimator());
		for(CoriTimeseriesEntry tse : this) {
			((CoriBootstrappedTimeseriesEntry) tse).sampleInto(ts, 1);
		}
		return ts;
	}

	
	

}
