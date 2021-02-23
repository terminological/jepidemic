package org.github.terminological.jepidemic.estimate;

import java.util.Optional;
import java.util.TreeSet;

import org.github.terminological.jepidemic.IncompleteTimeseriesException;
import org.github.terminological.jepidemic.Timeseries;

public class CoriTimeseries extends Timeseries<CoriTimeseriesEntry> implements Cloneable {

	private CoriEstimator estimator;
	private double[] infectivityProfile;
	private int maxTau;
	//double r;
	
	public CoriTimeseries(CoriEstimator estimator) { //, double initialGrowth) {
		super((c1,c2) -> c1.compareTo(c2));
		this.estimator= estimator;
		this.maxTau = estimator.getMaxTau();
		//this.r = initialGrowth;
	}
	
	public CoriTimeseries(CoriTimeseries ts) {
		this(ts.estimator); //, ts.r);
		this.infectivityProfile = ts.infectivityProfile;
		ts.forEach(e -> this.add(new CoriTimeseriesEntry(e, this)));
	}

	public void inferStart(Optional<Integer> fixedValue, int length, double r) {
		CoriTimeseriesEntry first = this.first();
		if (first.isInferred()) return; // must have already done it
		for (int i=0;i<length;i++) {
			CoriTimeseriesEntry prev;
			if (fixedValue.isPresent()) {
				prev = new CoriTimeseriesEntry(first.dateValue().minusDays(1), fixedValue.get(), this, true);
			} else {
				prev = new CoriTimeseriesEntry(first.dateValue().minusDays(1), first.value()/(1+r), this, true);
			}
			first = prev;
			this.add(prev);
		}	
	}
	
	public Optional<CoriTimeseriesEntry> start() {
		CoriTimeseriesEntry start = this.first();
		while(start.isInferred()) start = start.next().orElseThrow(() -> new RuntimeException("No non inferred values. Wierd."));
		return Optional.of(start);
	}

	public void checkComplete() throws IncompleteTimeseriesException {
		for (CoriTimeseriesEntry t: this) {
			if( !t.dateValue().equals(
					t.next().map(
							s -> s.dateValue().minusDays(1)
							).orElse(t.dateValue())) // this deals with the end of the time series
					) {
				throw new IncompleteTimeseriesException(t.dateValue().plusDays(1).toString());
			}
		};
	}
	
	public CoriEstimator getEstimator() {
		return estimator;
	}
	
	public CoriTimeseries estimateForProfile(double[] infectivityProfile) {
		this.infectivityProfile = infectivityProfile;
		this.maxTau = estimator.getMaxTau();
		this.forEach(CoriTimeseriesEntry::recache);
		return this;
	}
	
	public CoriTimeseries clone() {
		return new CoriTimeseries(this); 
	}

	public int getMaxTau() {
		return maxTau;
	}

	public double[] getInfectivityProfile() {
		return  this.infectivityProfile;
	}

	
}
