package org.github.terminological.jepidemic.estimate;

import java.util.Collection;
import java.util.Optional;

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
				prev = new CoriTimeseriesEntry(first.dateValue().minusDays(1), first.value()/Math.exp(r), this, true);
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
	
	@Override
	public boolean add(CoriTimeseriesEntry e) {
		if (this.contains(e)) throw new RuntimeException("Timeseries already contains this value: "+e.toString());
		return super.add(e);
	}

	@Override
	public boolean addAll(Collection<? extends CoriTimeseriesEntry> c) {
		return c.stream().map(this::add).anyMatch(b -> b.equals(Boolean.TRUE));
	}

	
}
