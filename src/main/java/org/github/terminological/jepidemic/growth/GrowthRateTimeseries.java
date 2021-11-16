package org.github.terminological.jepidemic.growth;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.github.terminological.jepidemic.Timeseries;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;
import org.github.terminological.jepidemic.growth.EstimateMetadata.GrowthMetadata;
import org.github.terminological.jepidemic.growth.EstimateMetadata.RtFromGrowthMetadata;

public class GrowthRateTimeseries extends Timeseries<GrowthRateTimeseriesEntry> implements Cloneable {

	private GrowthRateEstimator estimator;
	
	//double r;
	
	//Empty timeseries constructor
	public GrowthRateTimeseries(GrowthRateEstimator estimator) { //, double initialGrowth) {
		super((c1,c2) -> c1.compareTo(c2));
		this.estimator= estimator;
	}
	
	//Copy constructor
	public GrowthRateTimeseries(GrowthRateTimeseries ts) {
		this(ts.estimator); //, ts.r);
		ts.forEach(e -> this.add(new GrowthRateTimeseriesEntry(e, this)));
	}

	public void inferStart(Optional<Integer> fixedValue, int length, double r) {
		GrowthRateTimeseriesEntry first = this.first();
		if (first.isInferred()) return; // must have already done it
		for (int i=0;i<length;i++) {
			GrowthRateTimeseriesEntry prev;
			if (fixedValue.isPresent()) {
				prev = new GrowthRateTimeseriesEntry(first.dateValue().minusDays(1), fixedValue.get(), this, true);
			} else {
				prev = new GrowthRateTimeseriesEntry(first.dateValue().minusDays(1), first.value()/Math.exp(r), this, true);
			}
			first = prev;
			this.add(prev);
		}	
	}
	
	public Optional<GrowthRateTimeseriesEntry> start() {
		GrowthRateTimeseriesEntry start = this.first();
		while(start.isInferred()) start = start.next().orElseThrow(() -> new RuntimeException("No non inferred values. Wierd."));
		return Optional.of(start);
	}

	public GrowthRateEstimator getEstimator() {
		return estimator;
	}
	
	public GrowthRateTimeseries estimate() {
		this.forEach(GrowthRateTimeseriesEntry::getDerivedEstimates);
		return this;
	}
	
	public GrowthRateTimeseries clone() {
		return new GrowthRateTimeseries(this); 
	}

	@Override
	public boolean add(GrowthRateTimeseriesEntry e) {
		if (this.contains(e)) throw new RuntimeException("Timeseries already contains this value: "+e.toString());
		return super.add(e);
	}

	@Override
	public boolean addAll(Collection<? extends GrowthRateTimeseriesEntry> c) {
		return c.stream().map(this::add).anyMatch(b -> b.equals(Boolean.TRUE));
	}

	private ExtendedGammaDistribution firstLambda() {
		double initial = Math.max(1,this.first().value());
		return ExtendedGammaDistribution.fromMoments(initial, initial);
	}
	
	List<GrowthMetadata> growthSupport;
	
	public List<GrowthMetadata> getGrowthSupport() {
		if (growthSupport == null) {
		int maxTau = estimator.getMaxTau();
		int minTau = estimator.getMinTau();
		growthSupport = IntStream.range(minTau,maxTau).mapToObj(
				tau -> new GrowthMetadata(tau, estimator.getInitialIncidence().orElse(GrowthRateTimeseries.this.firstLambda())) 
				).collect(Collectors.toList());
		}
		return growthSupport;
	}

	List<RtFromGrowthMetadata> rtFromGrowthSupport;
	
	public List<RtFromGrowthMetadata> getRtFromGrowthSupport() {
		if (rtFromGrowthSupport == null) {
			List<GrowthMetadata> gs = getGrowthSupport();
			rtFromGrowthSupport = getEstimator().infProf.stream().flatMap(ip -> gs.stream().map(g -> g.extend(ip))).collect(Collectors.toList());
		}
		return rtFromGrowthSupport;
	}
}
