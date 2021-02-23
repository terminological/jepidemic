package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.github.terminological.jepidemic.Timeseries;

public class CoriEstimationResult extends Timeseries<CoriEstimationResultEntry> {
	
	private HashMap<LocalDate,TreeSet<CoriEstimationResultEntry>> byDateIndex = new HashMap<>();
	public CoriEstimator estimator;
	// TODO shoudl have infectivity profile here
	
	public CoriEstimationResult(CoriEstimator estimator) {
		super((e1,e2) -> e1.compareTo(e2));
		this.estimator = estimator;
		
	}
	
	public synchronized boolean addAll(Collection<? extends CoriEstimationResultEntry> entries) {
		boolean out = super.addAll(entries);
		entries.forEach(e -> {
			if (!byDateIndex.containsKey(e.dateValue())) byDateIndex.put(e.dateValue(),  new TreeSet<>());
			this.byDateIndex.get(e.dateValue()).add(e);
			e.setSeries(this);
		});
		//N.b. prev / next methods become odd after this point.
		return out;
	}
	
	public synchronized boolean add(CoriEstimationResultEntry e) {
		boolean out = super.add(e);
		if (!byDateIndex.containsKey(e.dateValue())) byDateIndex.put(e.dateValue(),  new TreeSet<>());
		this.byDateIndex.get(e.dateValue()).add(e);
		e.setSeries(this);
		return out;
	}
	
	public CoriEstimationSummary collectBootstraps() {
		if (estimator.combiningStrategy==null) throw new RuntimeException("Combining strategy is not set");
		CoriEstimationSummary out = new CoriEstimationSummary();
		byDateIndex.entrySet().stream().map(kv -> {
			List<DatedRtGammaEstimate> tmp = kv.getValue()
					.stream()
					.map(c -> c.selectSummaryForDayAndInfectivityProfile())
					.collect(Collectors.toList());
			return estimator.combiningStrategy.apply(kv.getKey(), tmp);
		}).forEach(out::add);
		return out;
	}
}
