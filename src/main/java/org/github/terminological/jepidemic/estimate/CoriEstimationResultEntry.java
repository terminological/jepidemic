package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.github.terminological.jepidemic.Timeseries;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.gamma.GammaParameters;
import org.github.terminological.jepidemic.gamma.WeightedGammaMixture;

import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;

public class CoriEstimationResultEntry implements TimeseriesEntry.Incidence, Comparable<CoriEstimationResultEntry> { 

	private CoriTimeseriesEntry tsEntry;
	//private double incidence;
	//private LocalDate date;
	private TreeSet<DatedRtGammaEstimate> posteriors;
	private int profId;
	private CoriEstimationResult series;
	private Optional<CoriEstimationResultEntry> next;
	private Optional<CoriEstimationResultEntry> prev;
	
	public CoriEstimationResultEntry(int profId, CoriTimeseriesEntry tsEntry, 
			List<DatedRtGammaEstimate> posteriorsForTaus, CoriEstimationResult series) {
		this.tsEntry = tsEntry;
		this.profId = profId;
		//this.incidence = tsEntry.value();
		//this.date = tsEntry.dateValue();
		this.posteriors = new TreeSet<>(posteriorsForTaus);
		this.series = series;
	}
	
	public CoriEstimationResultEntry(CoriEstimationResultEntry unfiltered, 
			DatedRtGammaEstimate selectedPosterior) {
		this.tsEntry = unfiltered.tsEntry;
		this.profId = unfiltered.profId;
		//this.incidence = tsEntry.value();
		//this.date = tsEntry.dateValue();
		this.posteriors = new TreeSet<DatedRtGammaEstimate>();
		this.posteriors.add(selectedPosterior);
		this.series = unfiltered.series;
	}

	public double getIncidence() {return tsEntry.value();}
	
	public LocalDate getDate() {return tsEntry.dateValue();}
	
	@Override
	public int compareTo(CoriEstimationResultEntry o) {
		//order by profile id and then date
		return (this.date().get().compareTo(o.date().get())+(this.profId-o.profId))*1000;
	}
	
	public Optional<CoriEstimationResultEntry> nextWithinBootstrap() {
		if (next == null) next = Optional.ofNullable(series.higher(this));
		if (!next.isPresent() || next.get().profId != this.profId) return Optional.empty();
		return next;
	}
	
	public Optional<CoriEstimationResultEntry> prevWithinBootstrap() {
		if (prev == null) prev = Optional.ofNullable(series.lower(this));
		if (!prev.isPresent() || prev.get().profId != this.profId) return Optional.empty();
		return prev;
	}
	
	public Optional<CoriEstimationResultEntry> lag(int i) {
		if (i==0) return Optional.of(this);
		return this.prevWithinBootstrap().flatMap(p->p.lag(i-1));
	}
	
	public DatedRtGammaEstimate selectWindowByCases(double minCases, int minWindow) {
		for (DatedRtGammaEstimate result: posteriors) {
			if (
					result.tau >= minWindow-1 &&
					this.tsEntry.casesInWindow(result.tau) > minCases
				) return result;
		}
		return posteriors.last();
	}
	
	public DatedRtGammaEstimate selectWindowByMinimumUncertainty(double timeVsRt, int minWindow) {
		timeVsRt = Math.sqrt(timeVsRt);
		DatedRtGammaEstimate out = CoriEstimationResultEntry.na(-1, this.getDate(), this.getIncidence(),this.profId);
		double score = Double.POSITIVE_INFINITY;
		for (DatedRtGammaEstimate result: posteriors) {
			double tmp = Math.pow(result.tau+1,timeVsRt)*
				Math.pow(result.convert().getSD(),1/timeVsRt);
			if (tmp < score && result.tau >= minWindow-1) {
				out = result;
				score = tmp;
			}
		}
		return out;
	}
	
	public DatedRtGammaEstimate selectMixtureCombination() {
		return selectMixtureCombination(x -> 1.0);
	}
	
	public DatedRtGammaEstimate selectMixtureCombination(Function<DatedRtGammaEstimate,Double> weighting) {
		List<DatedRtGammaEstimate> tmp = laterEstimatesForEffectiveDate();
		int tmp2 = (int) Math.round(tmp.stream().mapToInt(t -> t.tau).average().orElse(-1));
		return 
			new WeightedGammaMixture(tmp, weighting).gammaApproximation()
			.withDate(tmp2,this.dateValue(),this.getIncidence(),this.profId);
	}
	
	public DatedRtGammaEstimate selectSummaryForDayAndInfectivityProfile() {
		return series.estimator.posteriorSelectionStrategy.apply(this);
	}
	
	//public static DatedRtGammaEstimate NA = new GammaParameters(Double.NaN,Double.NaN).withDate(-1,null,0);  
	
	public static DatedRtGammaEstimate na(int window, LocalDate date, Double incidence, int profileId) {
		return new GammaParameters(Double.NaN,Double.NaN).withDate(window, date, incidence, profileId);
	}
	
	public DatedRtGammaEstimate forWindow(int window) {
		return posteriors.stream().filter(wg -> wg.tau==window-1).findFirst().orElse(na(window,this.getDate(),this.getIncidence(),this.profId));
	}
	
//	private List<DatedRtGammaEstimate> earlierEstimatesForEffectiveDate(LocalDate test) {
//		if (test.isBefore(getDate().minusDays(series.estimator.getMaxTau()/2))) return Collections.emptyList();
//		List<DatedRtGammaEstimate> out = new ArrayList<>();
//		out.addAll(this.todayEstimatesForEffectiveDate(test));
//		this.prevWithinBootstrap().ifPresent(p -> out.addAll(p.earlierEstimatesForEffectiveDate(test)));
//		return out;
//	}
	
	// This might work for multiple bootstraps... but only if we are in the first bootstrap...
	// Need bootstraps to be part of this class.... :-(
	public List<DatedRtGammaEstimate> laterEstimatesForEffectiveDate() {
		return laterEstimatesForEffectiveDate(getDate());
	}
	
	// This might work for multiple bootstraps... but only if we are in the first bootstrap...
	// Need bootstraps to be part of this class.... :-(
	private List<DatedRtGammaEstimate> laterEstimatesForEffectiveDate(LocalDate date) {
		if (this.getDate().isAfter(date.plusDays(series.estimator.getMaxTau()/2))) return Collections.emptyList();
		List<DatedRtGammaEstimate> out = new ArrayList<>();
		out.addAll(this.todayEstimatesForEffectiveDate(date));
		this.nextWithinBootstrap().ifPresent(p -> out.addAll(p.laterEstimatesForEffectiveDate(date)));
		return out;
	}
	
	@Override
	public RNumeric incidence() {
		return tsEntry.incidence();
	}
	
	
	
	public LocalDate dateValue() {
		return getDate();
	}
	
	@Override
	public RDate date() {
		return RConverter.convert(getDate());
	}
	
	public RInteger infectivityProfileId() {
		return RConverter.convert(profId);
	}

//	public DatedRtEstimate posterior(int window) {
//		return results.stream().filter(wg -> wg.tau==window-1).findFirst().get();
//	}
//	
	private List<DatedRtGammaEstimate> todayEstimatesForEffectiveDate(LocalDate effective) {
		return posteriors.stream().filter(wg -> wg.getEffectiveDate().equals(effective)).collect(Collectors.toList());
	}
	
	// prior estimates Rt
	public List<DatedRtGammaEstimate> adaptiveNextPrior(double factor) {
		return this.posteriors.stream().map(
				r -> r.wider(factor).orElse(series.estimator.defaultPrior(r.tau, r.date, r.profileId)))
				.collect(Collectors.toList());
	}
	
	
	public String toString() {
		return getDate().toString()+" "+getIncidence()+" "+forWindow(7).toString();
	}

	public void setSeries(CoriEstimationResult coriEstimationResult) {
		series=coriEstimationResult;
	}

	@Override
	public void setTimeseries(Timeseries<?> ts) {
		this.series = (CoriEstimationResult) ts;
	}

	@Override
	public Optional<CoriEstimationResultEntry> next() {
		return nextWithinBootstrap();
	}

	@Override
	public Optional<CoriEstimationResultEntry> prev() {
		return prevWithinBootstrap();
	}

	
	

	@Override
	public Optional<Timeseries<?>> getTimeseries() {
		return Optional.ofNullable(series);
	}

	

	public Integer getInfProfId() {
		return this.profId;
	}

	public TreeSet<DatedRtGammaEstimate> posteriors() {
		return this.posteriors;
	}

	

//	public List<DatedRtGammaEstimate> gammaMixturePrior() {
//		return results.stream().map(rEst -> {
//			List<DatedRtGammaEstimate> tmpList = this.earlierEstimatesForEffectiveDate(rEst.getEffectiveDate());
//			GammaParameters tmp = new WeightedGammaMixture(tmpList, x->(double) x.getWindow()).gammaApproximation();
//			DatedRtGammaEstimate out = tmp.withDate(-1, rEst.date, rEst.incidence).withPrior(rEst.prior);
//			if (!out.isDefined()) return series.estimator.defaultPrior(rEst.tau, rEst.date);
//			return out;
//		}).collect(Collectors.toList());
//	}

}