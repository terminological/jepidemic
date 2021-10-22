package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collector;
import java.util.stream.Collectors;

import org.github.terminological.jepidemic.RtTimeseriesEntry;
import org.github.terminological.jepidemic.Timeseries;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.gamma.StatSummary;

import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.types.RDataframe;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;

public class CoriEstimationSummaryEntry implements StatSummary, RtTimeseriesEntry, Comparable<CoriEstimationSummaryEntry> {

	CoriEstimationSummary ts;
	public LocalDate date;
	List<DatedRtGammaEstimate> posteriorsForProfiles = new ArrayList<>();
	StatSummary statSummary;
		
	public CoriEstimationSummaryEntry(List<DatedRtGammaEstimate> datedRtGammaEstimates, CoriEstimationSummary summary) {
		this.posteriorsForProfiles = datedRtGammaEstimates;
		List<LocalDate> dates = datedRtGammaEstimates.stream().map(de -> de.date).distinct().collect(Collectors.toList());
		if(dates.size() > 1) throw new RuntimeException("Summary has gone wrong");
		date = dates.get(0);
		this.ts = summary;
		if (ts.getEstimator().combiningStrategy != null) {
			this.statSummary = ts.getEstimator().combiningStrategy.apply(this);
		}
	}
	
	@Override
	public double getMean() {
		if (statSummary == null) throw new RuntimeException("Summary has not been calculated. No combining strategy defined.");
		return statSummary.getMean();
	}

	@Override
	public double getSD() {
		if (statSummary == null) throw new RuntimeException("Summary has not been calculated. No combining strategy defined.");
		return statSummary.getSD();
	}

	@Override
	public double quantile(double q) {
		if (statSummary == null) throw new RuntimeException("Summary has not been calculated. No combining strategy defined.");
		return statSummary.quantile(q);
	}
	
	public LocalDate getEndDate() {
		return date;
	}
	
	public LocalDate getStartDate() {
		if (getWindow()-1 > 0) return date.minusDays(getWindow()-1);
		else {return date;}
	}
	
	public double getIncidence() {
		return meanIncidence();
	}

	@Override
	public int compareTo(CoriEstimationSummaryEntry o) {
		return getEndDate().compareTo(o.getEndDate());
	}

	public int getWindow() {
		return meanWindow();
	}

	@Override
	public RNumeric incidence() {
		return RConverter.convert(getIncidence());
	}
	
	public RNumeric rate() {
		return RConverter.convert(getRate());
	}
	
	private double getRate() {
		return meanIncidence();
	}

	@Override
	public RDate date() {
		return RConverter.convert(date);
	}

	@Override
	public void setTimeseries(Timeseries<?> ts) {
		this.ts = (CoriEstimationSummary) ts;
	}

	@Override
	public Optional<? extends TimeseriesEntry> next() {
		if (this.ts == null) return Optional.empty();
		return Optional.ofNullable(ts.higher(this));
	}

	@Override
	public Optional<? extends TimeseriesEntry> prev() {
		if (this.ts == null) return Optional.empty();
		return Optional.ofNullable(ts.lower(this));
	}

	@Override
	public RInteger window() {
		return RConverter.convert(getWindow());
	}

	@Override
	public RDate endDate() {
		return RConverter.convert(getEndDate());
	}

	@Override
	public RDate startDate() {
		return RConverter.convert(getStartDate());
	}

	@Override
	public RNumeric mean() {
		return RConverter.convert(getMean());
	}

	@Override
	public RNumeric sd() {
		return RConverter.convert(getSD());
	}
	
	public List<DatedRtGammaEstimate> estimatesForDate() {
		return this.posteriorsForProfiles;
	}

//	public static Collector<CoriEstimationSummaryEntry, ?, RDataframe> collector(String dateCol, String incidenceCol, boolean epiEstimCompat) {
//		
//	}
	
	public static Collector<CoriEstimationSummaryEntry,?,RDataframe> collector(String dateCol, String incidenceCol, boolean epiEstimCompat) {
		return RConverter.flatteningDataframeCollector(
				RConverter.flatMapping(
						cere->cere.estimatesForDate().stream(),
						RConverter.mapping("Rt.StartDate", pp->pp.getStartDate()),
						RConverter.mapping("Rt.EndDate", pp->pp.getEndDate()),
						RConverter.mapping("Rt.Window",pp -> pp.getWindow()),
						RConverter.mapping(epiEstimCompat ? "Mean(R)": "Rt.Mean",pp -> pp.convert().getMean()),
						RConverter.mapping(epiEstimCompat ? "Std(R)": "Rt.SD",pp -> pp.convert().getSD()),
						RConverter.mapping("Rt.Shape",pp -> pp.getShape()),
						RConverter.mapping("Rt.Scale",pp -> pp.getScale()),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.025(R)": "Rt.Quantile.0.025",pp -> pp.convert().quantile(0.025)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.05(R)": "Rt.Quantile.0.05",pp -> pp.convert().quantile(0.05)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.25(R)": "Rt.Quantile.0.25",pp -> pp.convert().quantile(0.25)),
						RConverter.mapping(epiEstimCompat ? "Median(R)": "Rt.Quantile.0.5",pp -> pp.convert().quantile(0.5)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.75(R)": "Rt.Quantile.0.75",pp -> pp.convert().quantile(0.75)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.95(R)": "Rt.Quantile.0.95",pp -> pp.convert().quantile(0.95)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.975(R)": "Rt.Quantile.0.975",pp -> pp.convert().quantile(0.975)),
						RConverter.mapping("Rt.prior.Mean",pp -> pp.getPrior().map(x -> x.convert().getMean()).orElse(Double.NaN)),
						RConverter.mapping("Rt.prior.SD",pp -> pp.getPrior().map(x -> x.convert().getSD()).orElse(Double.NaN)),
						RConverter.mapping("Rt.ProfileId",pp -> pp.getProfileId())
				),
				RConverter.mapping(dateCol,cere -> cere.date().get()),
				RConverter.mapping(incidenceCol,cere -> cere.incidence().get())
				
				
		);
	}
	
//	@Deprecated
//	public List<DatedRtGammaEstimate> welchSatterthwaiteNextPrior(double factor) {
//		return results.stream().map(rEst -> {
//			List<DatedRtGammaEstimate> tmpList = this.earlierEstimatesForEffectiveDate(rEst.getEffectiveDate());
//			GammaParameters tmp = GammaParameters.welchSatterthwaiteCombination(tmpList);
//			DatedRtGammaEstimate out = tmp.withDate(-1, rEst.date, rEst.incidence).withPrior(rEst.prior);
//			if (!out.isDefined()) return series.estimator.defaultPrior(rEst.tau, rEst.date);
//			return out.wider(factor).get();
//		}).collect(Collectors.toList());
//	}
	
	public static Collector<CoriEstimationSummaryEntry, ?, RDataframe> summaryCollector(String dateCol, String incidenceCol, boolean epiEstimCompat) { 
		return RConverter.dataframeCollector(
				RConverter.mapping(dateCol,s -> s.getStartDate()),
				RConverter.mapping(incidenceCol,s -> s.getIncidence()),
				RConverter.mapping("Rt.StartDate",s -> s.getStartDate()),
				RConverter.mapping("Rt.EndDate",s -> s.getEndDate()),
				RConverter.mapping("Rt.Window",s -> s.getWindow()),
				RConverter.mapping(epiEstimCompat ? "Mean(R)": "Rt.Mean",s -> s.getMean()),
				RConverter.mapping(epiEstimCompat ? "Std(R)": "Rt.SD",s -> s.getSD()),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.025(R)": "Rt.Quantile.0.025",s -> s.quantile(0.025)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.05(R)": "Rt.Quantile.0.05",s -> s.quantile(0.05)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.25(R)": "Rt.Quantile.0.25",s -> s.quantile(0.25)),
				RConverter.mapping(epiEstimCompat ? "Median(R)": "Rt.Quantile.0.5",s -> s.quantile(0.5)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.75(R)": "Rt.Quantile.0.75",s -> s.quantile(0.75)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.95(R)": "Rt.Quantile.0.95",s -> s.quantile(0.95)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.975(R)": "Rt.Quantile.0.975",s -> s.quantile(0.975))						
		);
	}

	/*
	 * colnames(currentRtSlow$rt14)
 [1] "statistic"                            "type"                                 "code"                                 "codeType"                             "source"                               "subgroup"                             "ageCat"                              
 [8] "gender"                               "name"                                 "date"                                 "value"                                "population"                           "Anomaly"                              "Imputed.value"                       
[15] "Imputed"                              "RollMean.value"                       "Window.RollMean.value"                "t_start"                              "t_end"                                "Mean(R)"                              "Std(R)"                              
[22] "Quantile.0.025(R)"                    "Quantile.0.05(R)"                     "Quantile.0.25(R)"                     "Median(R)"                            "Quantile.0.75(R)"                     "Quantile.0.95(R)"                     "Quantile.0.975(R)"                   
[29] "errors"                               "Window.R"                             "Anomaly.R"                            "Est.log(value + 1)"                   "Est.SE.log(value + 1)"                "Slope.log(value + 1)"                 "Slope.SE.log(value + 1)"             
[36] "Ratio.log(value + 1)"                 "Ratio.SE.log(value + 1)"              "Window.log(value + 1)"                "Growth.value"                         "Growth.SE.value"                      "Growth.ProbPos.value"                 "interceptDate"                       
[43] "doublingTime"                         "doublingTime.Quantile.0.025"          "doublingTime.Quantile.0.25"           "doublingTime.Quantile.0.75"           "doublingTime.Quantile.0.975"          "Est.value"                            "Est.Quantile.0.025.value"            
[50] "Est.Quantile.0.25.value"              "Est.Quantile.0.5.value"               "Est.Quantile.0.75.value"              "Est.Quantile.0.975.value"             "Growth.windowed.value"                "Growth.windowed.SE.value"             "Growth.windowed.Window.value"        
[57] "doublingTime.windowed"                "doublingTime.windowed.Quantile.0.025" "doublingTime.windowed.Quantile.0.25"  "doublingTime.windowed.Quantile.0.75"  "doublingTime.windowed.Quantile.0.975" "Growth.windowed.ProbPos.value"        "Volatility.28.Mean(R)"               
> colnames(currentDataset$finalDataset)
	 */
	
	@Override
	public Optional<Timeseries<?>> getTimeseries() {
		return Optional.ofNullable(ts);
	}

	public static Collector<CoriEstimationSummaryEntry, ?, RDataframe> poissonSummaryCollector(String dateCol, String poissonRateCol, boolean epiEstimCompat) {
		return RConverter.dataframeCollector(
				RConverter.mapping(dateCol,s -> s.getStartDate()),
				RConverter.mapping(poissonRateCol,s -> s.getRate()),
				RConverter.mapping("Rt.StartDate",s -> s.getStartDate()),
				RConverter.mapping("Rt.EndDate",s -> s.getEndDate()),
				RConverter.mapping("Rt.Window",s -> s.getWindow()),
				RConverter.mapping(epiEstimCompat ? "Mean(R)": "Rt.Mean",s -> s.getMean()),
				RConverter.mapping(epiEstimCompat ? "Std(R)": "Rt.SD",s -> s.getSD()),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.025(R)": "Rt.Quantile.0.025",s -> s.quantile(0.025)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.05(R)": "Rt.Quantile.0.05",s -> s.quantile(0.05)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.25(R)": "Rt.Quantile.0.25",s -> s.quantile(0.25)),
				RConverter.mapping(epiEstimCompat ? "Median(R)": "Rt.Quantile.0.5",s -> s.quantile(0.5)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.75(R)": "Rt.Quantile.0.75",s -> s.quantile(0.75)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.95(R)": "Rt.Quantile.0.95",s -> s.quantile(0.95)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.975(R)": "Rt.Quantile.0.975",s -> s.quantile(0.975))						
		);
	}

	public static Collector<CoriEstimationSummaryEntry,?,RDataframe>  poissonCollector(String dateCol, String poissonRateCol, boolean epiEstimCompat) {
		return RConverter.flatteningDataframeCollector(
				RConverter.flatMapping(
						cere->cere.estimatesForDate().stream(),
						RConverter.mapping("Rt.StartDate", pp->pp.getStartDate()),
						RConverter.mapping("Rt.EndDate", pp->pp.getEndDate()),
						RConverter.mapping("Rt.Window",pp -> pp.getWindow()),
						RConverter.mapping(epiEstimCompat ? "Mean(R)": "Rt.Mean",pp -> pp.convert().getMean()),
						RConverter.mapping(epiEstimCompat ? "Std(R)": "Rt.SD",pp -> pp.convert().getSD()),
						RConverter.mapping("Rt.Shape",pp -> pp.getShape()),
						RConverter.mapping("Rt.Scale",pp -> pp.getScale()),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.025(R)": "Rt.Quantile.0.025",pp -> pp.convert().quantile(0.025)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.05(R)": "Rt.Quantile.0.05",pp -> pp.convert().quantile(0.05)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.25(R)": "Rt.Quantile.0.25",pp -> pp.convert().quantile(0.25)),
						RConverter.mapping(epiEstimCompat ? "Median(R)": "Rt.Quantile.0.5",pp -> pp.convert().quantile(0.5)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.75(R)": "Rt.Quantile.0.75",pp -> pp.convert().quantile(0.75)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.95(R)": "Rt.Quantile.0.95",pp -> pp.convert().quantile(0.05)),
						RConverter.mapping(epiEstimCompat ? "Quantile.0.975(R)": "Rt.Quantile.0.975",pp -> pp.convert().quantile(0.975)),
						RConverter.mapping("Rt.prior.Mean",pp -> pp.getPrior().map(x -> x.convert().getMean()).orElse(Double.NaN)),
						RConverter.mapping("Rt.prior.SD",pp -> pp.getPrior().map(x -> x.convert().getSD()).orElse(Double.NaN)),
						RConverter.mapping("Rt.ProfileId",pp -> pp.getProfileId())
				),
				RConverter.mapping(dateCol,cere -> cere.date().get()),
				RConverter.mapping(poissonRateCol,cere -> cere.rate().get())
				
		);
	}
	
	private int meanWindow() {
		return (int) Math.round(this.posteriorsForProfiles.stream().mapToInt(e ->e.getWindow()).average().orElse(-1));	
	}
	
	private double meanIncidence() {
		return this.posteriorsForProfiles.stream().mapToDouble(e ->e.getIncidence()).average().orElse(Double.NaN);
	}

	public void merge(CoriEstimationSummaryEntry t) {
		this.posteriorsForProfiles.addAll(t.posteriorsForProfiles);
		
	}
	
	
}
