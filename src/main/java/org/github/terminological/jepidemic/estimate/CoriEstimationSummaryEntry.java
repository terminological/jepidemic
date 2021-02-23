package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Optional;
import java.util.stream.Collector;

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

	Timeseries<CoriEstimationSummaryEntry> ts;
	StatSummary summary;
	LocalDate date;
	double incidence;
	int window;
	
	public CoriEstimationSummaryEntry(StatSummary summary, LocalDate endDate, int window, double incidence) {
		this.summary = summary;
		this.date = endDate;
		this.window = window;
		this.incidence = incidence;
	}
	
	public CoriEstimationSummaryEntry(DatedRtGammaEstimate datedRtGammaEstimate) {
		this(datedRtGammaEstimate, datedRtGammaEstimate.date, datedRtGammaEstimate.tau+1, datedRtGammaEstimate.incidence);
	}

	@Override
	public double getMean() {
		return summary.getMean();
	}

	@Override
	public double getSD() {
		return summary.getSD();
	}

	@Override
	public double quantile(double q) {
		return summary.quantile(q);
	}
	
	public LocalDate getEndDate() {
		return date;
	}
	
	public LocalDate getStartDate() {
		if (window-1 > 0) return date.minusDays(window-1);
		else {return date;}
	}
	
	public double getIncidence() {
		return incidence;
	}

	@Override
	public int compareTo(CoriEstimationSummaryEntry o) {
		return getEndDate().compareTo(o.getEndDate());
	}

	public int getWindow() {
		return window;
	}

	@Override
	public RNumeric incidence() {
		return RConverter.convert(incidence);
	}

	@Override
	public RDate date() {
		return RConverter.convert(date);
	}

	@SuppressWarnings("unchecked")
	@Override
	public void setTimeseries(Timeseries<?> ts) {
		this.ts = (Timeseries<CoriEstimationSummaryEntry>) ts;
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

//	public static Collector<CoriEstimationSummaryEntry, ?, RDataframe> collector(String dateCol, String incidenceCol, boolean epiEstimCompat) {
//		
//	}
	
	public static Collector<CoriEstimationSummaryEntry, ?, RDataframe> collector(String dateCol, String incidenceCol, boolean epiEstimCompat) { 
		return RConverter.dataframeCollector(
				RConverter.mapping(dateCol,s -> s.getStartDate()),
				RConverter.mapping(incidenceCol,s -> s.getIncidence()),
				RConverter.mapping("Rt.StartDate",s -> s.getStartDate()),
				RConverter.mapping("Rt.EndDate",s -> s.getEndDate()),
				RConverter.mapping("Rt.Window",s -> s.getWindow()),
				RConverter.mapping(epiEstimCompat ? "Mean(R)": "Rt.Mean",s -> s.getMean()),
				RConverter.mapping(epiEstimCompat ? "Std(R)": "Rt.SD",s -> s.getSD()),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.025(R)": "Rt.Quantile.0.025",s -> s.quantile(0.025)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.25(R)": "Rt.Quantile.0.25",s -> s.quantile(0.25)),
				RConverter.mapping(epiEstimCompat ? "Median(R)": "Rt.Quantile.0.5",s -> s.quantile(0.5)),
				RConverter.mapping(epiEstimCompat ? "Quantile.0.75(R)": "Rt.Quantile.0.75",s -> s.quantile(0.75)),
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
	
}
