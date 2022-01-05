package org.github.terminological.jepidemic.growth;

import java.util.stream.Collector;

import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.types.RDataframe;

public class Formatter {

	public static Collector<GrowthRateTimeseriesEntry, ?, RDataframe> rtAndGrowthCollector(String dateCol, String incidenceCol) { 
		return RConverter.dataframeCollector(
				RConverter.mapping(dateCol,s -> s.getDate()),
				RConverter.mapping(incidenceCol,s -> s.getIncidence()),
				RConverter.mapping("Rt.Mean",s -> s.getSummaryEstimates().getRt().getMean()),
				RConverter.mapping("Rt.SD",s -> s.getSummaryEstimates().getRt().getSD()),
				RConverter.mapping("Rt.Quantile.0.025",s -> s.getSummaryEstimates().getRt().quantile(0.025)),
				RConverter.mapping("Rt.Quantile.0.05",s -> s.getSummaryEstimates().getRt().quantile(0.05)),
				RConverter.mapping("Rt.Quantile.0.25",s -> s.getSummaryEstimates().getRt().quantile(0.25)),
				RConverter.mapping("Rt.Quantile.0.5",s -> s.getSummaryEstimates().getRt().quantile(0.5)),
				RConverter.mapping("Rt.Quantile.0.75",s -> s.getSummaryEstimates().getRt().quantile(0.75)),
				RConverter.mapping("Rt.Quantile.0.95",s -> s.getSummaryEstimates().getRt().quantile(0.95)),
				RConverter.mapping("Rt.Quantile.0.975",s -> s.getSummaryEstimates().getRt().quantile(0.975)),
				RConverter.mapping("Growth.value",s -> s.getSummaryEstimates().getGrowth().getMean()),
				RConverter.mapping("Growth.SE.value",s -> s.getSummaryEstimates().getGrowth().getSD()),
				RConverter.mapping("Growth.Quantile.0.025",s -> s.getSummaryEstimates().getGrowth().quantile(0.025)),
				RConverter.mapping("Growth.Quantile.0.05",s -> s.getSummaryEstimates().getGrowth().quantile(0.05)),
				RConverter.mapping("Growth.Quantile.0.25",s -> s.getSummaryEstimates().getGrowth().quantile(0.25)),
				RConverter.mapping("Growth.Quantile.0.5",s -> s.getSummaryEstimates().getGrowth().quantile(0.5)),
				RConverter.mapping("Growth.Quantile.0.75",s -> s.getSummaryEstimates().getGrowth().quantile(0.75)),
				RConverter.mapping("Growth.Quantile.0.95",s -> s.getSummaryEstimates().getGrowth().quantile(0.95)),
				RConverter.mapping("Growth.Quantile.0.975",s -> s.getSummaryEstimates().getGrowth().quantile(0.975)),
				RConverter.mapping("Lambda.value",s -> s.getSummaryEstimates().getLambda().getMean()),
				RConverter.mapping("Lambda.SE.value",s -> s.getSummaryEstimates().getLambda().getSD()),
				RConverter.mapping("Est.value",s -> s.getSummaryEstimates().getIncidence().getMean()),
				RConverter.mapping("Est.SE.value",s -> s.getSummaryEstimates().getIncidence().getSD()),
				RConverter.mapping("Est.Quantile.0.025",s -> s.getSummaryEstimates().getIncidence().quantile(0.025)),
				RConverter.mapping("Est.Quantile.0.05",s -> s.getSummaryEstimates().getIncidence().quantile(0.05)),
				RConverter.mapping("Est.Quantile.0.25",s -> s.getSummaryEstimates().getIncidence().quantile(0.25)),
				RConverter.mapping("Est.Quantile.0.5",s -> s.getSummaryEstimates().getIncidence().quantile(0.5)),
				RConverter.mapping("Est.Quantile.0.75",s -> s.getSummaryEstimates().getIncidence().quantile(0.75)),
				RConverter.mapping("Est.Quantile.0.95",s -> s.getSummaryEstimates().getIncidence().quantile(0.95)),
				RConverter.mapping("Est.Quantile.0.975",s -> s.getSummaryEstimates().getIncidence().quantile(0.975))
		);
	}
	
	public static Collector<GrowthRateTimeseriesEntry, ?, RDataframe> growthCollector(String dateCol, String incidenceCol) { 
		return RConverter.dataframeCollector(
				RConverter.mapping(dateCol,s -> s.getDate()),
				RConverter.mapping(incidenceCol,s -> s.getIncidence()),
				RConverter.mapping("Growth.value",s -> s.getSummaryEstimates().getGrowth().getMean()),
				RConverter.mapping("Growth.SE.value",s -> s.getSummaryEstimates().getGrowth().getSD()),
				RConverter.mapping("Growth.Quantile.0.025",s -> s.getSummaryEstimates().getGrowth().quantile(0.025)),
				RConverter.mapping("Growth.Quantile.0.05",s -> s.getSummaryEstimates().getGrowth().quantile(0.05)),
				RConverter.mapping("Growth.Quantile.0.25",s -> s.getSummaryEstimates().getGrowth().quantile(0.25)),
				RConverter.mapping("Growth.Quantile.0.5",s -> s.getSummaryEstimates().getGrowth().quantile(0.5)),
				RConverter.mapping("Growth.Quantile.0.75",s -> s.getSummaryEstimates().getGrowth().quantile(0.75)),
				RConverter.mapping("Growth.Quantile.0.95",s -> s.getSummaryEstimates().getGrowth().quantile(0.95)),
				RConverter.mapping("Growth.Quantile.0.975",s -> s.getSummaryEstimates().getGrowth().quantile(0.975)),
				RConverter.mapping("Lambda.value",s -> s.getSummaryEstimates().getLambda().getMean()),
				RConverter.mapping("Lambda.SE.value",s -> s.getSummaryEstimates().getLambda().getSD()),
				RConverter.mapping("Est.value",s -> s.getSummaryEstimates().getIncidence().getMean()),
				RConverter.mapping("Est.SE.value",s -> s.getSummaryEstimates().getIncidence().getSD()),
				RConverter.mapping("Est.Quantile.0.025",s -> s.getSummaryEstimates().getIncidence().quantile(0.025)),
				RConverter.mapping("Est.Quantile.0.05",s -> s.getSummaryEstimates().getIncidence().quantile(0.05)),
				RConverter.mapping("Est.Quantile.0.25",s -> s.getSummaryEstimates().getIncidence().quantile(0.25)),
				RConverter.mapping("Est.Quantile.0.5",s -> s.getSummaryEstimates().getIncidence().quantile(0.5)),
				RConverter.mapping("Est.Quantile.0.75",s -> s.getSummaryEstimates().getIncidence().quantile(0.75)),
				RConverter.mapping("Est.Quantile.0.95",s -> s.getSummaryEstimates().getIncidence().quantile(0.95)),
				RConverter.mapping("Est.Quantile.0.975",s -> s.getSummaryEstimates().getIncidence().quantile(0.975))
		);
	}
	
}
