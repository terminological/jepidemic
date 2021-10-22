package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Optional;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.github.terminological.jepidemic.Timeseries;
import org.github.terminological.jepidemic.TimeseriesEntry;

import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;

public class CoriBootstrappedTimeseriesEntry extends CoriTimeseriesEntry implements TimeseriesEntry.RateIncidence {

	private CoriBootstrappedTimeseries rtSeries;
	private PoissonDistribution pois;
	private int sample = 0;
	private double rate;
	
	public CoriBootstrappedTimeseriesEntry(TimeseriesEntry.Rate e, CoriBootstrappedTimeseries rtSeries) {
		// use poisson rate as the incidence
		this(e.date().get(), e.rate().javaPrimitive(),  rtSeries, false);
	}
	
	public CoriBootstrappedTimeseriesEntry(CoriBootstrappedTimeseriesEntry rate, int incidence, int sample, CoriBootstrappedTimeseries rtSeries) {
		this(rate.getDate(), incidence, rtSeries, rate.isInferred());
		this.rate = rate.rate;
		this.sample = sample;
		
	}
		
	public double rateValue() {
		return rate;
	}

	public CoriBootstrappedTimeseriesEntry(LocalDate date, double value,
			CoriTimeseries series, boolean inferred) {
		super(date, value, series, inferred);
		this.rate = value;
		if(rate > 0) {
			this.pois = new PoissonDistribution(rate);
		} else {
			this.pois = null;
		}
		this.sample = 0;
	}

	@Override
	public Optional<Timeseries<?>> getTimeseries() {
		return Optional.ofNullable(rtSeries);
	}

	

	

	public int compareTo(CoriBootstrappedTimeseriesEntry c2) {
		return this.getDate().compareTo(c2.getDate());
	}
	
	public CoriTimeseries sampleInto(CoriBootstrappedTimeseries ts, int samples) {
		int[] s;
		if (pois != null) {
			s = pois.sample(samples);
		} else {
			s = new int[samples];
		}
		for(int i=0; i < samples; i++) {
			CoriBootstrappedTimeseriesEntry tmp = new CoriBootstrappedTimeseriesEntry(this,s[i], i,ts);
			ts.add(tmp);
		}
		return(ts);
		
	}

	@Override
	public RInteger sample() {
		return RConverter.convert(sample);
	}

	@Override
	public RNumeric rate() {
		return RConverter.convert(rate);
	}
	
	

}
