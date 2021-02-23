package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import org.github.terminological.jepidemic.Timeseries;
import org.github.terminological.jepidemic.TimeseriesEntry;
import org.github.terminological.jepidemic.gamma.GammaParameters;

import uk.co.terminological.rjava.RConverter;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RNumeric;

/**
 * A timeseries entry in a ordered straucture with accessors for the previous and next items in the time series. <br>
 * 
 * This is where the component calculations for the Cori method are performed. The calculations for the cori method are not
 * performed until the recache() method is called, which is done by the CoriEstimator::estimateForProfile method when considering a new infectivity profile.
 *  
 *
 */
public class CoriTimeseriesEntry implements TimeseriesEntry, Comparable<CoriTimeseriesEntry> {

	private LocalDate date;
	private CoriTimeseries series;
	private Optional<CoriTimeseriesEntry> next;
	private Optional<CoriTimeseriesEntry> prev;
	private boolean inferred;
	private double[] I_s = null; // the incidence at s timepoints age
	private double[] Lambda_s = null;
	
	private double value;
	
	public CoriTimeseriesEntry(TimeseriesEntry e, CoriTimeseries series) {
		this(e.date().get(), e.incidence().javaPrimitive(),  series, false);	
	}
	
	public CoriTimeseriesEntry(CoriTimeseriesEntry e, CoriTimeseries series) {
		this(e.dateValue(), e.value(), series, e.isInferred());
		this.inferred = e.inferred;
	}
	
	public CoriTimeseriesEntry(LocalDate date, double value, CoriTimeseries series, boolean inferred) {
		this.value = value;
		this.date = date;
		this.setTimeseries(series);
		this.inferred = inferred;
	}

	@Override
	public void setTimeseries(Timeseries<?> ts) {
		this.series = (CoriTimeseries) ts;
	}
	
	public Optional<CoriTimeseriesEntry> next() {
		if (next == null) next = Optional.ofNullable(series.higher(this));
		return next;
	}
	
	public Optional<CoriTimeseriesEntry> prev() {
		if (prev == null) prev = Optional.ofNullable(series.lower(this));
		return prev;
	}
	
	public double value() { return value;}
	public LocalDate dateValue() { return date;}
	
	@Override
	public RNumeric incidence() {return RConverter.convert(this.value);}

	@Override
	public RDate date() {return RConverter.convert(this.date);}

	@Override
	public int compareTo(CoriTimeseriesEntry o) {
		return this.dateValue().compareTo(o.dateValue());
	}

	public String toString() {
		return dateValue()+"\t"+value();
	}
	
	
	
	
	
	
	
	// TODO: Lots of potential strategies here based on windowing
	// could select prior based on previous step posterior with 
	
	// calculate posteriors basd on priors
	
	public DatedRtGammaEstimate posterior(DatedRtGammaEstimate prior, int tau) {
		return new GammaParameters(
				prior.getShape()+sumI_s(tau),
				1/(1/prior.getScale()+sumLambda_s(tau)))
				.withDate(tau, date, incidence().javaPrimitive(Double.NaN))
				.withPrior(prior);
		
	}
	
	// must have a prior for every tau?	
	public List<DatedRtGammaEstimate> results(List<DatedRtGammaEstimate> priors) {
		return priors.stream().map(prior -> this.posterior(prior, prior.tau)
				).collect(Collectors.toList());
	}
	
	
	// This sets up the caches.
	// It should work in any order but makes most sense to start at the beginning of the timeseries.
	
	protected void recache() {
		cacheI_s();
		cacheLambda_s();
	}
	
	public double casesInWindow(int window) {
		return sumI_s(window);
	}
	
	// I_s is the sum of incidence over the window.
	// Here we calculate it for all windows as we need that anyway for 
	// calculating the max size window eventually
	
	private double sumI_s(int tau) {
		if (this.inferred) return Double.NaN;
		if (tau < 0) return 0;
		if (I_s == null) cacheI_s();
		return I_s[tau];
	}
	
	
		
	private void cacheI_s() {
		I_s = new double[series.getMaxTau()];
		for (int tau = 0; tau<series.getMaxTau(); tau++) {
			final int tmpTau = tau;
			// sumI_s calculated here
			I_s[tau] = this.value + prev().map(ts -> ts.sumI_s(tmpTau-1)).orElse(Double.NaN);
		}
	}
	
	// Lambda_s is the sum of Lambda_t over the window.
	// Here we calculate it for all windows as we need that anyway for 
	// calculating the max size window eventually
		
	
	private double sumLambda_s(int tau) {
		if (tau < 0) return 0;
		return Lambda_s[tau];
	}
	
	
	
	private void cacheLambda_s() {
		Lambda_s = new double[series.getMaxTau()];
		for (int tau = 0; tau<series.getMaxTau(); tau++) {
			final int tmpTau = tau;
			// sumI_s calculated here
			Lambda_s[tau] = this.Lambda_t() + prev().map(ts -> ts.sumLambda_s(tmpTau-1)).orElse(Double.NaN);
		}
	}
	
	// Lambda_t is the sum of incidence * infectivity on day t over the whole infectivity profile.
	// We just calculate that for each day as none of that info is shared between days.
		
	
	private double Lambda_t = -1;
	
	private double Lambda_t() {
		if (Lambda_t < 0) Lambda_t = Lambda_tComponent(0);
		return Lambda_t;
	}
	
	private double Lambda_tComponent(int s) {
		double[] omega = this.series.getInfectivityProfile();
		if (s >= omega.length) return 0;
		return this.value*omega[s] + prev().map(ts -> ts.Lambda_tComponent(s+1)).orElse(Double.NaN);
	}

	public void setInferred() {
		this.inferred = true;
	}

	public boolean isInferred() {
		return inferred;
	}

	@Override
	public Optional<Timeseries<?>> getTimeseries() {
		return Optional.ofNullable(series);
	}

	
}
