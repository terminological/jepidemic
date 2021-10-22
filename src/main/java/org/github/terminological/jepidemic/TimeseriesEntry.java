package org.github.terminological.jepidemic;

import java.time.LocalDate;
import java.util.Optional;
import java.util.stream.Stream;

import uk.co.terminological.rjava.RName;
import uk.co.terminological.rjava.types.RDate;
import uk.co.terminological.rjava.types.RInteger;
import uk.co.terminological.rjava.types.RNumeric;

/**
 * A timeseries entry is simply an adaptor place holder for a dataframe with numeric columns "I" and "date"
 *
 */
public interface TimeseriesEntry {
	
	
	@RName("date") RDate date();
	
	public default LocalDate getDate() {return date().get();} 
	
	
	void  setTimeseries(Timeseries<?> ts);
	public Optional<Timeseries<?>> getTimeseries();
	public Optional<? extends TimeseriesEntry> next();
	public Optional<? extends TimeseriesEntry> prev();
	
	public default Stream<? extends TimeseriesEntry> laggingWindow(int count) {
		if (count == 1) return this.prev().stream();
		return Stream.concat(
				this.prev().stream(),
				this.prev().stream().flatMap(p -> p.laggingWindow(count-1))
				);
	};
	
	public default Stream<? extends TimeseriesEntry> laggingWindowInclusive(int count) {
		if (count == 0) return Stream.of(this);
		return Stream.concat(
				Stream.of(this),
				this.prev().stream().flatMap(p -> p.laggingWindowInclusive(count-1))
				);
	};
	
	
	public static interface Incidence extends TimeseriesEntry {
		@RName("I") RNumeric incidence();
		public default double getIncidence() {return incidence().get();}
	}
	
	public static interface Rate extends TimeseriesEntry {
		@RName("pois") RNumeric rate();
		public default double getPoissonRate() {return rate().get();}
	}
	
	public static interface RateIncidence extends TimeseriesEntry.Rate, TimeseriesEntry.Incidence {
		@RName("sample") RInteger sample();
		public default int getSample() {return sample().get();}
	}
	
}
