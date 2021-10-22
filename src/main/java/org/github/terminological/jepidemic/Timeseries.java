package org.github.terminological.jepidemic;

import java.util.Comparator;
import java.util.TreeSet;

public class Timeseries<T extends TimeseriesEntry> extends TreeSet<T> {

	public Timeseries(Comparator<T> comparator) {
		super(comparator);
	};
	
	public static <X extends TimeseriesEntry> Timeseries<X> of(Comparator<X> comparator) {
		return new Timeseries<X>(comparator);
	}
	
	public void checkComplete() throws IncompleteTimeseriesException {
		for (TimeseriesEntry t: this) {
			if( !t.getDate().equals(
					t.next().map(
							s -> s.getDate().minusDays(1)
							).orElse(t.getDate())) // this deals with the end of the time series
					) {
				throw new IncompleteTimeseriesException(t.getDate().plusDays(1).toString());
			}
		};
	}

	
	
	
	
}
