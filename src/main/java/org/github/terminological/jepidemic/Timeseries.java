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
	
	
	
}
