package org.github.terminological.jepidemic;

import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;

import uk.co.terminological.rjava.RClass;
import uk.co.terminological.rjava.RMethod;
import uk.co.terminological.rjava.types.RNumericVector;

/**
 * A holder of a discrete probability distribution of the infectivity profile 
 *
 */
@RClass
public class InfectivityProfile extends EnumeratedIntegerDistribution {

	RNumericVector infProf;
	private int id;
	
	/**
	 * @param discretePdf - A discrete pdf of probabilities over time of secondary infection given priamry infection where probability at time 0 is 0
	 * @param id - a numeric index or id for this profile
	 */
	@RMethod
	public InfectivityProfile(RNumericVector discretePdf, int id) {
		super(
				IntStream.rangeClosed(1, discretePdf.size()).toArray(),  //sequence of ints starting at 1
				discretePdf.javaPrimitive(0) // the 
				);
		this.infProf = discretePdf;
		this.id = id;
	}
	
	public double[] profile() {
		return infProf.javaPrimitive(Double.NaN);
	}

	public int getId() {
		return id;
	}

	public int length() {
		return infProf.size();
	}
	
}
