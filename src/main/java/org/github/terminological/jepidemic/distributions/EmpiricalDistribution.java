package org.github.terminological.jepidemic.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.AbstractRealDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.EnumeratedRealDistribution;
import org.apache.commons.math3.distribution.IntegerDistribution;
import org.apache.commons.math3.distribution.RealDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;


public class EmpiricalDistribution extends AbstractRealDistribution {

	
	SummaryStatistics data;
	List<Double> values;
	List<Double> invCdf;
	List<Double> cdf;
	List<Double> pdf;
	
	/**
     * RNG instance used to generate samples from the distribution.
     * @since 3.1
     */
    protected final RandomGenerator random;
	
	public EmpiricalDistribution(List<Double> values) {
		super(Random.RNG);
		this.values = values;
		Collections.sort(this.values);
		this.random = new Well19937c();
		this.data = new SummaryStatistics();
		values.stream().forEach(data::addValue);
	}
	
	public EmpiricalDistribution(double[] values) {
		this(Arrays.stream(values).mapToObj(a->a).collect(Collectors.toList()));
	}
	
	@Override
	public double probability(double x) {
		return(density(x));
	}

	@Override
	public double density(double x) {
		if (x < getSupportLowerBound()) return 0;
		if (x > getSupportUpperBound()) return 0;
		if (this.pdf == null) {
			double range = data.getMax() - data.getMin();
			int elements = cdf.size();
			double binSize = range/((double) elements);
			pdf = SavitzkyGolay.convolute(cdf, SavitzkyGolay.filter(
					Math.min(25, (int) cdf.size()/10), 
					1, 1, binSize), false);
		}
		double index = (x-data.getMin())/(data.getMax() - data.getMin())*pdf.size();
		if (index<0 || index > pdf.size()) return 0;
		return fracGet(index,pdf);
	}

	@Override
	public double cumulativeProbability(double x) {
		if (x < getSupportLowerBound()) return 0;
		if (x > getSupportUpperBound()) return 1;
		if(this.cdf == null) {
			double range = data.getMax() - data.getMin();
			int elements = 1000; //Math.min(1000, values.size()/50);
			double binSize = range/((double) elements);
			cdf = new ArrayList<>(elements);
			int i = 0;
			for (double binEnd= data.getMin(); binEnd <= data.getMax(); binEnd+=binSize) {
				while (values.get(i) < binEnd) {
					i+=1;
				}
				cdf.add((double) i/((double) values.size()));
			}
		}
		double index = (x-data.getMin())/(data.getMax() - data.getMin())*cdf.size();
		if (index<0) return 0;
		if (index > values.size()) return 1;
		return fracGet(index,cdf);
	}

	private double fracGet(double index, List<Double> list) {
		int lower_x = (int) Math.floor(index);
		int upper_x = (int) Math.ceil(index);
		if (lower_x==upper_x) return list.get(upper_x);
		double upper_y = list.get(upper_x);
		double lower_y = list.get(lower_x);
		return lower_y+(upper_y-lower_y)*(index-(double) lower_x);
	}

	@Override
	public double cumulativeProbability(double x0, double x1) throws NumberIsTooLargeException {
		return cumulativeProbability(x1)-cumulativeProbability(x0);
	}

	@Override
	public double inverseCumulativeProbability(double p) throws OutOfRangeException {
		if(p<0 || p>1) throw new OutOfRangeException(0,p,1);
		return fracGet(values.size()*p, values);
	}

	@Override
	public double getNumericalMean() {
		return data.getMean();
	}

	@Override
	public double getNumericalVariance() {
		return data.getVariance();
	}

	@Override
	public double getSupportLowerBound() {
		return data.getMin();
	}

	@Override
	public double getSupportUpperBound() {
		return data.getMax();
	}

	@Override
	public boolean isSupportLowerBoundInclusive() {
		return true;
	}

	@Override
	public boolean isSupportUpperBoundInclusive() {
		return true;
	}

	@Override
	public boolean isSupportConnected() {
		return true;
	}

	public void reseedRandomGenerator(long seed) {
        random.setSeed(seed);
    }

	@Override
	public double sample() {
		int item = this.random.nextInt(this.values.size());
		return this.values.get(item);
	}

	@Override
	public double[] sample(int sampleSize) {
		if(sampleSize == values.size()) return values.stream().mapToDouble(i->i).toArray();
		double[] out = new double[sampleSize];
		for (int i = 0; i<sampleSize; i++) {
			out[i] = sample();
		}
		return out;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((values == null) ? 0 : values.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof EmpiricalDistribution))
			return false;
		EmpiricalDistribution other = (EmpiricalDistribution) obj;
		if (values == null) {
			if (other.values != null)
				return false;
		} else if (!values.equals(other.values))
			return false;
		return true;
	}

	public EmpiricalDistribution from(List<Double> input) {
		return new EmpiricalDistribution(input);
	}
	
	public EmpiricalDistribution from(double[] input) {
		return new EmpiricalDistribution(input);
	}
	
	public EnumeratedIntegerDistribution fromDiscrete(List<Integer> input) {
		int[] values = new int[input.size()];
		for (int i=0; i<input.size(); i++) values[i] = input.get(i);
		return fromDiscrete(values);
	}
	
	public EnumeratedIntegerDistribution fromDiscrete(int[] input) {
		return new EnumeratedIntegerDistribution(input);
	}
	
	public static EnumeratedRealDistribution bySamplingDiscrete(
			IntegerDistribution x,
			IntegerDistribution y,
			BiFunction<Integer,Integer,Double> fn_x_y,
			int sampleSize
			) {
		return bySamplingDiscrete(Collections.singletonList(x),Collections.singletonList(y),fn_x_y,sampleSize);
	}
	
	public static EnumeratedRealDistribution bySamplingDiscrete(
			List<? extends IntegerDistribution> x,
			List<? extends IntegerDistribution> y,
			BiFunction<Integer,Integer,Double> fn_x_y,
			int sampleSize // for each distribution...
			) {
		double[] tmp = new double[sampleSize*x.size()];
		if (x.size() != y.size()) throw new DimensionMismatchException(x.size(),y.size());
		for (int i = 0; i<x.size(); i++) {
			int[] xs = x.get(i).sample(sampleSize);
			int[] ys = y.get(i).sample(sampleSize);
			for (int j=0; j<xs.length; j++) {
				double z = fn_x_y.apply(xs[j], ys[j]);
				tmp[i*x.size()+j] = z;
			}
		}
		return new EnumeratedRealDistribution(tmp);
	}
	
	public static EmpiricalDistribution bySampling(RealDistribution x,
			RealDistribution y,
			BiFunction<Double,Double,Double> fn_x_y,
			int sampleSize
			) {
		return bySampling(Collections.singletonList(x),Collections.singletonList(y),fn_x_y,sampleSize);
	}
	
	public static EmpiricalDistribution bySampling(
			List<RealDistribution> x,
			List<RealDistribution> y,
			BiFunction<Double,Double,Double> fn_x_y,
			int sampleSize // for each distribution...
			) {
		List<Double> tmp = new ArrayList<>();
		if (x.size() != y.size()) throw new DimensionMismatchException(x.size(),y.size());
		for (int i = 0; i<x.size(); i++) {
			double[] xs = x.get(i).sample(sampleSize);
			double[] ys = y.get(i).sample(sampleSize);
			for (int j=0; j<xs.length; j++) {
				double z = fn_x_y.apply(xs[j], ys[j]);
				tmp.add(z);
			}
		}
		return new EmpiricalDistribution(tmp);
	}
	
}
