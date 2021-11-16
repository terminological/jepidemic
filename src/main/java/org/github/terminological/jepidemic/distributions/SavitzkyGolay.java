package org.github.terminological.jepidemic.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.Function;

public class SavitzkyGolay {

	//TODO: https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
	
	//https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Appendix

	static double[] smooth_5_cubic() {return new double[]{-3D/35,12D/35,17D/35,12D/35,-3D/35};}
	
	static double[] smooth_7_cubic() {return new double[]{-2D/21,3D/21,6D/21,7D/21,6D/21,3D/21,-2D/21};}
	static double[] smooth_9_cubic() {return new double[]{-21D/231,14D/231,39D/231,54D/231,59D/231,54D/231,39D/231,14D/231,-21D/231};}

	static double[] derivative_5_quad(double w) {return new double[]{-2D/(10*w),-1D/(10*w),0,1D/(10*w),2D/(10*w)};}
	static double[] derivative_7_quad(double w) {return new double[]{-3D/(28*w),-2D/(28*w),-1D/(28*w),0,1D/(28*w),2D/(28*w),3D/(10*w)};}

	static double[] derivative_5_quartic(double w) {return new double[]{1D/(12*w),-8D/(12*w),0,8D/(12*w),-1D/(12*w)};}
	static double[] derivative_7_quartic(double w) {return new double[]{22D/(252*w),-67D/(252*w),-58D/(252*w),0,58D/(252*w),67D/(252*w),-22D/(252*w)};}

	static double[] smooth_N_sliding(int n) { 
		double[] tmp = new double[n];
		Arrays.fill(tmp, 1D/n);
		return tmp;}
	
	//https://www.people.iup.edu/jford/courses/CHEM421/Resources/GeneralLeastSquaresSmoothingAndDifferentiationByTheConvolutionMethod.pdf
	
	/*
	 * Calculates the Gram Polynomial (s=0) or its s'th derivative evaluated at i, order k, over 2m+1 points 
	 */
	private static double gramPoly(int i, int m, int k, int s) {
		if (k<=0) return (k==0 && s==0) ? 1 : 0;
		return (4D*k-2D)/(k*(2D*m-k+1D))*(i*gramPoly(i,m,k-1,s)
				+s*gramPoly(i,m,k-1,s-1))
				-((k-1)*(2D*m+k))/(k*(2D*m-k+1))*gramPoly(i,m,k-2,s);
	}
	
	/*
	 * Calculates the generalised factorial (a)(a+1)..(a-b+1)
	 */
	private static double genFact(int a, int b) {
		double out = 1D;
		for (int i=a-b+1; i<=a; i++) {
			out *= i;
		}
		return out;
	}
	
	/*
	 * Calculates the weight of the i'th data point for the t'th Least-Square point of the s'th derivative, over 2m+1 points, order n 
	 */
	private static double weight(int i, int t, int m, int n, int s) {
		double sum = 0D;
		for (int k = 0; k<=n; k++) {
			sum += (2D*k+1D)*(genFact(2*m,k)/genFact(2*m+k+1,k+1))*gramPoly(i,m,k,0)*gramPoly(t,m,k,s);
		}
		return sum;
	}
	
	/*
	 * generate a Savitzky Golay filter for arbitrary set of parameters 
	 */
	static double[] filter(int width, int order, int derivative, double interval) {
		int m = (width-1)/2;
		double[] out = new double[2*m+1];
		for (int i = -m; i<=m; i++) {
			double tmp = weight(i,0,m,order,derivative);
			out[m+i] = tmp/Math.pow(interval,derivative);
		}
		return out;
	}
	
	
	static <X> Double convolute(List<X> input, double[] filter, boolean circular, int position, Function<X,Double> toDouble) {
		List<X> window = circular ? circular(input, filter.length, position) : tailed(input, filter.length, position);
		Double collect = 0D;
		for (int i=0; i<filter.length; i=i+1) {
			collect += toDouble.apply(window.get(i))*filter[i];
			//TODO: deal with null etc
		}
		return collect;
	}
	
	static Double convolute(List<Double> input, double[] filter, boolean circular, int position) {
		List<Double> window = circular ? circular(input, filter.length, position) : tailed(input, filter.length, position);
		Double collect = 0D;
		for (int i=0; i<filter.length; i++) {
			collect += window.get(i)*filter[i];
			//TODO: deal with null etc
		}
		return collect;
	}
	
	static List<Double> convolute(List<Double> input, double[] filter, boolean circular) {
		int size = filter.length;
		List<Double> out = new ArrayList<>();
		Iterator<List<Double>> window = circular ? circular(input, size) : tailed(input, size);
		while (window.hasNext()) {
			List<Double> tmp = window.next();
			Double collect = 0D;
			for (int i=0; i<filter.length; i++) {
				collect += tmp.get(i)*filter[i];
				//TODO: deal with null etc
			}
			out.add(collect);
		}
		return out;
	}

	static <X> Iterator<List<X>> circular(final List<X> input, int width) { 

		if (input.size() < width) throw new ArrayIndexOutOfBoundsException("window width larger than list size");
		return new Iterator<List<X>>() {

			int position = 0;

			@Override
			public boolean hasNext() {
				return position < input.size();
			}

			@Override
			public List<X> next() {
				return circular(input,width,position++);
			}

		};
	}
	
	static <X> List<X> circular(final List<X> input, int width, int position) {
		if (position > input.size() || position < 0) throw new ArrayIndexOutOfBoundsException("window position outside of list");
		int i = width/2;
		int start = position-i;
		int end = position+i+1;
		List<X> out = new ArrayList<>();
		if (start < 0) {
			out.addAll(input.subList(input.size()+start, input.size()));
			out.addAll(input.subList(0, end));
		} else if (end > input.size()) {
			out.addAll(input.subList(start, input.size()));
			out.addAll(input.subList(0, end-input.size()));
		} else {
			out.addAll(input.subList(start, end));
		}
		return out;
	}

	static <X> Iterator<List<X>> symmetric(final List<X> input, int width) { 
		if (input.size() < width) throw new ArrayIndexOutOfBoundsException("window width larger than list size");
		return new Iterator<List<X>>() {

			int position = 0;
			

			@Override
			public boolean hasNext() {
				return position < input.size();
			}

			@Override
			public List<X> next() {
				if (!hasNext()) throw new NoSuchElementException();
				return symmetric(input,width,position++);
			}

		};
	}
	
	static <X> List<X> symmetric(final List<X> input, int width, int position) {
		if (position > input.size() || position < 0) throw new ArrayIndexOutOfBoundsException("window position outside of list");
		int i = width/2;
		int start = position-i;
		int end = position+i+1;
		List<X> out = new ArrayList<>();
		if (start < 0) {
			out.addAll(input.subList(1, 1-start));
			Collections.reverse(out);
			out.addAll(input.subList(0, end));
		} else if (end > input.size()) {
			List<X> tmp = new ArrayList<>(input.subList(2*input.size()-end-1, input.size()-1));
			Collections.reverse(tmp);
			out.addAll(input.subList(start, input.size()));
			out.addAll(tmp);
		} else {
			out.addAll(input.subList(start, end));
		}
		start = start == input.size() ? 0 : start+1;
		end = end == input.size() ? 0 : end+1;
		return out;
	}

	public static <X> Iterator<List<X>> tailed(final List<X> input, int width) {
		if (input.size() < width) throw new ArrayIndexOutOfBoundsException("window width larger than list size");
		return new Iterator<List<X>>() {

			int position = 0;
			

			@Override
			public boolean hasNext() {
				return position < input.size();
			}

			@Override
			public List<X> next() {
				if (!hasNext()) throw new NoSuchElementException();
				return tailed(input,width,position++);
			}

		};
	}
	
	static <X> List<X> tailed(final List<X> input, int width, int position) {
		if (position > input.size() || position < 0) throw new ArrayIndexOutOfBoundsException("window position outside of list");
		int i = width/2;
		int start = position-i;
		int end = position+i+1;
		List<X> out = new ArrayList<>();
		int tmp = start;
		while (tmp<0) {
			out.add(input.get(0)); 
			tmp++;
		}
		out.addAll(input.subList(tmp, end>input.size()-1 ? input.size()-1 : end));
		tmp = end;
		while (tmp > input.size()-1) {
			out.add(input.get(input.size()-1));
			tmp-=1;
		}
		start = start == input.size() ? 0 : start+1;
		end = end == input.size() ? 0 : end+1;
		return out;
	}

}