package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Optional;

public class DatedGammaEstimate extends GammaParameters implements Comparable<DatedGammaEstimate> {

	int tau;
	LocalDate date;
	double incidence;
	double incidenceInWindow;
	DatedGammaEstimate prior = null;
	
	private DatedGammaEstimate(double shape, double scale, int tau, LocalDate date, double incidence, double incidenceInWindow) {
		super(shape, scale);
		this.tau = tau;
		this.date = date;
		this.incidence = incidence;
		this.incidenceInWindow = incidenceInWindow;
		
	}
	
	public DatedGammaEstimate(GammaParameters p, int tau, LocalDate date, double incidence, double incidenceInWindow) {
		super(p.getShape(), p.getScale());
		this.tau = tau;
		this.date = date;
		this.incidence = incidence;
		this.incidenceInWindow = incidenceInWindow;
		
	}

	@Override
	public int compareTo(DatedGammaEstimate o) {
		return this.date.compareTo(o.date)*1000+(this.tau-o.tau);
	}

	public int getWindow() {return tau+1;}

//	public LocalDate getEffectiveDate() {
////		if (tau > 0)
////			return date.minusDays(tau/2);
//		return date;
//	}
	
	public LocalDate getStartDate() {
		if (tau > 0)
			return date.minusDays(tau);
		return date;
	}
	
	public LocalDate getEndDate() {
		return date;
	}
	
	public DatedGammaEstimate withPrior(DatedGammaEstimate prior) {
		this.prior = prior;
		return this;
	}
	
	public Optional<DatedGammaEstimate> getPrior() {
		return Optional.ofNullable(prior);
	}
	
	Optional<? extends DatedGammaEstimate> wider(double factor) {
		GammaParameters out = this.convert().wider(factor).convert();
		if (!(Double.isFinite(out.getShape()) && Double.isFinite(out.getScale()))) return Optional.empty();
		return Optional.of(
			out
				.withDate(tau, date, incidence, incidenceInWindow)
				.withPrior(prior));
				//.withProfileId(profileId));
	}
	
	

	
	
	public double getIncidence() {return incidence;}
	
	public String toString() {return this.getEndDate().toString()+": "+super.toString();}

	public DatedRtEstimate withProfileId(int profId) {
		return new DatedRtEstimate(this, profId);
	}

//	public CoriEstimationSummaryEntry toStatSummary() {
//		return new CoriEstimationSummaryEntry(this.getStartDate(), this.getWindow(), this.getIncidence());
//	}
//	
//	public CoriEstimationSummaryEntry toStatSummary(LocalDate d) {
//		return this.withDate(d, this.getWindow(), this.getIncidence());
//	}
}
