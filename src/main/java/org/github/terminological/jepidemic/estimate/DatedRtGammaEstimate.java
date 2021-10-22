package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Optional;

import org.github.terminological.jepidemic.gamma.GammaParameters;

public class DatedRtGammaEstimate extends GammaParameters implements Comparable<DatedRtGammaEstimate> {

	int tau;
	LocalDate date;
	double incidence;
	int profileId;
	DatedRtGammaEstimate prior = null;
	
	private DatedRtGammaEstimate(double shape, double scale, int tau, LocalDate date, double incidence, int profileId) {
		super(shape, scale);
		this.tau = tau;
		this.date = date;
		this.incidence = incidence;
		this.profileId = profileId;
	}
	
	public DatedRtGammaEstimate(GammaParameters p, int tau, LocalDate date, double incidence, int profileId) {
		super(p.getShape(), p.getScale());
		this.tau = tau;
		this.date = date;
		this.incidence = incidence;
		this.profileId = profileId;
	}

	@Override
	public int compareTo(DatedRtGammaEstimate o) {
		return this.date.compareTo(o.date)*1000+(this.tau-o.tau);
	}

	public int getWindow() {return tau+1;}

	public LocalDate getEffectiveDate() {
		if (tau > 0)
			return date.minusDays(tau/2);
		return date;
	}
	
	public LocalDate getStartDate() {
		if (tau > 0)
			return date.minusDays(tau);
		return date;
	}
	
	public LocalDate getEndDate() {
		return date;
	}
	
	public DatedRtGammaEstimate withPrior(DatedRtGammaEstimate prior) {
		this.prior = prior;
		return this;
	}
	
	public Optional<DatedRtGammaEstimate> getPrior() {
		return Optional.ofNullable(prior);
	}
	
	Optional<DatedRtGammaEstimate> wider(double factor) {
		GammaParameters out = this.convert().wider(factor).convert();
		if (!(Double.isFinite(out.getShape()) && Double.isFinite(out.getScale()))) return Optional.empty();
		return Optional.of(out.withDate(tau, date, incidence, profileId).withPrior(prior));
	}
	
	public int getProfileId() {
		return profileId;
	}
	
	public double getIncidence() {return incidence;}
	
	public String toString() {return this.getEffectiveDate().toString()+": "+super.toString();}

//	public CoriEstimationSummaryEntry toStatSummary() {
//		return new CoriEstimationSummaryEntry(this.getStartDate(), this.getWindow(), this.getIncidence());
//	}
//	
//	public CoriEstimationSummaryEntry toStatSummary(LocalDate d) {
//		return this.withDate(d, this.getWindow(), this.getIncidence());
//	}
}
