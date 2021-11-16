package org.github.terminological.jepidemic.estimate;

import java.time.LocalDate;
import java.util.Optional;

public class DatedRtEstimate extends DatedGammaEstimate {
	
	int profileId;

	public DatedRtEstimate(
			GammaParameters p, 
			int tau, 
			LocalDate date, 
			double incidence, 
			double incidenceInWindow,
			int profileId
			) {
		super(p, tau, date, incidence, incidenceInWindow);
		this.profileId = profileId;
	}
	
	public DatedRtEstimate(DatedGammaEstimate g, int profId) {
		this(g, g.tau, g.date, g.incidence, g.incidenceInWindow, profId);
	}

	public int getProfileId() {
		return profileId;
	}
	
	Optional<DatedRtEstimate> wider(double factor) {
		return super.wider(factor)
			.map(d -> d.withProfileId(profileId));
	}

}
