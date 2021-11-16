package org.github.terminological.jepidemic.growth;

import org.github.terminological.jepidemic.InfectivityProfile;
import org.github.terminological.jepidemic.distributions.ExtendedGammaDistribution;

public abstract class EstimateMetadata {

	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + tau;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof EstimateMetadata))
			return false;
		EstimateMetadata other = (EstimateMetadata) obj;
		if (tau != other.tau)
			return false;
		return true;
	}

	int tau;
	
	public EstimateMetadata(int tau) {
		this.tau = tau;
	}
	
	public int getTau() {return tau;}
	public abstract boolean isTwoSided();
	
	// ================================================
	
	
	public static class RtMetadata extends EstimateMetadata {

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = super.hashCode();
			result = prime * result + ((infectivityProfile == null) ? 0 : infectivityProfile.hashCode());
			result = prime * result + ((priorR0 == null) ? 0 : priorR0.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (!super.equals(obj))
				return false;
			if (!(obj instanceof RtMetadata))
				return false;
			RtMetadata other = (RtMetadata) obj;
			if (infectivityProfile == null) {
				if (other.infectivityProfile != null)
					return false;
			} else if (!infectivityProfile.equals(other.infectivityProfile))
				return false;
			if (priorR0 == null) {
				if (other.priorR0 != null)
					return false;
			} else if (!priorR0.equals(other.priorR0))
				return false;
			return true;
		}

		InfectivityProfile infectivityProfile;
		ExtendedGammaDistribution priorR0;
		
		public RtMetadata(int tau, InfectivityProfile infectivityProfile, ExtendedGammaDistribution priorR0) {
			super(tau);
			this.infectivityProfile = infectivityProfile;
			this.priorR0 = priorR0;
		}
		
		public InfectivityProfile getInfectivityProfile() {return infectivityProfile;}

		public ExtendedGammaDistribution defaultR0Prior() {return priorR0;}

		@Override
		public boolean isTwoSided() {
			return false;
		}

	}
	
	// ================================================
	
	public static class GrowthMetadata extends EstimateMetadata {

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = super.hashCode();
			result = prime * result + ((priorLambda == null) ? 0 : priorLambda.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (!super.equals(obj))
				return false;
			if (!(obj instanceof GrowthMetadata))
				return false;
			GrowthMetadata other = (GrowthMetadata) obj;
			if (priorLambda == null) {
				if (other.priorLambda != null)
					return false;
			} else if (!priorLambda.equals(other.priorLambda))
				return false;
			return true;
		}

		ExtendedGammaDistribution priorLambda;
		
		public GrowthMetadata(GrowthMetadata copy) {
			this(copy.tau, copy.priorLambda);
		}
		
		public GrowthMetadata(int tau, ExtendedGammaDistribution priorLambda) {
			super(tau);
			this.priorLambda = priorLambda;
		}

		public ExtendedGammaDistribution defaultLambdaPrior() {return priorLambda;}

		public RtFromGrowthMetadata extend(InfectivityProfile infProf) {
			return new RtFromGrowthMetadata(this,infProf);
		}

		@Override
		public boolean isTwoSided() {
			return true;
		}
		
//		public GrowthRate getUnknownGR(GrowthRateTimeseriesEntry forEntry) {
//			return new GrowthRate((GrowthRateDistribution) null,forEntry,this);
//		}

	}

	// ================================================
	
	
	public static class RtFromGrowthMetadata extends GrowthMetadata {

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = super.hashCode();
			result = prime * result + ((infectivityProfile == null) ? 0 : infectivityProfile.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (!super.equals(obj))
				return false;
			if (!(obj instanceof RtFromGrowthMetadata))
				return false;
			RtFromGrowthMetadata other = (RtFromGrowthMetadata) obj;
			if (infectivityProfile == null) {
				if (other.infectivityProfile != null)
					return false;
			} else if (!infectivityProfile.equals(other.infectivityProfile))
				return false;
			return true;
		}

		InfectivityProfile infectivityProfile;
		
		public RtFromGrowthMetadata(GrowthMetadata meta, InfectivityProfile infProf) {
			super(meta);
			this.infectivityProfile = infProf;
			
		}
		
		public InfectivityProfile getInfectivityProfile() {return infectivityProfile;}

		

	}

	

}
