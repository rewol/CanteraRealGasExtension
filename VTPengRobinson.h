#pragma once

#include "cantera/thermo/MixtureFugacityTP.h"
#include "cantera/base/Array.h"

namespace Cantera
{

	class VTPengRobinson : public MixtureFugacityTP
	{
	public:
		VTPengRobinson();
		VTPengRobinson(const std::string& infile, const std::string& id = "");
		int deitersSolver(double temp, double pressure, doublereal a, doublereal b, double Vroot);
		virtual double densityCalc(double temp, double pressure, int phaseRequested, double rhoguess);
	private:
		vector_fp getCoeff(const std::string& iName);			// Returns a vector that contains a_k, b_k and w
		void setSpeciesCoeffs();								// Sets a_k and b_k
		void calculateKappa();									// Calculates kappa for species using acentric factor
		doublereal calculateAlpha(const std::string& iName);	// Return alpha_k
		doublereal calculateSpeciesCritTemp(size_t k);			// Returns species critical temperature
		void calculateSpeciesAlpha();							// Calculates a_k * alpha
		void calculateCrossSpecies();					// Calculates sqrt(a_i * a_j)
		void calculateAB();										// Calculates a_m and b_m
		double GibbsFreeEnergyChange(double* Z, doublereal A, doublereal B);	// Calculates gibbs free energy change to determine phase
		virtual void updateMixingExpressions();


	private:
		const doublereal m_a0 = 0.45724;
		const doublereal m_b0 = 0.0778;
		double m_Vroot[3];				// Molar volumes returned by cubic solver
		doublereal m_a, m_b;			// a_mix and b_mix
		vector_fp a_k, b_k, w_k;		// species values. This will only be calculated during initialization !
		vector_fp kappa_k;				// kappa depends on w. This will only be calculated during initialization !

	private:
		vector_fp alpha_k;				// DO NOT FORGET TO RESERVE DURING CONSTRUCTION
		Array2D a_alpha;				// Contains cross species coefficients as per van der Waals mixture rules
	};

}


