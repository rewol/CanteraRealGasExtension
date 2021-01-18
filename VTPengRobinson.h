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
private:
	vector_fp getCoeff(const std::string& iName);		// Returns a vector that contains a_k, b_k and w
	void setSpeciesCoeffs();							// Sets a_k and b_k


private:
	const doublereal m_a0 = 0.45724;
	const doublereal m_b0 = 0.0778;
	double m_Vroot[3];				// Molar volumes returned by cubic solver
	doublereal m_a, m_b;			// a_mix and b_mix
	vector_fp a_k, b_k, w_k;		// a_k values. This will only be calculated during initialization !

};

}


