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
	double m_Vroot[3];		// Molar volumes returned by cubic solver
	doublereal m_a, m_b;	// a_mix and b_mix

};

}


