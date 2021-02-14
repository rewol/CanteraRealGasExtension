#include "VTPengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
namespace Cantera
{
	VTPengRobinson::VTPengRobinson() : m_a(0), m_b(0)
	{
		std::fill_n(m_Vroot, 3, 0.0);
	}

	VTPengRobinson::VTPengRobinson(const std::string& infile, const std::string& id_) : m_a(0), m_b(0)
	{
		std::fill_n(m_Vroot, 3, 0.0);
		initThermoFile(infile, id_);
	}

	vector_fp VTPengRobinson::getCoeff(const std::string& iName)
	{
		vector_fp spCoeff{ 0, 0, 0, 0, 0 };
		// Get number of species in the database
		// open xml file critProperties.xml
		XML_Node* doc = get_XML_File("critProperties.xml");
		size_t nDatabase = doc->nChildren();
		// Loop through all species in the database and attempt to match supplied
		// species to each. If present, calculate pureFluidParameters a_k and b_k
		// based on crit properties T_c and P_c
		for (size_t isp = 0; isp < nDatabase; isp++) {
			XML_Node& acNodeDoc = doc->child(isp);
			std::string iNameLower = toLowerCopy(iName);
			std::string dbName = toLowerCopy(acNodeDoc.attrib("name"));
			// Attempt to match provided species iName to current database species
			//  dbName:
			if (iNameLower == dbName) {
				// Read from database and calculate a and b coefficients
				double vParams;
				double T_crit = 0.0, P_crit = 0.0, w_ac = 0.0, V_crit = 0.0, Z_crit = 0;
				if (acNodeDoc.hasChild("Tc")) {
					vParams = 0.0;
					XML_Node& xmlChildCoeff = acNodeDoc.child("Tc");
					if (xmlChildCoeff.hasAttrib("value")) {
						std::string critTemp = xmlChildCoeff.attrib("value");
						vParams = strSItoDbl(critTemp);
					}
					if (vParams <= 0.0) { //Assuming that Pc and Tc are non zero.
						throw CanteraError("PengRobinson::getCoeff", "Critical Temperature must be positive");
					}
					T_crit = vParams;
				}
				if (acNodeDoc.hasChild("Pc")) {
					vParams = 0.0;
					XML_Node& xmlChildCoeff = acNodeDoc.child("Pc");
					if (xmlChildCoeff.hasAttrib("value")) {
						std::string critPressure = xmlChildCoeff.attrib("value");
						vParams = strSItoDbl(critPressure);
					}
					if (vParams <= 0.0) { //Assuming that Pc and Tc are non zero.
						throw CanteraError("PengRobinson::getCoeff", "Critical Pressure must be positive");
					}
					P_crit = vParams;
				}
				if (acNodeDoc.hasChild("omega")) {
					vParams = 0.0;
					XML_Node& xmlChildCoeff = acNodeDoc.child("omega");
					if (xmlChildCoeff.hasChild("value")) {
						std::string acentric_factor = xmlChildCoeff.attrib("value");
						vParams = strSItoDbl(acentric_factor);
					}
					w_ac = vParams;
				}
				if (acNodeDoc.hasChild("Vc")) {
					vParams = 0.0;
					XML_Node& xmlChildCoeff = acNodeDoc.child("Vc");
					if (xmlChildCoeff.hasChild("value")) {
						std::string critVolume = xmlChildCoeff.attrib("value");
						vParams = strSItoDbl(critVolume);
					}
					V_crit = vParams;
				}
				if (acNodeDoc.hasChild("Zc")) {
					vParams = 0.0;
					XML_Node& xmlChildCoeff = acNodeDoc.child("Zc");
					if (xmlChildCoeff.hasChild("value")) {
						std::string critComp = xmlChildCoeff.attrib("value");
						vParams = strSItoDbl(critComp);
					}
					Z_crit = vParams;
				}

				spCoeff[0] = m_a0 * (GasConstant * GasConstant) * (T_crit * T_crit) / P_crit;
				spCoeff[1] = m_b0 * GasConstant * T_crit / P_crit;
				spCoeff[2] = w_ac;
				spCoeff[3] = V_crit;
				spCoeff[4] = Z_crit;
				break;
			}
		}
		return spCoeff;
	}

	void VTPengRobinson::setSpeciesCoeffs()
	{
		for (auto& item : speciesNames())
		{
			size_t k = speciesIndex(item);
			vector_fp coeffs = getCoeff(item);
			a_k[k] = coeffs.at(0);
			b_k[k] = coeffs.at(1);
			w_k[k] = coeffs.at(2);
			vc_k[k] = coeffs.at(3);
			zc_k[k] = coeffs.at(4);
		}
	}

	void VTPengRobinson::calculateKappa()
	{
		kappa_k.reserve(m_kk);
		for (size_t k = 0; k < m_kk; k++)
		{
			kappa_k.push_back(0.37464 + 1.54226 * w_k.at(k) - 0.26992 * w_k.at(k) * w_k.at(k));
		}
	}

	doublereal VTPengRobinson::calculateSpeciesCritTemp(size_t k)
	{
		return (a_k.at(k) * m_b0) / (b_k.at(k) * GasConstant * m_a0);
	}

	doublereal VTPengRobinson::calculateAlpha(const std::string& iName)
	{
		size_t k = speciesIndex(iName);
		doublereal tempcrit = calculateSpeciesCritTemp(k);
		doublereal temp = (1 + kappa_k.at(k) * (1 - pow((temperature() / tempcrit), 0.5)));
		return temp * temp;
	}

	void VTPengRobinson::calculateSpeciesAlpha()
	{
		for (size_t k = 0; k < m_kk; k++)
		{
			alpha_k[k] = a_k.at(k) * calculateAlpha(speciesName(k));
			//alpha_k.push_back(a_k.at(k) * calculateAlpha(speciesName(k)));
		}
	}

	void VTPengRobinson::calculateCrossSpecies()
	{
		for (size_t j = 0; j < m_kk; j++)
		{
			for (size_t i = 0; i < m_kk; i++)
			{
				a_alpha(i, j) = pow(alpha_k.at(i) * alpha_k.at(j), 0.5);
			}
		}
	}

	void VTPengRobinson::calculateAB()
	{
		m_a = 0;
		m_b = 0;
		
		for (size_t i = 0; i < m_kk; i++)
		{
			m_b = m_b + moleFractions_.at(i) * b_k.at(i);
			for (size_t j = 0; j < m_kk; j++)
			{
				m_a = m_a + moleFractions_.at(i) * moleFractions_.at(j) * a_alpha(j, i);
			}
		}
	}

	double VTPengRobinson::GibbsFreeEnergyChange(double* Z, doublereal A, doublereal B, bool* phasecheck)
	{
		double maximum = 0;
		double minimum = 0;
		double final = 0;
		if (Z[0] > Z[1])
		{
			maximum = Z[0];
			minimum = Z[1];
		}
		else
		{
			maximum = Z[1];
			minimum = Z[0];
		}


		if (Z[2] > maximum)
			maximum = Z[2];
		else if (Z[2] < minimum)
			minimum = Z[2];

		/*std::cout << std::endl;
		std::cout << "Liquid-like Z is = " << minimum << std::endl;
		std::cout << "Vapor-like Z is  = " << maximum << std::endl;
		std::cout << std::endl;*/

		// Check Gibbs Free Energy change:
		double dg = 0;
		double s1 = 1 + pow(2, 0.5);
		double s2 = 1 - pow(2, 0.5);
		double zh = maximum;
		double zl = minimum;

		double term1 = log((zl - B) / (zh - B));
		double term2 = A / (B * (s2 - s1));
		double t1 = zl + s1 * B;
		double t2 = zl + s2 * B;
		double t3 = zh + s1 * B;
		double t4 = zh + s2 * B;
		double term3 = log((t1 / t2) * (t4 / t3));

		dg = (zh - zl) + term1 - term2 * term3;

		// For positive dg, phase is liquid-like. For negative dg, phase is vapor like.
		//if (dg > 0)
		//{
		//	//std::cout << "Phase is liquid-like !" << std::endl;
		//	*phasecheck = true;
		//}			//std::cout << "Phase is liquid-like !" << std::endl;
		//else
		//{
		//	//std::cout << "Phase is vapor-like" << std::endl;
		//	*phasecheck = false;
		//}

		*phasecheck = (dg > 0) ? true : false;
			

		final = (dg > 0) ? zl : zh;
		return final;
	}

	int VTPengRobinson::deitersSolver(double temp, double pressure, doublereal a, doublereal b, double* Vroot, double tcrit, double vcrit)
	{
		bool isliquid = false;
		// Pressure should be in Pascal
		double Z[3] = { 0, 0, 0 };
		double m[4];
		double R = GasConstant;
		double A = (a * pressure) / (R * R * temp * temp);
		double B = (b * pressure) / (R * temp);
		m[0] = -1 * (A * B - pow(B, 2) - pow(B, 3));
		m[1] = A - 3 * pow(B, 2) - 2 * B;
		m[2] = -1 * (1 - B);
		m[3] = 1;
		double x_infl = -0.3333 * m[2];
		double y = m[0] + x_infl * (m[1] + x_infl * (m[2] + x_infl));
		double x = 0;
		double c1 = 0;
		double c0 = 0;
		int nor = 0;		// Number of roots
		double finalroot = 0;

		if (y == 0)
		{
			Z[0] = x_infl;
			c1 = Z[0] + m[2];
			c0 = c1 * Z[0] + m[1];
			double delta = c1 * c1 * 0.25 - c0;
			double beta = -0.5 * c1;
			if (delta <= 0)
			{
				//std::cout << "other roots cant be used, there is only one root" << std::endl;
				//std::cout << "Z[0] = " << Z[0] << std::endl;
				nor = 1;
				finalroot = Z[0];
			}
			else
			{
				//std::cout << "there are three usable roots" << std::endl;
				nor = 2;
				Z[1] = -0.5 * c1 - pow(c1 * c1 * 0.25 - c0, 0.5);
				Z[2] = -0.5 * c1 + pow(c1 * c1 * 0.25 - c0, 0.5);
				//for (int k = 0; k < 3; k++)
					//std::cout << "Z[" << k << "] = " << Z[k] << std::endl;
				finalroot = GibbsFreeEnergyChange(Z, A, B, &isliquid);
				//cout << "Final root is = " << finalroot << std::endl;

			}
		}

		double D = pow(m[2], 2) - 3 * m[1];
		if (D == 0)
		{
			Z[0] = x_infl - pow(y, 0.333);
			finalroot = Z[0];
			nor = 1;
		}

		else
		{
			bool cnd = (y > 0);
			if (D > 0)
			{
				x = x_infl + ((cnd) ? 0.66666 : -0.66666) * pow(D, 0.5);
			}
			else if (D < 0)
			{
				x = x_infl;
			}
			// Iteration
			double errorlimit = 1e-6;
			double dx;
			int iter = 0;
			double y1, y2, y3;
			do
			{
				y1 = m[0] + x * (m[1] + x * (m[2] + x));
				y2 = m[1] + 2 * x * m[2] + 3 * x * x;
				y3 = 2 * m[2] + 6 * x;
				dx = (y1 * y2) / (y2 * y2 - 0.5 * y1 * y3);
				x = x - dx;
				iter += 1;
				//std::cout << "Iteration No: " << iter << std::endl;

			} while (fabs(dx) > errorlimit * fabs(x));

			if (D > 0)
			{
				Z[0] = x;
				c1 = Z[0] + m[2];
				c0 = c1 * Z[0] + m[1];
				double delta = c1 * c1 * 0.25 - c0;
				double beta = -0.5 * c1;
				if (delta <= 0)
				{
					//std::cout << "other roots cant be used, there is only one root" << endl;
					//std::cout << "Z[0] = " << Z[0] << std::endl;
					finalroot = Z[0];
					nor = 1;
				}
				else
				{
					//std::cout << "there are three usable roots" << std::endl;
					nor = 2;
					Z[1] = -0.5 * c1 - pow(c1 * c1 * 0.25 - c0, 0.5);
					Z[2] = -0.5 * c1 + pow(c1 * c1 * 0.25 - c0, 0.5);
					/*for (int k = 0; k < 3; k++)
						cout << "Z[" << k << "] = " << Z[k] << std::endl;*/
					finalroot = GibbsFreeEnergyChange(Z, A, B, &isliquid);
					//cout << "Final root is = " << finalroot << std::endl;
				}
			}
			else
			{
				Z[0] = x;
				//std::cout << "There is only one usable root" << std::endl;
				//std::cout << "Z[0] = " << Z[0] << std::endl;
				finalroot = Z[0];
				nor = 1;
			}

		}
		m_Vroot[0] = (finalroot * R * temp / pressure);

		if (nor == 1)
		{
			if (temp < tcrit)
			{
				if (m_Vroot[0] < vcrit)
					nor = -1; // Liquid phase
			}
		}
		else
		{
			nor = (isliquid == true) ? -2 : 2; // Liquid phase is more stable where 2 phases are found
		}

		return nor;
	}

	double VTPengRobinson::densityCalc(double temp, double pressure, int phaseRequested, double rhoguess)
	{
		double density = 0;
		setTemperature(temp);
		double tcrit = critTemperature();
		double vcrit = critVolume();
		double mmw = meanMolecularWeight();
		double nor = deitersSolver(temp, pressure, m_a, m_b, m_Vroot, tcrit, vcrit);
		double vut = m_Vroot[0];
		if (nor < 0)
		{
			double vut = m_Vroot[0];
			density = mmw / volumeTranslation(vut);
		}
		else
		{
			density = mmw / vut;
		}

		return density;
	}

	void VTPengRobinson::updateMixingExpressions()
	{
		calculateSpeciesAlpha();
		calculateCrossSpecies();
		calculateAB();
	}

	double VTPengRobinson::volumeTranslation(double& vut)
	{
		double tcm, pcm, vcm, wm;
		pseudoCritProperties(tcm, pcm, vcm, wm);
		double cm = 0;
		for (size_t k = 0; k < m_kk; k++)
		{
			cm = cm * moleFractions_.at(k) * (0.4266 * zc_k.at(k) - 0.1101);
		}
		double vcmpr = (GasConstant * tcm / pcm) * m_vc - vcm;		
		
		
		
		return 0;
	}

	void VTPengRobinson::pseudoCritProperties(double& tcrit, double& pcrit, double& vcrit, double& wm)
	{
		vcrit = 0;
		wm = 0;
		double sum = 0; // This will be used to calculate tcrit !
		// Calculate pseudo-crit volume
		for (size_t k = 0; k < m_kk; k++)
		{
			vcrit = vcrit + moleFractions_.at(k) * vc_k.at(k);
			sum = sum + moleFractions_.at(k) * pow(vc_k.at(k), 0.666);
			wm = wm + moleFractions_.at(k) * w_k.at(k);
		}

		tcrit = 0;
		for (size_t k = 0; k < m_kk; k++)
		{
			tcrit = tcrit + moleFractions_.at(k) * pow(vc_k.at(k), 0.666) * calculateSpeciesCritTemp(k);
		}
		tcrit = tcrit / sum;

		pcrit = (0.2905 - wm * 0.085) * GasConstant * tcrit / vcrit;
	}

	

	void VTPengRobinson::setTemperature(const doublereal temp)
	{
		Phase::setTemperature(temp);
		_updateReferenceStateThermo();
		updateMixingExpressions();
		//iState_ = phaseState(true);
	}

	void VTPengRobinson::calcCriticalConditions(double& tcrit, double& pcrit, double& vcrit)
	{
		tcrit = m_a * m_b0 / (m_b * m_a0 * GasConstant);
		pcrit = m_b0 * GasConstant * tcrit / m_b;
		vcrit = m_vc * GasConstant * tcrit / pcrit;
	}

	double VTPengRobinson::critTemperature()
	{
		double tc, pc, vc;
		calcCriticalConditions(tc, pc, vc);
		return tc;
	}

	double VTPengRobinson::critPressure()
	{
		double tc, pc, vc;
		calcCriticalConditions(tc, pc, vc);
		return pc;
	}

	double VTPengRobinson::critVolume()
	{
		double tc, pc, vc;
		calcCriticalConditions(tc, pc, vc);
		return vc;
	}





	
}
