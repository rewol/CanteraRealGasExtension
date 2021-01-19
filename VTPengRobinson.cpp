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
		vector_fp spCoeff{ NAN, NAN, NAN };
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
				double T_crit = 0.0, P_crit = 0.0, w_ac = 0.0;
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

				spCoeff[0] = m_a0 * (GasConstant * GasConstant) * (T_crit * T_crit) / P_crit;
				spCoeff[1] = m_b0 * GasConstant * T_crit / P_crit;
				spCoeff[2] = w_ac;
				break;
			}
		}
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
			
			alpha_k.push_back(a_k.at(k) * calculateAlpha(speciesName(k)));
		}
	}

	void VTPengRobinson::calculateCrossSpecies(Array2D& array2D)
	{
		for (size_t j = 0; j < m_kk; j++)
		{
			for (size_t i = 0; i < m_kk; i++)
			{
				array2D(i, j) = pow(alpha_k.at(i) * alpha_k.at(j), 0.5);
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
}