
#include "CellPolarityOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"
#include "RandomNumberGenerator.hpp"

CellPolarityOdeSystem::CellPolarityOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<CellPolarityOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Polarity Angle for this cell
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 0.0); // soon overwritten

    this->mParameters.push_back(0.0);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

CellPolarityOdeSystem::~CellPolarityOdeSystem()
{
}

void CellPolarityOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    std::vector<double> Alpha = rY;
    double dVpdAlpha = this->mParameters[0]; // Shorthand for "this->mParameter("dVpdAlpha");"
    
    double x = RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, (10.0^(-3.0)));
    // The next line define the ODE system by Nissen et al.
    rDY[0] = -0.1*dVpdAlpha + x;  // d[V_i]/dAlpha_i
}

template<>
void CellwiseOdeSystemInformation<CellPolarityOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Polarity Angle");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    // If this is ever not the first parameter change the line
    // double mean_delta = this->mParameters[0]; in EvaluateYDerivatives().
    this->mParameterNames.push_back("dVpdAlpha");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellPolarityOdeSystem)
