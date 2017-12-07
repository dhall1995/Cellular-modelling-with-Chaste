
#include "CellPolaritySrnModel.hpp"
#include "Debug.hpp"

#include <cassert>

CellPolaritySrnModel::CellPolaritySrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : DHALLAbstractOdeSrnModel(1, pOdeSolver)
{
//    TRACE("Now attempting to initialise the Srn Model");
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<CellPolaritySrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<CellPolaritySrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

CellPolaritySrnModel::CellPolaritySrnModel(const CellPolaritySrnModel& rModel)
    : DHALLAbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
//    	TRACE("Now constructing srn Model. We have the following state variables for our ODE");
//    	PRINT_VECTOR(rModel.GetOdeSystem()->rGetStateVariables());
    	assert(rModel.GetOdeSystem());
    	SetOdeSystem(new CellPolarityOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* CellPolaritySrnModel::CreateSrnModel()
{
    return new CellPolaritySrnModel(*this);
}

void CellPolaritySrnModel::SimulateToCurrentTime()
{
//    TRACE("Now attempting SimulateToCurrentTime within CellPolaritySrnModel");

	// Custom behaviour
    UpdatedVpdAlpha();

    // Run the ODE simulation as needed
    DHALLAbstractOdeSrnModel::SimulateToCurrentTime();
}

void CellPolaritySrnModel::Initialise()
{
//    TRACE("Now attempting CellPolaritySrnModel::initialise within CellPolaritySrnModel");
	
	DHALLAbstractOdeSrnModel::Initialise(new CellPolarityOdeSystem);
}

void CellPolaritySrnModel::UpdatedVpdAlpha()
{
//    TRACE("Now attempting UpdatePolarityAngle within CellPolaritySrnModel");
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);

    double dVpdAlpha = mpCell->GetCellData()->GetItem("dVpdAlpha");
    
    mpOdeSystem->SetParameter("dVpdAlpha", dVpdAlpha);
}

double CellPolaritySrnModel::GetPolarityAngle()
{
//    TRACE("Now attempting GetPolarityAngle within CellPolaritySrnModel");
	assert(mpOdeSystem != NULL);
//	TRACE("OdeSystem Exists");
    double polarity_angle = mpOdeSystem->rGetStateVariables()[0];
//    PRINT_VARIABLE(polarity_angle);
    return polarity_angle;
}

void CellPolaritySrnModel::SetPolarityAngle(double polarityAngle)
{
    assert(mpOdeSystem != NULL);
    mpOdeSystem->rGetStateVariables()[0] = polarityAngle;
}

double CellPolaritySrnModel::GetdVpdAlpha()
{
//    TRACE("Now attempting GetdVdAlpha within CellPolaritySrnModel");
	assert(mpOdeSystem != NULL);
    double dVpdAlpha = mpOdeSystem->GetParameter("dVpdAlpha");
    return dVpdAlpha;
}

void CellPolaritySrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    DHALLAbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellPolaritySrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CellPolaritySrnModel)
