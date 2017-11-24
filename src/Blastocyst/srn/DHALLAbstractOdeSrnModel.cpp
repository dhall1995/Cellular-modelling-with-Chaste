
#include "DHALLAbstractOdeSrnModel.hpp"
#include "Debug.hpp"

DHALLAbstractOdeSrnModel::DHALLAbstractOdeSrnModel(unsigned stateSize, boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractSrnModel(),
      CellCycleModelOdeHandler(SimulationTime::Instance()->GetTime(), pOdeSolver),
      mStateSize(stateSize)
{
}

DHALLAbstractOdeSrnModel::DHALLAbstractOdeSrnModel(const DHALLAbstractOdeSrnModel& rModel)
    : AbstractSrnModel(rModel),
      CellCycleModelOdeHandler(rModel),
      mInitialConditions(rModel.mInitialConditions),
      mStateSize(rModel.mStateSize)
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
}

DHALLAbstractOdeSrnModel::~DHALLAbstractOdeSrnModel()
{
}

void DHALLAbstractOdeSrnModel::SimulateToCurrentTime()
{
    TRACE("Now attempting simulate to current time within DHALLAbstractOdeSrnModel");
    assert(mpOdeSystem != nullptr);
    assert(SimulationTime::Instance()->IsStartTimeSetUp());

    double current_time = SimulationTime::Instance()->GetTime();

    // Run ODEs if needed
    if (current_time > mLastTime)
    {
        if (!this->mFinishedRunningOdes)
        {
            // Update whether a stopping event has occurred
            this->mFinishedRunningOdes = SolveOdeToTime(current_time);
        }
        else
        {
            // ODE model finished, just increasing time...
        }
    }

    // Update the SimulatedToTime value
    mLastTime = current_time;
    SetSimulatedToTime(current_time);
}

void DHALLAbstractOdeSrnModel::Initialise(AbstractOdeSystem* pOdeSystem)
{
    TRACE("Now within AbstractOdeSrnModel::Initialise");
    assert(mpOdeSystem == nullptr);
    assert(mpCell != nullptr);

    mpOdeSystem = pOdeSystem;
    TRACE("ODE system has been initialised");
    if (mInitialConditions == std::vector<double>())
    {
        mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());
        TRACE("Initial conditions have been set");
    }
    else
    {
        mpOdeSystem->SetStateVariables(mInitialConditions);
        TRACE("Initial conditions have been set");
    }

    SetLastTime(mSimulatedToTime);
}

void DHALLAbstractOdeSrnModel::ResetForDivision()
{
    AbstractSrnModel::ResetForDivision();
    assert(mLastTime == mSimulatedToTime);
    mFinishedRunningOdes = false;
}

void DHALLAbstractOdeSrnModel::SetInitialConditions(std::vector<double> initialConditions)
{
    assert(initialConditions.size() == mStateSize);
    mInitialConditions = initialConditions;
}

void DHALLAbstractOdeSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSrnModel::OutputSrnModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DHALLAbstractOdeSrnModel)
