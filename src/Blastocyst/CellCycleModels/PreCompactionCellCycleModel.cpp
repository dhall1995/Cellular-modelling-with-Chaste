
#include "PreCompactionCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TrophectodermCellProliferativeType.hpp"

PreCompactionCellCycleModel::PreCompactionCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDuration(18.0), // Hours
      mMaxCellCycleDuration(22.0)  // Hours
{
}

PreCompactionCellCycleModel::PreCompactionCellCycleModel(const PreCompactionCellCycleModel& rModel)
: AbstractSimpleCellCycleModel(rModel),
mMinCellCycleDuration(rModel.mMinCellCycleDuration),
mMaxCellCycleDuration(rModel.mMaxCellCycleDuration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* PreCompactionCellCycleModel::CreateCellCycleModel()
{
    return new PreCompactionCellCycleModel(*this);
}

void PreCompactionCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    //Check that the cell actually exists
    assert(mpCell != NULL);
    
    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()) 
    {
        mCellCycleDuration = DBL_MAX;
	//Differentiated Cells shouldn't divide
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TrophectodermCellProliferativeType>())
    {
	mCellCycleDuration = 0.5*(mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf());
	//Trophectoderm Cell Cycle time should be half that of ICM
    }
    else
    {
        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf();
        // U[MinCCD,MaxCCD]
    }
}

double PreCompactionCellCycleModel::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void PreCompactionCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double PreCompactionCellCycleModel::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void PreCompactionCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double PreCompactionCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double PreCompactionCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

void PreCompactionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";
    
    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PreCompactionCellCycleModel)

