
#include "UndifferentiatedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

UndifferentiatedCellCycleModel::UndifferentiatedCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDuration(16.0), // Hours
      mMaxCellCycleDuration(18.0)  // Hours
{
}

UndifferentiatedCellCycleModel::UndifferentiatedCellCycleModel(const UndifferentiatedCellCycleModel& rModel)
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

AbstractCellCycleModel* UndifferentiatedCellCycleModel::CreateCellCycleModel()
{
    return new UndifferentiatedCellCycleModel(*this);
}

void UndifferentiatedCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    
    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
            mCellCycleDuration = DBL_MAX;
    }
    else
    {
        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf();
        // U[MinCCD,MaxCCD]
    }
}

double UndifferentiatedCellCycleModel::GetMinCellCycleDuration()
{
    return mMinCellCycleDuration;
}

void UndifferentiatedCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
    mMinCellCycleDuration = minCellCycleDuration;
}

double UndifferentiatedCellCycleModel::GetMaxCellCycleDuration()
{
    return mMaxCellCycleDuration;
}

void UndifferentiatedCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
    mMaxCellCycleDuration = maxCellCycleDuration;
}

double UndifferentiatedCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

double UndifferentiatedCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDuration + mMaxCellCycleDuration);
}

void UndifferentiatedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";
    
    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UndifferentiatedCellCycleModel)
