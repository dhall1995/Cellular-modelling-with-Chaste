
#ifndef PRECOMPACTIONCELLCYCLEMODEL_HPP_
#define PRECOMPACTIONCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

class PreCompactionCellCycleModel : public AbstractSimpleCellCycleModel
{
    friend class TestSimpleCellCycleModels;
    
private:
    
    double mMinCellCycleDuration;
    
    double mMaxCellCycleDuration;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        
        //Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mMinCellCycleDuration;
        archive & mMaxCellCycleDuration;
    }
    
protected:
    
    PreCompactionCellCycleModel(const PreCompactionCellCycleModel& rModel);
    
public:
    
    PreCompactionCellCycleModel();
    
    void SetCellCycleDuration();
    
    AbstractCellCycleModel* CreateCellCycleModel();
    
    double GetMinCellCycleDuration();
    
    void SetMinCellCycleDuration(double minCellCycleDuration);
    
    double GetMaxCellCycleDuration();
    
    void SetMaxCellCycleDuration(double maxCellCycleDuration);
    
    double GetAverageTransitCellCycleTime();
    
    double GetAverageStemCellCycleTime();
    
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PreCompactionCellCycleModel)

#endif /* PRECOMPACTIONCELLCYCLEMODEL_HPP_ */


