
#ifndef EPIBLASTCELLCYCLEMODEL_HPP_
#define UNDIFFERENTATEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

class EpiblastCellCycleModel : public AbstractSimpleCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:
    
    /**
     * The minimum cell cycle duration. Used to define the uniform distribution.
     * Defaults to 16 hours.
     */
    double mMinCellCycleDuration;


    /**
     * The maximum cell cycle duration. Used to define the uniform distribution.
     * Defaults to 18 hours.
     */
    double mMaxCellCycleDuration;
    
    //Needed for Serialization
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
    
    EpiblastCellCycleModel(const EpiblastCellCycleModel& rModel);
    
public:
    
    EpiblastCellCycleModel();
    
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
CHASTE_CLASS_EXPORT(EpiblastCellCycleModel)

#endif /* EPIBLASTCELLCYCLEMODEL_HPP_ */
