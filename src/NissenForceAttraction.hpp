
#ifndef NISSENFORCEATTRACTION_HPP_
#define NISSENFORCEATTRACTION_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class NissenForceAttraction : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mS;
        archive & mBeta;
        archive & mGrowthDuration;
    }

    double mS;
    double mBeta;
    double mGrowthDuration;

public:

    NissenForceAttraction();

    virtual ~NissenForceAttraction();

    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                           unsigned nodeBGlobalIndex,
                                                           AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    double GetS();
    void SetS(double s);

    double GetBeta();
    void SetBeta(double beta);
    
    double GetGrowthDuration();
    void SetGrowthDuration(double GrowthDuration);

    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NissenForceAttraction)

#endif /*NISSENFORCEATTRACTION_HPP_*/
