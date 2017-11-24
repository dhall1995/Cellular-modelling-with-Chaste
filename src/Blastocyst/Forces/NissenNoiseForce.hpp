
#ifndef NISSENNOISEFORCE_HPP_
#define NISSENNOISEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A 'diffusion force' class to model the random movement of nodes.
 *
 * This class works with all off-lattice cell populations.
 */
template<unsigned DIM>
class NissenNoiseForce : public AbstractForce<DIM>
{
private :

    /**
     * Standard Deviation governing noise parameter.
     */
    double mNoiseStandardDev;

    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mNoiseStandardDev;
    }

public :

    /**
     * Constructor.
     */
    NissenNoiseForce();

    /**
     * Destructor.
     */
    ~NissenNoiseForce();

    /**
     * Set the Noise Standard Deviation which sets the Noise in the system
     */
    void SetNoiseStandardDev(double newValue);

    /**
     * Get Noise Standard Deviation.
     *
     * @return mNoiseStandardDev
     */
    double GetNoiseStandardDev();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NissenNoiseForce)

#endif /*NISSENNOISEFORCE_HPP_*/
