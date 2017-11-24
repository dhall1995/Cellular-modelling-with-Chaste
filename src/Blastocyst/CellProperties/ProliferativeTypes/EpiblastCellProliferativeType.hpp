
#ifndef EPIBLASTCELLPROLIFERATIVETYPE_HPP_
#define EPIBLASTCELLPROLIFERATIVETYPE_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProliferativeType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Epiblast cell ProliferativeType.
 *
 * Each Cell owns a CellProliferativeType, which may include a shared pointer
 * to an object of this type.
 *
 * The EpiblastCellProliferativeType object keeps track of the number of cells that are epiblast, as well
 * as what colour should be used by the visualizer to display cells that are apoptotic.
 */
class EpiblastCellProliferativeType : public AbstractCellProliferativeType
{
private:
    
   /** Needed for serialization. */
   friend class boost::serialization::access;
   /**
    * Archive the member variables.
    *
    * @param archive the archive
    * @param version the current version of this class
    */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProliferativeType>(*this);
    }

public:

    /**
     * Constructor.
     */
    EpiblastCellProliferativeType();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(EpiblastCellProliferativeType)

#endif /* EPIBLASTCELLPROLIFERATIVETYPE_HPP_ */
