
#ifndef PRECELLPROLIFERATIVETYPE_HPP_
#define PRECELLPROLIFERATIVETYPE_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProliferativeType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * PrE cell proliferative type.
 *
 * Each Cell owns a shared pointer
 * to an object of this type.
 *
 * The PrECellProliferative type object keeps track of the number of cells that are Primitive Endoderm, as well
 * as what colour should be used by the visualizer to display cells that are PrE.
 */
class PrECellProliferativeType : public AbstractCellProliferativeType
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
    PrECellProliferativeType();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PrECellProliferativeType)

#endif /* PRECELLPROLIFERATIVETYPE_HPP_ */
