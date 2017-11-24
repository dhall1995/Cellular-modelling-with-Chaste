
#ifndef TROPHECTODERMCELLPROLIFERATIVETYPE_HPP_
#define TROPHECTODERMCELLPROLIFERATIVETYPE_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProliferativeType.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Trophectoderm cell proliferative type.
 *
 * Subclass of AbstractCellProliferativeType defining a trophectoderm cell
 * what colour should be used by the visualizer to display cells that are apoptotic.
 */
class TrophectodermCellProliferativeType : public AbstractCellProliferativeType
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
    TrophectodermCellProliferativeType();
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TrophectodermCellProliferativeType)

#endif /* TROPHECTODERMCELLPROLIFERATIVETYPE_HPP_ */
