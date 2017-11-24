
#ifndef POLARITYCELLPROPERTY_HPP_
#define POLARITYCELLPROPERTY_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Polarity cell property.
 *
 * Each Cell owns a CellPropertyCollection, which may include a shared pointer
 * to an object of this type.
 *
 * The PolarityCellProperty object keeps track of the number of cells that are polar, as well
 * as what colour should be used by the visualizer to display cells that are apoptotic.
 */
class PolarityCellProperty : public AbstractCellProperty
{
private:

   /**
    * Colour for use by visualizer.
    */
   unsigned mColour;

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
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

    public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this property should be in the visualizer (defaults to 5)
     */
    PolarityCellProperty(unsigned colour=5);

    /**
     * Destructor.
     */
    virtual ~PolarityCellProperty();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PolarityCellProperty)

#endif /* POLARITYCELLPROPERTY_HPP_ */
