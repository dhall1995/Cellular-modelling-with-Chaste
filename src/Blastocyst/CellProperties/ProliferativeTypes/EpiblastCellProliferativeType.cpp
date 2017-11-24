
#include "EpiblastCellProliferativeType.hpp"

EpiblastCellProliferativeType::EpiblastCellProliferativeType()
    : AbstractCellProliferativeType(3)
{
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(EpiblastCellProliferativeType)
