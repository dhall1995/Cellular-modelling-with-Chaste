
#include "PrECellProliferativeType.hpp"

PrECellProliferativeType::PrECellProliferativeType()
    : AbstractCellProliferativeType(4)
{
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PrECellProliferativeType)
