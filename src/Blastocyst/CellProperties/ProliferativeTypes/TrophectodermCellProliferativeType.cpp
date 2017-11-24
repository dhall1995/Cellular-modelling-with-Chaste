
#include "TrophectodermCellProliferativeType.hpp"

TrophectodermCellProliferativeType::TrophectodermCellProliferativeType()
    : AbstractCellProliferativeType(2)
{
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TrophectodermCellProliferativeType)
