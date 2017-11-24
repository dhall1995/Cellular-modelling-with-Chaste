
#include "PolarityCellProperty.hpp"

PolarityCellProperty::PolarityCellProperty(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

PolarityCellProperty::~PolarityCellProperty()
{
}

unsigned PolarityCellProperty::GetColour() const
{
	return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(PolarityCellProperty)
