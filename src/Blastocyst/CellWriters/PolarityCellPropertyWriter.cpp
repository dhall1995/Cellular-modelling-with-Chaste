
#include "PolarityCellPropertyWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "PolarityCellProperty.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PolarityCellPropertyWriter<ELEMENT_DIM, SPACE_DIM>::PolarityCellPropertyWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizlabels")
{
    this->mVtkCellDataName = "Polarity Cells";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PolarityCellPropertyWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double label = 0.0;
    if (pCell->HasCellProperty<PolarityCellProperty>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<PolarityCellProperty>();
        boost::shared_ptr<PolarityCellProperty> p_label = boost::static_pointer_cast<PolarityCellProperty>(collection.GetProperty());
        label = p_label->GetColour();
    }
    return label;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolarityCellPropertyWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned label = 0;
    if (pCell->HasCellProperty<PolarityCellProperty>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<PolarityCellProperty>();
        boost::shared_ptr<PolarityCellProperty> p_label = boost::static_pointer_cast<PolarityCellProperty>(collection.GetProperty());
        label = p_label->GetColour();
    }

    *this->mpOutStream << " " << label;

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << location_index;

    c_vector<double, SPACE_DIM> coords = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << " " << coords[i];
    }
}

// Explicit instantiation
template class PolarityCellPropertyWriter<1,1>;
template class PolarityCellPropertyWriter<1,2>;
template class PolarityCellPropertyWriter<2,2>;
template class PolarityCellPropertyWriter<1,3>;
template class PolarityCellPropertyWriter<2,3>;
template class PolarityCellPropertyWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PolarityCellPropertyWriter)
