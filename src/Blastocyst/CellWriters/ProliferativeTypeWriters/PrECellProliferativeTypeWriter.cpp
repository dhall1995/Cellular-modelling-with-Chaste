
#include "PrECellProliferativeTypeWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "PrECellProliferativeType.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PrECellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::PrECellProliferativeTypeWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizlabels")
{
    this->mVtkCellDataName = "PrE Cells";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double PrECellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double label = 0.0;
    if (pCell->HasCellProperty<PrECellProliferativeType>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<PrECellProliferativeType>();
        boost::shared_ptr<PrECellProliferativeType> p_label = boost::static_pointer_cast<PrECellProliferativeType>(collection.GetProperty());
        label = p_label->GetColour();
    }
    return label;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PrECellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned label = 0;
    if (pCell->HasCellProperty<PrECellProliferativeType>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<PrECellProliferativeType>();
        boost::shared_ptr<PrECellProliferativeType> p_label = boost::static_pointer_cast<PrECellProliferativeType>(collection.GetProperty());
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
template class PrECellProliferativeTypeWriter<1,1>;
template class PrECellProliferativeTypeWriter<1,2>;
template class PrECellProliferativeTypeWriter<2,2>;
template class PrECellProliferativeTypeWriter<1,3>;
template class PrECellProliferativeTypeWriter<2,3>;
template class PrECellProliferativeTypeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PrECellProliferativeTypeWriter)
