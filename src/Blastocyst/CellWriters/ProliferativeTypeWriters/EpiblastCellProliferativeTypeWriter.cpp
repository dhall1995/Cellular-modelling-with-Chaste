
#include "EpiblastCellProliferativeTypeWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "EpiblastCellProliferativeType.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
EpiblastCellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::EpiblastCellProliferativeTypeWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizlabels")
{
    this->mVtkCellDataName = "Epiblast Cells";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double EpiblastCellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double label = 0.0;
    if (pCell->HasCellProperty<EpiblastCellProliferativeType>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<EpiblastCellProliferativeType>();
        boost::shared_ptr<EpiblastCellProliferativeType> p_label = boost::static_pointer_cast<EpiblastCellProliferativeType>(collection.GetProperty());
        label = p_label->GetColour();
    }
    return label;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EpiblastCellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned label = 0;
    if (pCell->HasCellProperty<EpiblastCellProliferativeType>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<EpiblastCellProliferativeType>();
        boost::shared_ptr<EpiblastCellProliferativeType> p_label = boost::static_pointer_cast<EpiblastCellProliferativeType>(collection.GetProperty());
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
template class EpiblastCellProliferativeTypeWriter<1,1>;
template class EpiblastCellProliferativeTypeWriter<1,2>;
template class EpiblastCellProliferativeTypeWriter<2,2>;
template class EpiblastCellProliferativeTypeWriter<1,3>;
template class EpiblastCellProliferativeTypeWriter<2,3>;
template class EpiblastCellProliferativeTypeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(EpiblastCellProliferativeTypeWriter)
