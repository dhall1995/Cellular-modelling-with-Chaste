
#include "TrophectodermCellProliferativeTypeWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "TrophectodermCellProliferativeType.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
TrophectodermCellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::TrophectodermCellProliferativeTypeWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizlabels")
{
    this->mVtkCellDataName = "Trophectoderm Cells";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double TrophectodermCellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double label = 0.0;
    if (pCell->HasCellProperty<TrophectodermCellProliferativeType>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<TrophectodermCellProliferativeType>();
        boost::shared_ptr<TrophectodermCellProliferativeType> p_label = boost::static_pointer_cast<TrophectodermCellProliferativeType>(collection.GetProperty());
        label = p_label->GetColour();
    }
    return label;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void TrophectodermCellProliferativeTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned label = 0;
    if (pCell->HasCellProperty<TrophectodermCellProliferativeType>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<TrophectodermCellProliferativeType>();
        boost::shared_ptr<TrophectodermCellProliferativeType> p_label = boost::static_pointer_cast<TrophectodermCellProliferativeType>(collection.GetProperty());
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
template class TrophectodermCellProliferativeTypeWriter<1,1>;
template class TrophectodermCellProliferativeTypeWriter<1,2>;
template class TrophectodermCellProliferativeTypeWriter<2,2>;
template class TrophectodermCellProliferativeTypeWriter<1,3>;
template class TrophectodermCellProliferativeTypeWriter<2,3>;
template class TrophectodermCellProliferativeTypeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TrophectodermCellProliferativeTypeWriter)
