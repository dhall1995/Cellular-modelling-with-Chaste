
#include "CellPolarityTrackingModifier.hpp"
#include "CellPolaritySrnModel.hpp"
#include "TrophectodermCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CellPolarityTrackingModifier<DIM>::CellPolarityTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellPolarityTrackingModifier<DIM>::~CellPolarityTrackingModifier()
{
}

template<unsigned DIM>
void CellPolarityTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
//    TRACE("Now attempting UpdateAtEndOfTimeStep within the CellPolarityTrackingModifier");
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellPolarityTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
//    TRACE("Now attempting SetupSolve within the CellPolarityTrackingModifier");
     /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellPolarityTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    NodeBasedCellPopulation<DIM>* p_population = static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);
    //TRACE("Now attempting to update cell data within CellPolarityTrackingModifier");
    // Make sure the cell population is updated
    p_population.Update();

    // First recover each cell's polarity angle from the ODEs and store in CellData. Keep the polarity angle as a variable outside the scope of the for loop
    double this_alpha = 0.0;
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_population.Begin();
         cell_iter != p_population.End();
         ++cell_iter)
    {
//		TRACE("Now attempting to get a cell's polarity angle within UpdateCellData within CellPolarityTrackingModifier");
//		bool variable = cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>();
//		TRACE("Are we working with a trophectoderm cell??");
//		PRINT_VARIABLE(variable);
		CellPolaritySrnModel* p_srn_model = static_cast<CellPolaritySrnModel*>(cell_iter->GetSrnModel());

		// NOTE: Here we assert that the cell does actually have the right SRN model
		assert(p_srn_model != nullptr);

        this_alpha = p_srn_model->GetPolarityAngle();

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        cell_iter->GetCellData()->SetItem("Polarity Angle", this_alpha);
    }

    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_population.Begin();
         cell_iter != p_population.End();
         ++cell_iter)
    {
        if (cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>() == true)
        {
            // Get the set of neighbouring location indices
            //std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
            std::set<unsigned> neighbour_indices = p_population.GetNodesWithinNeighbourhoodRadius(p_population.GetLocationIndexUsingCell(*cell_iter),1.25);
            
	    // Compute this trophectoderm cell's neighbouring trophectoderm cells and store in
            // CellData the sin of the angle differences
            if (!neighbour_indices.empty())
            {
//                TRACE("Now working out the polarity potential for a trophectoderm cell")
		        double sum_sin_angles = 0.0;
                
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                     iter != neighbour_indices.end();
                     ++iter)
                {
                    CellPtr p_cell = p_population.GetCellUsingLocationIndex(*iter);
                    
                    if (cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>() == true)
                    {
                        double alpha_p_cell = p_cell->GetCellData()->GetItem("Polarity Angle");
                        sum_sin_angles += sin(this_alpha - alpha_p_cell);
                    }
                }
                cell_iter->GetCellData()->SetItem("dVpdAlpha", sum_sin_angles);
            }
            else
            {
                // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
                cell_iter->GetCellData()->SetItem("dVpdAlpha", 0.0);
            }
        }
        else
        {
            // For non-trophectoderm cells we just set the polarity potential to zero - we don't care what happens to their polarity angle so we just let it evolve via random noise
            cell_iter->GetCellData()->SetItem("dVpdAlpha",0.0);
        }
    }
}

template<unsigned DIM>
void CellPolarityTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CellPolarityTrackingModifier<1>;
template class CellPolarityTrackingModifier<2>;
template class CellPolarityTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellPolarityTrackingModifier)
