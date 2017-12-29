
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
    //TRACE("Now attempting to update cell data within CellPolarityTrackingModifier");
    // Make sure the cell population is updated
    rCellPopulation.Update();
    
    //for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         //cell_iter != rCellPopulation.End();
         //++cell_iter)
    //{
	//CellPolaritySrnModel* p_srn_model = static_cast<CellPolaritySrnModel*>(cell_iter->GetSrnModel());

	// NOTE: Here we assert that the cell does actually have the right SRN model
	//assert(p_srn_model != nullptr);

        //this_alpha = p_srn_model->GetPolarityAngle();

        // Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        //cell_iter->GetCellData()->SetItem("Polarity Angle", this_alpha);
    //}

    //Iterate over the population to compute and store each cell's neighbouring polarity angle in the srn model
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        //For our cell get the srn model
	CellPolaritySrnModel* p_srn_model_A = static_cast<CellPolaritySrnModel*>(cell_iter->GetSrnModel());
	
	// NOTE: Here we assert that the cell does actually have the right SRN model
	assert(p_srn_model_A != nullptr);
	//Use the srn model to get the polarity angle for 
        double this_alpha = p_srn_model_A->GetPolarityAngle();
	cell_iter->GetCellData()->SetItem("Polarity Angle", this_alpha);
	    
	//Test if the cell is trophectoderm
	if (cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>() == true)
        {
	    //If the cell is trophectoderm then initialise the sum of neighbouring polarities to update the current polarity
	    double sum_sin_angles = 0.0;
	    //Get cell index for the cell in question	
	    unsigned cell_A_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
		
	    //Iterate through the cells again to first check if a taget cell is close and then if it is trophectoderm  
	    for (typename AbstractCellPopulation<DIM>::Iterator cell_B_iter = rCellPopulation.Begin();
		 cell_B_iter != rCellPopulation.End(); ++cell_B_iter)
	    {
		    unsigned cell_B_index = rCellPopulation.GetLocationIndexUsingCell(*cell_B_iter);

    		    // Assign labels to each node in the pair
    		    Node<DIM>* p_node_A = rCellPopulation.GetNode(cell_A_index);
   		    Node<DIM>* p_node_B = rCellPopulation.GetNode(cell_B_index);
		    
		    // Find locations of each node in the pair
    		    const c_vector<double, DIM>& r_node_A_location = p_node_A->rGetLocation();
    		    const c_vector<double, DIM>& r_node_B_location = p_node_B->rGetLocation();

    		    // Work out the vector from node A to node B and use the GetVector method from rGetMesh
    		    c_vector<double, DIM> vector_from_A_to_B;
                    vector_from_A_to_B = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_A_location, r_node_B_location);

    	            // Distance between the two nodes
    		    double d = norm_2(vector_from_A_to_B);
		    if (d < 2.5)
		    {
			    if(cell_B_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>() == true)
			    {
				  CellPolaritySrnModel* p_srn_model_B = static_cast<CellPolaritySrnModel*>(cell_B_iter->GetSrnModel());
				  // NOTE: Here we assert that the cell does actually have the right SRN model
	    			  assert(p_srn_model_B != nullptr);
	    			  //Use the srn model to get the polarity angle for cell B
            			  double alpha_B_cell = p_srn_model_B->GetPolarityAngle();
				  cell_B_iter->GetCellData()->SetItem("Polarity Angle", alpha_B_cell);
				  //Work out the correct polarity force
                        	  sum_sin_angles += sin(this_alpha - alpha_B_cell); 
			    }
		    }
	    }
	    //whatever our final sum_sin_angles is must be applied to the current trophectoderm cell
	    cell_iter->GetCellData()->SetItem("dVpdAlpha", sum_sin_angles);
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
