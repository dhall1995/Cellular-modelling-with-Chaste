/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

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
    	TRACE("Now attempting UpdateAtEndOfTimeStep within the CellPolarityTrackingModifier");
	UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellPolarityTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
     TRACE("Now attempting SetupSolve within the CellPolarityTrackingModifier");
     /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CellPolarityTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    TRACE("Now attempting to update cell data within CellPolarityTrackingModifier");
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // First recover each cell's polarity angle from the ODEs and store in CellData. Keep the polarity angle as a variable outside the scope of the for loop
    double this_alpha = 0.0;
    
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

		TRACE("Now attempting to get a cell's polarity angle within UpdateCellData within CellPolarityTrackingModifier");
		bool variable = cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>();
		TRACE("Are we working with a trophectoderm cell??");
		PRINT_VARIABLE(variable);
		CellPolaritySrnModel* p_model = static_cast<CellPolaritySrnModel*>(cell_iter->GetSrnModel());
        	this_alpha = p_model->GetPolarityAngle();

        	// Note that the state variables must be in the same order as listed in DeltaNotchOdeSystem
        	cell_iter->GetCellData()->SetItem("Polarity Angle", this_alpha);

    }

    // Next iterate over the population to compute and store each cell's neighbouring Delta concentration in CellData
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        if(cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>() == true)
        {
            // Get the set of neighbouring location indices
            std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(rCellPopulation.GetLocationIndexUsingCell(*cell_iter));
            
            // Compute this trophectoderm cell's neighbouring trophectoderm cells and store in
            // CellData the sin of the angle differences
            if (!neighbour_indices.empty())
            {
                TRACE("Now working out the polarity potential for a trophectoderm cell")
		double sum_sin_angles = 0.0;
                
                for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                     iter != neighbour_indices.end();
                     ++iter)
                {
                    CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                    
                    if (cell_iter->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>() == true)
                    {
                        double alpha_p_cell = p_cell->GetCellData()->GetItem("Polarity Angle");
                        sum_sin_angles += sin(this_alpha - alpha_p_cell);
                    }
                }
                cell_iter->GetCellData()->SetItem("dVp/dAlpha", sum_sin_angles);
            }
            else
            {
                // If this cell has no neighbours, such as an isolated cell in a CaBasedCellPopulation, store 0.0 for the cell data
                cell_iter->GetCellData()->SetItem("dVp/dAlpha", 0.0);
            }
        }
	else
	{
		// for non-trophectoderm cells we just set the polarity potential to zero - we don't care what happens to their polarity angle so we just let it evolve via random noise
		cell_iter->GetCellData()->SetItem("dVp/dAlpha",0.0);
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
