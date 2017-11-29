
#include "NissenForce.hpp"
#include "PolarityCellProperty.hpp"
#include "TrophectodermCellProliferativeType.hpp"
#include "EpiblastCellProliferativeType.hpp"
#include "PrECellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForce<ELEMENT_DIM,SPACE_DIM>::NissenForce()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mS_ICM_ICM(0.6), // ICM-ICM interaction strength - NOTE: Before TE specification all cells are considered ICM-like in their adhesion properties
     mS_TE_ICM(0.6),  // TE-ICM interaction strength
     mS_TE_EPI(0.6),  // TE-EPI interaction strength
     mS_TE_PrE(0.4),  // TE-PrE interaction strength
     mS_TE_TE(-1.4),  // TE-TE interaction strength - NOTE: This is just a prefactor and polarity effects will be included
     mS_PrE_PrE(0.4), // PrE-PrE interaction strength
     mS_PrE_EPI(0.4), // Pre-EPI interaction strength
     mS_PrE_ICM(0.4), // PrE-ICM interaxction strength
     mS_EPI_EPI(0.6), // EPI-EPI interaction strength
     mS_EPI_ICM(0.6), // EPI-ICM interaction strength
     mGrowthDuration(3.0)
{
}

// NOTE: TROPHECTODERM CUTOFF IS 2.5 CELL RADII (Essentially the cutoff for polarity-polarity interactions)

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForce<ELEMENT_DIM,SPACE_DIM>::~NissenForce()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> NissenForce<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                            unsigned nodeBGlobalIndex,
                                                                                            AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);
    
    // Assign labels to each node in the pair
    Node<SPACE_DIM>* p_node_A = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_B = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Find locations of each node in the pair
    const c_vector<double, SPACE_DIM>& r_node_A_location = p_node_A->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_B_location = p_node_B->rGetLocation();

    // Work out the vector from node A to node B and use the GetVector method from rGetMesh
    c_vector<double, SPACE_DIM> unit_vector_from_A_to_B;
    unit_vector_from_A_to_B = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_A_location, r_node_B_location);

    // Distance between the two nodes
    double d = norm_2(unit_vector_from_A_to_B);
    
    // Normalise the vector between A and B
    unit_vector_from_A_to_B /= d;
    
    // NISSEN DISTANCES ARE GIVEN IN UNITS OF CELL RADII
    d = 2.0*d;
    
    // Get ages of cells
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();
    
    // Check that the cells actually have ages
    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));
    
    // Work out the actual force acting between A and B (up to a constant defining the adhesion for a cell-cell pair)
    c_vector<double, SPACE_DIM> potential_gradient;
    potential_gradient = -exp(-d/5.0)*unit_vector_from_A_to_B*5.0;
    potential_gradient_repulsion = - exp(-d)*unit_vector_from_A_to_B;
    c_vector<double, SPACE_DIM> force;
    c_vector<double, SPACE_DIM> zeroes;
    
    /*
     * FIRST WE DEAL WITH TROPHECTODERM INTERACTIONS
     *  - WE TEST WHETHER BOTH CELLS ARE TROPHECTODERM AND THEN IF THEY ARE POLAR. IF THEY AREN'T THEN WE ASSUME
     *    THAT THE INTERACTIONS TAKE THE FORM OF TE-ICM INTERACTIONS FOR THE TIME BEING. AT A LATER TIME I WILL NEED TO IMPLEMENT
     *    TE-PRE INTERACTIONS WHICH WILL INVOLVE USING A DIFFERENT ADHESION FACTOR
     *  - IF BOTH CELLS ARE TROPHECTODERM THEN WE TEST WHETHER BOTH CELLS ARE POLAR. IF BY SOME CHANCE THEY AREN'T BOTH POLAR THEN
     *    INTERACTIONS ARE ASSUMED TO BE MODELLED ON ICM-ICM INTERACTIONS
     *  - IF BOTH CELLS ARE POLAR THEN WE PROCEED TO WORK OUT THE POLARITY VECTORS AND THE ASSOCIATED ADHESION FACTOR
     *  - SANITY CHECKS ARE INCLUDED TO MAKE SURE THAT THE POLARITY VECTORS ARE NON-ZERO AND UNIT VECTORS
     */
    
    // CASE 1: Cell A is trophectoderm
    if(p_cell_A->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
    {
       //CASE 1-1: Cell B is also trophectoderm
       if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
          // For POLAR throphectoderm cells we restrict the distance of interaction to 2.5 cell radii (half of for normal cells)
            // No cells should ever interact beyond the cutoff length
            if (this->mUseCutOffLength)
            {
                if (d >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            /*
            * NOTE WE WANT POLAR TROPHECTODERM INTERACTIONS TO ONLY HAPPEN BETWEEN NEAREST NEIGHBOUR CELLS. 
            * IN ORDER TO AVOID 'BUNCHING' WE SAY THAT TWO TROPHECTODERM CELLS HAVE A POLAR INTERACTION IF THERE IS NO
            * CELL 'BETWEEN' THEM WITHIN THE INTERACTION DISTANCE I.E. FOR CELLS A AND B THERE DOES NOT EXIST A CELL C 
            * SUCH THAT (DISTANCE_FROM_A_TO_C)
            */
            
            // Fill vectors using the polarity_vector data which should be stored when specifiying trophectoderm (See TestNodeBasedMorula.hpp)
            double angle_A = p_cell_A->GetCellData()->GetItem("Polarity Angle");
            double angle_B = p_cell_B->GetCellData()->GetItem("Polarity Angle");
            
         
            double cell_difference_angle = atan2(unit_vector_from_A_to_B[1],unit_vector_from_A_to_B[0]);
          
            double polarity_factor = -0.5*cos(angle_A - angle_B) + 0.5*(angle_A + angle_B - 2.0*cell_difference_angle);
          
            /*
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
               
                 //If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 // with the age of the cells.
                 
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_TE - 5.0)*growth_factor;
                force = -potential_gradient*polarity_factor*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*polarity_factor*s_TE_TE + potential_gradient_repulsion;
                return force;
            }
            */
            force = force = -potential_gradient*polarity_factor*mS_TE_TE + potential_gradient_repulsion;
                return force;
       }
       //CASE 1-2: Cell B is epiblast
       else if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_EPI - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_TE_EPI + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 1-3: Cell B is Undetermined ICM
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
       {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_TE_ICM + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 1-4: Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
       {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_PrE - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_TE_PrE + potential_gradient_repulsion;
                return force;
            }
       }
       else
       {
          return zeroes;
       }
    }
   
   /*
    * Now we deal with interactions of other cells. For the time being cells have been considered in ordered pairs for clarity.
    * Individual attraction factors should be set at the top of the document. Listing cell pairs allows easy and precise 
    * manipulation of attractions between cell lineages. 
    *
    * NOTE: whilst all cells carry a CellPolaritySrn Model it should have zero effect unless we specify in this force law
    */
    //CASE 2: Cell A is Undertermined ICM
    else if(p_cell_A->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
    {
       //CASE 2-1: Cell B is Undertermined ICM
       if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
        {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_ICM_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_ICM_ICM + potential_gradient_repulsion;
                return force;
            }
        }
       //CASE 2-2: Cell B is Epiblast
       else if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_EPI_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_EPI_ICM + potential_gradient_repulsion;
                return force;
            }
        }
       //CASE 2-3: Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
        {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_PrE_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_PrE_ICM + potential_gradient_repulsion;
                return force;
            }
        }
       //CASE 2-4: Cell B is Trophectoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
        {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_TE_ICM + potential_gradient_repulsion;
                return force;
            }
        }
       else
       {
          return zeroes;
       }
    }
    //CASE 3 Cell A is Epiblast
    else if(p_cell_A->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
    {
       //CASE 3-1 Cell B is Epiblast
       if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_EPI_EPI - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_EPI_EPI + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 3-2 Cell B is Undetermined ICM
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_EPI_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_EPI_ICM + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 3-3 Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_PrE_EPI - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_PrE_EPI + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 3-4 Cell B is Trophectoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_EPI - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_TE_EPI + potential_gradient_repulsion;
                return force;
            }
       }
       else
       {
          return zeroes;
       }
    }
    //CASE 4 Cell A is Primitive Endoderm
    else if(p_cell_A->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
    {
       //CASE 4-1 Cell B is Undetermined ICM
       if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_PrE_ICM - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_PrE_ICM + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 4-2 Cell B is Epiblast
       else if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_PrE_EPI - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_PrE_EPI + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 4-3 Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_PrE_PrE - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_PrE_PrE + potential_gradient_repulsion;
                return force;
            }
       }
       //CASE 4-4 Cell B is Trophectoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
          // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion length between the cells grows linearly
                 * with the age of the cells.
                 */
                double growth_factor = std::min(ageA, ageB)/(mGrowthDuration);
                double s = 5.0 + (mS_TE_PrE - 5.0)*growth_factor;
                force = -potential_gradient*s + potential_gradient_repulsion;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                force = -potential_gradient*mS_TE_PrE + potential_gradient_repulsion;
                return force;
            }
       }
       else
       {
          return zeroes;
       }
    }
    else
    {
       return zeroes;
    }
            
          
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_ICM_ICM()
{
    return mS_ICM_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_ICM_ICM(double s)
{
    mS_ICM_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_TE_ICM()
{
    return mS_TE_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_TE_ICM(double s)
{
    mS_TE_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_TE_EPI()
{
    return mS_TE_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_TE_EPI(double s)
{
    mS_TE_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_TE_PrE()
{
    return mS_TE_PrE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_TE_PrE(double s)
{
    mS_TE_PrE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_TE_TE()
{
    return mS_TE_TE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_TE_TE(double s)
{
    mS_TE_TE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_PrE()
{
    return mS_PrE_PrE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_PrE(double s)
{
    mS_PrE_PrE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_EPI()
{
    return mS_PrE_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_EPI(double s)
{
    mS_PrE_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_ICM()
{
    return mS_PrE_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_ICM(double s)
{
    mS_PrE_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_EPI_EPI()
{
    return mS_EPI_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_EPI_EPI(double s)
{
    mS_EPI_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS_EPI_ICM()
{
    return mS_EPI_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS_EPI_ICM(double s)
{
    mS_EPI_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetGrowthDuration(double GrowthDuration)
{
    mGrowthDuration = GrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<S_TE_TE>" << mS_TE_TE << "</S_TE_TE>\n";
    *rParamsFile << "\t\t\t<S_TE_ICM>" << mS_TE_ICM << "</S_TE_ICM>\n";
    *rParamsFile << "\t\t\t<S_TE_EPI>" << mS_TE_EPI << "</S_TE_EPI>\n";
    *rParamsFile << "\t\t\t<S_TE_PrE>" << mS_TE_PrE << "</S_TE_PrE>\n";
    *rParamsFile << "\t\t\t<S_ICM_ICM>" << mS_ICM_ICM << "</S_ICM_ICM>\n";
    *rParamsFile << "\t\t\t<S_PrE_ICM>" << mS_PrE_ICM << "</S_PrE_ICM>\n";
    *rParamsFile << "\t\t\t<S_EPI_ICM>" << mS_EPI_ICM << "</EPI_ICM>\n";
    *rParamsFile << "\t\t\t<S_PrE_EPI>" << mS_PrE_EPI << "</S_PrE_EPI>\n";
    *rParamsFile << "\t\t\t<S_PrE_PrE>" << mS_PrE_PrE << "</S_PrE_PrE>\n";
    *rParamsFile << "\t\t\t<S_ICM_ICM>" << mS_ICM_ICM << "</S_ICM>\n";
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

//Explicit Instantiation of the Force
template class NissenForce<1,1>;
template class NissenForce<1,2>;
template class NissenForce<2,2>;
template class NissenForce<1,3>;
template class NissenForce<2,3>;
template class NissenForce<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NissenForce)
