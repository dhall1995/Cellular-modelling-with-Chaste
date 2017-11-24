
#include "NissenForceAttractionTest.hpp"
#include "PolarityCellProperty.hpp"
#include "TrophectodermCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::NissenForceAttractionTest()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mS_ICM_ICM(0.7), // ICM-ICM interaction strength - NOTE: Before TE specification all cells are considered ICM-like in their adhesion properties
     mS_TE_ICM(0.7),  // TE-ICM interaction strength
     mS_TE_EPI(0.7),  // TE-EPI interaction strength
     mS_TE_PrE(0.4),  // TE-PrE interaction strength
     mS_TE_TE(-1.5),  // TE-TE interaction strength - NOTE: This is just a prefactor and polarity effects will be included
     mS_PrE_PrE(0.4), // PrE-PrE interaction strength
     mS_PrE_EPI(0.4), // Pre-EPI interaction strength
     mS_PrE_ICM(0.4), // PrE-ICM interaxction strength
     mS_EPI_EPI(0.7), // EPI-EPI interaction strength
     mS_EPI_ICM(0.7), // EPI-ICM interaction strength
     mBeta(2.0),
     mGrowthDuration(3.0)
{
}

// NOTE: TROPHECTODERM CUTOFF IS 2.5 CELL RADII (Essentially the cutoff for polarity-polarity interactions)

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::~NissenForceAttractionTest()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
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
    potential_gradient = -exp(-d/mBeta)*unit_vector_from_A_to_B/mBeta;
    c_vector<double, SPACE_DIM> force;
    
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
    
    // Test Whether nodes are both of the trophectoderm cell proliferative type
    if(p_cell_A->GetCellProliferativeType() == p_cell_B->GetCellProliferativeType())
    {
        // Test whether both nodes have polarity
        if (p_cell_A->HasCellProperty<PolarityCellProperty>() && p_cell_B->HasCellProperty<PolarityCellProperty>())
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
            
            // Inititalise vectors to work out polarity interaction
            c_vector<double, SPACE_DIM> polarity_factor_node_A;
            c_vector<double, SPACE_DIM> polarity_factor_node_B;
            
            // Fill vectors using the polarity_vector data which should be stored when specifiying trophectoderm (See TestNodeBasedMorula.hpp)
            double angle_A = p_cell_A->GetCellData()->GetItem("Polarity Angle");
            double angle_B = p_cell_B->GetCellData()->GetItem("Polarity Angle");

//            double polarity_factor_node_A[0] = p_node_A->GetCellData()->GetItem("polarity_vector_x_value");
//            double polarity_factor_node_A[1] = p_node_A->GetCellData()->GetItem("polarity_vector_y_value");
//
//            double polarity_factor_node_B[0] = p_node_B->GetCellData()->GetItem("polarity_vector_x_value");
//            double polarity_factor_node_B[1] = p_node_B->GetCellData()-GetItem("polarity_vector_y_value");
//
//            assert((polarity_factor_node_A[0]*polarity_factor_node_A[0] + polarity_factor_node_A[1]*polarity_factor_node_A[1]) == 1.0);
//            assert((polarity_factor_node_B[0]*polarity_factor_node_B[0] + polarity_factor_node_B[1]*polarity_factor_node_B[1]) == 1.0);
//
            double cell_difference_angle = atan(unit_vector_from_A_to_B[1]/unit_vector_from_A_to_B[0]);

            
            // Need to work out (polarity_factor_node_A X r_AB) . (polarity_factor_node_B X r_BA)
            // We only need the z-component so we start with (polarity_factor_node_A X r_AB)[2]
//            double node_A_factor = polarity_factor_node_A[0]*unit_vector_from_A_to_B[1] - polarity_vector_node_A[1]*unit_vector_from_A_to_B[0];
            
            // Now we store the z-component of (polarity_factor_node_B X r_BA)
//            double node_B_factor = polarity_factor_node_B[1]*unit_vector_from_A_to_B[0] - polarity_vector_node_B[0]*unit_vector_from_A_to_B[1];
            
            // Finally we store the polarity factor for the pair of TE cells
//            double polarity_factor = node_A_factor*node_B_factor;
            
            double polarity_factor = -0.5*cos(angle_A - angle_B) + 0.5*(angle_A + angle_B - 2.0*cell_difference_angle);
            
            force = -potential_gradient*mS_TE_TE*polarity_factor;
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion between the cells grows linearly
                 * with the age of the cells.
                 */
                force = (std::min(ageA, ageB)/mGrowthDuration)*force;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                return force;
            }
        }
        else
        {
            // No cells should ever interact beyond the cutoff length OF 5.0 Cell Radii
            if (this->mUseCutOffLength)
            {
                if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
                {
                    return force;
                }
            }
            
            // If one or more of the TE Cells isn't polar for any reason then the attraction reverts back to that of undifferentiated cells
            force = -potential_gradient*mS_ICM_ICM;
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
                /*
                 * If the cells are both newly divided, then the repulsion between the cells grows linearly
                 * with the age of the cells.
                 */
                force = (std::min(ageA, ageB)/mGrowthDuration)*force;
                return force;
            }
            else // if no other conditions are met then return the force
            {
                return force;
            }
        }
    }

    // WE NOW DEAL WITH INTERACTIONS NOT INVOLVING TROPHECTODERM CELLS. FOR THE TIME BEING THIS ONLY INVOLVES ICM-ICM INTERACTIONS

    // No cells should ever interact beyond the cutoff length
    if (this->mUseCutOffLength)
    {
        if (d/2.0 >= this->GetCutOffLength())  // remember chaste distances given in DIAMETERS
        {
            return force;
        }
    }
    
    if (ageA < mGrowthDuration && ageB < mGrowthDuration)
    {
        force -= potential_gradient*mS_ICM_ICM;
        /*
         * If the cells are both newly divided, then the repulsion between the cells grows linearly
         * with the age of the cells.
         */
        force = (std::min(ageA, ageB)/mGrowthDuration)*force;
        return force;
    }
    else // if no other conditions are met then return the force from Nissen
    {
        force -= potential_gradient*mS_ICM_ICM;
        return force;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_ICM_ICM()
{
    return mS_ICM_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_ICM_ICM(double s)
{
    mS_ICM_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_TE_ICM()
{
    return mS_TE_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_TE_ICM(double s)
{
    mS_TE_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_TE_EPI()
{
    return mS_TE_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_TE_EPI(double s)
{
    mS_TE_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_TE_PrE()
{
    return mS_TE_PrE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_TE_PrE(double s)
{
    mS_TE_PrE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_TE_TE()
{
    return mS_TE_TE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_TE_TE(double s)
{
    mS_TE_TE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_PrE()
{
    return mS_PrE_PrE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_PrE(double s)
{
    mS_PrE_PrE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_EPI()
{
    return mS_PrE_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_EPI(double s)
{
    mS_PrE_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_ICM()
{
    return mS_PrE_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_ICM(double s)
{
    mS_PrE_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_EPI_EPI()
{
    return mS_EPI_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_EPI_EPI(double s)
{
    mS_EPI_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetS_EPI_ICM()
{
    return mS_EPI_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetS_EPI_ICM(double s)
{
    mS_EPI_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetBeta()
{
    return mBeta;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetBeta(double beta)
{
    mBeta = beta;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::SetGrowthDuration(double GrowthDuration)
{
    mBeta = GrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttractionTest<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Beta>" << mBeta << "</Beta>\n";
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

//Explicit Instantiation of the Force
template class NissenForceAttractionTest<1,1>;
template class NissenForceAttractionTest<1,2>;
template class NissenForceAttractionTest<2,2>;
template class NissenForceAttractionTest<1,3>;
template class NissenForceAttractionTest<2,3>;
template class NissenForceAttractionTest<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NissenForceAttractionTest)
