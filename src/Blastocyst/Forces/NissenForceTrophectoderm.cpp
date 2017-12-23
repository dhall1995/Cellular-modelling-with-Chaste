
#include "NissenForceTrophectoderm.hpp"
#include "PolarityCellProperty.hpp"
#include "TrophectodermCellProliferativeType.hpp"
#include "EpiblastCellProliferativeType.hpp"
#include "PrECellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellPolaritySrnModel.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::NissenForceTrophectoderm()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mS_TE_ICM(0.6),  // TE-ICM interaction strength
     mS_TE_EPI(0.6),  // TE-EPI interaction strength
     mS_TE_PrE(0.4),  // TE-PrE interaction strength
     mS_TE_TE(-1.4),  // TE-TE interaction strength - NOTE: This is just a prefactor and polarity effects will be included
     mGrowthDuration(3.0)
{
}

// NOTE: TROPHECTODERM CUTOFF IS 2.5 CELL RADII (Essentially the cutoff for polarity-polarity interactions)

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::~NissenForceTrophectoderm()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
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
    c_vector<double, SPACE_DIM> potential_gradient_repulsion;
    potential_gradient = exp(-d/5.0)*unit_vector_from_A_to_B/5.0;
    potential_gradient_repulsion = -exp(-d)*unit_vector_from_A_to_B;
    c_vector<double, SPACE_DIM> zeroes;
    zeroes[0] = 0.0;
    zeroes[1] = 0.0;
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
    
    // CASE 1: Cell A is trophectoderm
    if(p_cell_A->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
    {
       //Initialise vectors to hold each focus of trophectoderm A
       c_vector<double, SPACE_DIM> p_cell_A_first_focus;
       c_vector<double, SPACE_DIM> p_cell_A_second_focus;
       
       // Fill vectors using the polarity_vector data which should be stored when specifiying trophectoderm (See TestNodeBasedMorula.hpp)
       CellPolaritySrnModel* p_srn_model_A = static_cast<CellPolaritySrnModel*>(p_cell_A->GetSrnModel());
       double angle_A = p_srn_model_A->GetPolarityAngle();
       //initialise the perpendicular vector to the polarity of the trophectoderm cell A to define the axis of our pseudo-elipse
       c_vector<double, SPACE_DIM> perp_polarity_vector_A;
       perp_polarity_vector_A[0] = sin(angle_A);
       perp_polarity_vector_A[1] = -cos(angle_A);
          
       p_cell_A_first_focus = r_node_A_location + 0.25*perp_polarity_vector_A;
       p_cell_A_second_focus = r_node_A_location -0.25*perp_polarity_vector_A;
       
       //CASE 1-1: Cell B is also trophectoderm
       if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
            TRACE("calculating TE-TE interaction");
           
            //First thing we want to do is get the polarity angle for the trophectoderm cell B
            CellPolaritySrnModel* p_srn_model_B = static_cast<CellPolaritySrnModel*>(p_cell_B->GetSrnModel());
            double angle_B = p_srn_model_B->GetPolarityAngle();
          
            double s = mS_TE_TE;
            
            //if each of the cells is young then we treat them as spherical with radius 2.0 whilst they grow
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
               double cell_difference_angle = atan2(unit_vector_from_A_to_B[1],unit_vector_from_A_to_B[0]);
               double polarity_factor = -sin(cell_difference_angle - angle_A)*sin(cell_difference_angle - angle_B);
               
               potential_gradient = exp(-d/10.0)*unit_vector_from_A_to_B/5.0;
               potential_gradient_repulsion = -exp(-d/2.0)*unit_vector_from_A_to_B;
               
               AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

               std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

               if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
               {
                   //Spring rest length increases from a small value to the normal rest length over 1 hour
                   if(polarity_factor < 0.0)
                   {
                     s = -5.0 + (mS_TE_TE + 5.0) * ageA/mGrowthDuration;
                   }
               }
               if (ageA + SimulationTime::Instance()->GetTimeStep() >= mGrowthDuration)
               {
                  // This spring is about to go out of scope
                  p_static_cast_cell_population->UnmarkSpring(cell_pair);
               }
               force = potential_gradient*polarity_factor*s + potential_gradient_repulsion;
               return force;
            }

            //otherwise we retrieve the focii of cell B
            c_vector<double, SPACE_DIM> p_cell_B_first_focus;
            c_vector<double, SPACE_DIM> p_cell_B_second_focus;
          
            //initialise the perpendicular vector to the polarity of the trophectoderm cell B
            c_vector<double, SPACE_DIM> perp_polarity_vector_B;
            perp_polarity_vector_B[0] = sin(angle_B);
            perp_polarity_vector_B[1] = -cos(angle_B);
            
            //Define the two focii for cellB
            p_cell_B_first_focus = r_node_B_location + 0.25*perp_polarity_vector_B;
            p_cell_B_second_focus = r_node_B_location -0.25*perp_polarity_vector_B;
       
            //Initialise the distances between the focii
            double d_A1_B1;
            double d_A2_B1;
            double d_A1_B2;
            double d_A2_B2;
          
            //Initialise the forces between each focus
            c_vector<double, SPACE_DIM> force_first_A_focus_first_B_focus;
            c_vector<double, SPACE_DIM> force_second_A_focus_second_B_focus;
            c_vector<double, SPACE_DIM> force_first_A_focus_second_B_focus;
            c_vector<double, SPACE_DIM> force_second_A_focus_first_B_focus;
          
            //Set up vectors between the various focii
            c_vector<double, SPACE_DIM> unit_vector_from_A1_to_B1 = -p_cell_A_first_focus + p_cell_B_first_focus;
            c_vector<double, SPACE_DIM> unit_vector_from_A2_to_B1 = -p_cell_A_second_focus + p_cell_B_first_focus;
            c_vector<double, SPACE_DIM> unit_vector_from_A2_to_B2 = -p_cell_A_second_focus + p_cell_B_second_focus;
            c_vector<double, SPACE_DIM> unit_vector_from_A1_to_B2 = -p_cell_A_first_focus + p_cell_B_second_focus;
          
            //set the distances between the various focii
            d_A1_B1 = norm_2(unit_vector_from_A1_to_B1);
            d_A2_B1 = norm_2(unit_vector_from_A2_to_B1);
            d_A1_B2 = norm_2(unit_vector_from_A1_to_B2);
            d_A2_B2 = norm_2(unit_vector_from_A2_to_B2);
          
            PRINT_VARIABLE(d_A1_B2);
            PRINT_VARIABLE(d_A2_B1);
            PRINT_VARIABLE(d_A1_B1);
            PRINT_VARIABLE(d_A2_B2);
          
            //normalise our vectors between the focii
            unit_vector_from_A1_to_B1 /= d_A1_B1;
            unit_vector_from_A2_to_B1 /= d_A2_B1;
            unit_vector_from_A1_to_B2 /= d_A1_B2;
            unit_vector_from_A2_to_B2 /= d_A2_B2;
          
            //Nissen distances given in radii
            d_A1_B1 *= 2.0;
            d_A1_B2 *= 2.0;
            d_A2_B1 *= 2.0;
            d_A2_B2 *= 2.0;
          
          
            //keep track of how many interactions are non-zero (we want the force normalised) as if it was the action of a single cell
            double number_of_active_forces = 0.0;
          
            //Now consider the forces generated by the action of each focus on each other
            if(d_A1_B1 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A1_B1/5.0)*unit_vector_from_A1_to_B1/5.0;
               potential_gradient_repulsion = -exp(-d_A1_B1)*unit_vector_from_A1_to_B1;
               
               double cell_difference_angle_A1_B1 = atan2(unit_vector_from_A1_to_B1[1],unit_vector_from_A1_to_B1[0]);
               double polarity_factor_A1_B1 = -sin(cell_difference_angle_A1_B1 - angle_A)*sin(cell_difference_angle_A1_B1 - angle_B);
               
               force_first_A_focus_first_B_focus = potential_gradient*polarity_factor_A1_B1*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A1_B2 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A1_B2/5.0)*unit_vector_from_A1_to_B2/5.0;
               potential_gradient_repulsion = -exp(-d_A1_B2)*unit_vector_from_A1_to_B2;
               
               double cell_difference_angle_A1_B2 = atan2(unit_vector_from_A1_to_B2[1],unit_vector_from_A1_to_B2[0]);
               double polarity_factor_A1_B2 = -sin(cell_difference_angle_A1_B2 - angle_A)*sin(cell_difference_angle_A1_B2 - angle_B);
               
               force_first_A_focus_first_B_focus = potential_gradient*polarity_factor_A1_B2*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A2_B1 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A2_B1/5.0)*unit_vector_from_A2_to_B1/5.0;
               potential_gradient_repulsion = -exp(-d_A2_B1)*unit_vector_from_A2_to_B1;
               
               double cell_difference_angle_A2_B1 = atan2(unit_vector_from_A2_to_B1[1],unit_vector_from_A2_to_B1[0]);
               double polarity_factor_A2_B1 = -sin(cell_difference_angle_A2_B1 - angle_A)*sin(cell_difference_angle_A2_B1 - angle_B);
               
               force_first_A_focus_first_B_focus = potential_gradient*polarity_factor_A2_B1*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A2_B2 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A2_B2/5.0)*unit_vector_from_A2_to_B2/5.0;
               potential_gradient_repulsion = -exp(-d_A2_B2)*unit_vector_from_A2_to_B2;
               
               double cell_difference_angle_A2_B2 = atan2(unit_vector_from_A2_to_B2[1],unit_vector_from_A2_to_B2[0]);
               double polarity_factor_A2_B2 = -sin(cell_difference_angle_A2_B2 - angle_A)*sin(cell_difference_angle_A2_B2 - angle_B);
               
               force_first_A_focus_first_B_focus = potential_gradient*polarity_factor_A2_B2*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_first_A_focus_first_B_focus + force_first_A_focus_second_B_focus + force_second_A_focus_first_B_focus + force_second_A_focus_second_B_focus);
               return force;
            }
          
             
       }
       //CASE 1-2: Cell B is epiblast
       else if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {
            TRACE("calculating TE-EPI interaction");
          
            //Initialise the distances between the focii and the centre of cell B
            double d_A1_B;
            double d_A2_B;
          
            //Initialise the forces between each focus and the centre of cell B
            c_vector<double, SPACE_DIM> force_first_A_focus_B;
            c_vector<double, SPACE_DIM> force_second_A_focus_B;
          
            //Set up vectors between the various focii and the centre of cell B
            c_vector<double, SPACE_DIM> unit_vector_from_A1_to_B = -p_cell_A_first_focus + r_node_B_location;
            c_vector<double, SPACE_DIM> unit_vector_from_A2_to_B = -p_cell_A_second_focus + r_node_B_location;

            //set the distances between the various focii
            d_A1_B = norm_2(unit_vector_from_A1_to_B);
            d_A2_B = norm_2(unit_vector_from_A2_to_B);
          
            PRINT_VARIABLE(d_A1_B);
            PRINT_VARIABLE(d_A2_B);
          
            //normalise our vectors between the focii
            unit_vector_from_A1_to_B /= d_A1_B;
            unit_vector_from_A2_to_B /= d_A2_B;
          
            //nissen distances given in radii
            d_A1_B *= 2.0;
            d_A1_B *= 2.0;
          
            double s = mS_TE_EPI;
            double number_of_active_forces = 0.0;
          
            if(d_A1_B/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A1_B/5.0)*unit_vector_from_A1_to_B/5.0;
               potential_gradient_repulsion = -exp(-d_A1_B)*unit_vector_from_A1_to_B;
               
               force_first_A_focus_B = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A2_B/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A2_B/5.0)*unit_vector_from_A2_to_B/5.0;
               potential_gradient_repulsion = -exp(-d_A2_B)*unit_vector_from_A2_to_B;
               
               force_second_A_focus_B = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
          
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_first_A_focus_B + force_second_A_focus_B)/number_of_active_forces;
               return force;
            }
 
       }
       //CASE 1-3: Cell B is Undetermined ICM
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
       {
            /*
            TRACE("calculating TE-ICM interaction");
          
            //Initialise the distances between the focii and the centre of cell B
            double d_A1_B;
            double d_A2_B;
    
            //Initialise the forces between each focus and the centre of cell B
            c_vector<double, SPACE_DIM> force_first_A_focus_B;
            c_vector<double, SPACE_DIM> force_second_A_focus_B;
          
            //Set up vectors between the various focii and the centre of cell B
            c_vector<double, SPACE_DIM> unit_vector_from_A1_to_B = -p_cell_A_first_focus + r_node_B_location;
            c_vector<double, SPACE_DIM> unit_vector_from_A2_to_B = -p_cell_A_second_focus + r_node_B_location;

            //set the distances between the various focii
            d_A1_B = norm_2(unit_vector_from_A1_to_B);
            d_A2_B = norm_2(unit_vector_from_A2_to_B);
          
            PRINT_VARIABLE(d_A1_B);
            PRINT_VARIABLE(d_A2_B);
          
            //normalise our vectors between the focii
            unit_vector_from_A1_to_B /= d_A1_B;
            unit_vector_from_A2_to_B /= d_A2_B;
          
            //distances given in terms of radii
            d_A1_B *= 2.0;
            d_A1_B *= 2.0;
          
            double s = mS_TE_ICM;
            double number_of_active_forces = 0.0;
          
            if(d_A1_B/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A1_B/5.0)*unit_vector_from_A1_to_B/5.0;
               potential_gradient_repulsion = -exp(-d_A1_B)*unit_vector_from_A1_to_B;
               
               force_first_A_focus_B = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
               PRINT_VARIABLE(number_of_active_forces);
            }
            if(d_A2_B/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A2_B/5.0)*unit_vector_from_A2_to_B/5.0;
               potential_gradient_repulsion = -exp(-d_A2_B)*unit_vector_from_A2_to_B;
               
               force_second_A_focus_B = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
               PRINT_VARIABLE(number_of_active_forces);
            }
          
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_first_A_focus_B + force_second_A_focus_B)/number_of_active_forces;
               return force;
            }
*/
       }
       //CASE 1-4: Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
       {          
            TRACE("calculating TE-PrE interaction");
          
            //Initialise the distances between the focii and the centre of cell B
            double d_A1_B;
            double d_A2_B;
          
            //Initialise the forces between each focus and the centre of cell B
            c_vector<double, SPACE_DIM> force_first_A_focus_B;
            c_vector<double, SPACE_DIM> force_second_A_focus_B;
          
            //Set up vectors between the various focii and the centre of cell B
            c_vector<double, SPACE_DIM> unit_vector_from_A1_to_B = -p_cell_A_first_focus + r_node_B_location;
            c_vector<double, SPACE_DIM> unit_vector_from_A2_to_B = -p_cell_A_second_focus + r_node_B_location;

            //set the distances between the various focii
            d_A1_B = norm_2(unit_vector_from_A1_to_B);
            d_A2_B = norm_2(unit_vector_from_A2_to_B);
          
            PRINT_VARIABLE(d_A1_B);
            PRINT_VARIABLE(d_A2_B);
          
            //normalise our vectors between the focii
            unit_vector_from_A1_to_B /= d_A1_B;
            unit_vector_from_A2_to_B /= d_A2_B;
          
            //nissen distances given in radii
            d_A1_B *= 2.0;
            d_A1_B *= 2.0;
          
            double s = mS_TE_PrE;
            double number_of_active_forces = 0.0;
          
            if(d_A1_B/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A1_B/5.0)*unit_vector_from_A1_to_B/5.0;
               potential_gradient_repulsion = -exp(-d_A1_B)*unit_vector_from_A1_to_B;
               
               force_first_A_focus_B = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A2_B/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A2_B/5.0)*unit_vector_from_A2_to_B/5.0;
               potential_gradient_repulsion = -exp(-d_A2_B)*unit_vector_from_A2_to_B;
               
               force_second_A_focus_B = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
          
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_first_A_focus_B + force_second_A_focus_B)/number_of_active_forces;
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
       /*
       //CASE 2-1: Cell B is Trophectoderm
       if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
            TRACE("calculating ICM-TE interaction");
            //First thing we want to do is get the polarity angle for the trophectoderm cell B
            CellPolaritySrnModel* p_srn_model_B = static_cast<CellPolaritySrnModel*>(p_cell_B->GetSrnModel());
            double angle_B = p_srn_model_B->GetPolarityAngle();
          
            //initialise the perpendicular vector to the polarity of the trophectoderm cell B
            c_vector<double, SPACE_DIM> perp_polarity_vector_B;
            perp_polarity_vector_B[0] = sin(angle_B);
            perp_polarity_vector_B[1] = -cos(angle_B);
            
            //Define the two focii for cellB
            c_vector<double, SPACE_DIM> p_cell_B_first_focus = r_node_B_location + 0.25*perp_polarity_vector_B;
            c_vector<double, SPACE_DIM> p_cell_B_second_focus = r_node_B_location -0.25*perp_polarity_vector_B;
          
            //Initialise the distances between the focii and the centre of cell B
            double d_A_B1;
            double d_A_B2;
          
            //Initialise the forces between each focus and the centre of cell B
            c_vector<double, SPACE_DIM> force_A_first_B_focus;
            c_vector<double, SPACE_DIM> force_A_second_B_focus;
          
            //Set up vectors between the various focii and the centre of cell B
            c_vector<double, SPACE_DIM> unit_vector_from_A_to_B1 = -r_node_A_location + p_cell_B_first_focus;
            c_vector<double, SPACE_DIM> unit_vector_from_A_to_B2 = -r_node_A_location + p_cell_B_second_focus ;

            //set the distances between the various focii
            d_A_B1 = norm_2(unit_vector_from_A_to_B1);
            d_A_B2 = norm_2(unit_vector_from_A_to_B2);
          
            PRINT_VARIABLE(d_A_B1);
            PRINT_VARIABLE(d_A_B2);
          
            //normalise our vectors between the focii
            unit_vector_from_A_to_B1 /= d_A_B1;
            unit_vector_from_A_to_B2 /= d_A_B2;
          
            //Nissen distances given in radii
            d_A_B1 *= 2.0;
            d_A_B2 *= 2.0;
          
            double s = mS_TE_ICM;
            double number_of_active_forces = 0.0;
          
            if(d_A_B1/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A_B1/5.0)*unit_vector_from_A_to_B1/5.0;
               potential_gradient_repulsion = -exp(-d_A_B1)*unit_vector_from_A_to_B1;
               
               force_A_first_B_focus = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A_B2/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A_B2/5.0)*unit_vector_from_A_to_B2/5.0;
               potential_gradient_repulsion = -exp(-d_A_B2)*unit_vector_from_A_to_B2;
               
               force_A_second_B_focus = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
          
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_A_first_B_focus + force_A_second_B_focus)/number_of_active_forces;
               return force;
            }
          
       }
       else
       {
          return zeroes;
       }
       */
    }
    //CASE 3 Cell A is Epiblast
    else if(p_cell_A->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
    {   
       //CASE 3-1 Cell B is Trophectoderm
       if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
            TRACE("calculating EPI-TE interaction");
          
            //First thing we want to do is get the polarity angle for the trophectoderm cell B
            CellPolaritySrnModel* p_srn_model_B = static_cast<CellPolaritySrnModel*>(p_cell_B->GetSrnModel());
            double angle_B = p_srn_model_B->GetPolarityAngle();
          
            //initialise the perpendicular vector to the polarity of the trophectoderm cell B
            c_vector<double, SPACE_DIM> perp_polarity_vector_B;
            perp_polarity_vector_B[0] = sin(angle_B);
            perp_polarity_vector_B[1] = -cos(angle_B);
            
            //Define the two focii for cellB
            c_vector<double, SPACE_DIM> p_cell_B_first_focus = r_node_B_location + 0.25*perp_polarity_vector_B;
            c_vector<double, SPACE_DIM> p_cell_B_second_focus = r_node_B_location -0.25*perp_polarity_vector_B;
          
            //Initialise the distances between the focii and the centre of cell B
            double d_A_B1;
            double d_A_B2;
          
            //Initialise the forces between each focus and the centre of cell B
            c_vector<double, SPACE_DIM> force_A_first_B_focus;
            c_vector<double, SPACE_DIM> force_A_second_B_focus;
          
            //Set up vectors between the various focii and the centre of cell B
            c_vector<double, SPACE_DIM> unit_vector_from_A_to_B1 = -r_node_A_location + p_cell_B_first_focus;
            c_vector<double, SPACE_DIM> unit_vector_from_A_to_B2 = -r_node_A_location + p_cell_B_second_focus ;

            //set the distances between the various focii
            d_A_B1 = norm_2(unit_vector_from_A_to_B1);
            d_A_B2 = norm_2(unit_vector_from_A_to_B2);
          
            PRINT_VARIABLE(d_A_B1);
            PRINT_VARIABLE(d_A_B2);
          
            //normalise our vectors between the focii
            unit_vector_from_A_to_B1 /= d_A_B1;
            unit_vector_from_A_to_B2 /= d_A_B2;
          
            //Nissen distances given in radii
            d_A_B1 *= 2.0;
            d_A_B2 *= 2.0;
          
            double s = mS_TE_EPI;
            double number_of_active_forces = 0.0;
          
            if(d_A_B1/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A_B1/5.0)*unit_vector_from_A_to_B1/5.0;
               potential_gradient_repulsion = -exp(-d_A_B1)*unit_vector_from_A_to_B1;
               
               force_A_first_B_focus = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A_B2/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A_B2/5.0)*unit_vector_from_A_to_B2/5.0;
               potential_gradient_repulsion = -exp(-d_A_B2)*unit_vector_from_A_to_B2;
               
               force_A_second_B_focus = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
          
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_A_first_B_focus + force_A_second_B_focus)/number_of_active_forces;
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
       TRACE("calculating PrE-TE interaction");
       
       //CASE 4-1 Cell B is Trophectoderm
       if(p_cell_B->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
       {
            //First thing we want to do is get the polarity angle for the trophectoderm cell B
            CellPolaritySrnModel* p_srn_model_B = static_cast<CellPolaritySrnModel*>(p_cell_B->GetSrnModel());
            double angle_B = p_srn_model_B->GetPolarityAngle();
          
            //initialise the perpendicular vector to the polarity of the trophectoderm cell B
            c_vector<double, SPACE_DIM> perp_polarity_vector_B;
            perp_polarity_vector_B[0] = sin(angle_B);
            perp_polarity_vector_B[1] = -cos(angle_B);
            
            //Define the two focii for cellB
            c_vector<double, SPACE_DIM> p_cell_B_first_focus = r_node_B_location + 0.25*perp_polarity_vector_B;
            c_vector<double, SPACE_DIM> p_cell_B_second_focus = r_node_B_location -0.25*perp_polarity_vector_B;
          
            //Initialise the distances between the focii and the centre of cell B
            double d_A_B1;
            double d_A_B2;
          
            //Initialise the forces between each focus and the centre of cell B
            c_vector<double, SPACE_DIM> force_A_first_B_focus;
            c_vector<double, SPACE_DIM> force_A_second_B_focus;
          
            //Set up vectors between the various focii and the centre of cell B
            c_vector<double, SPACE_DIM> unit_vector_from_A_to_B1 = -r_node_A_location + p_cell_B_first_focus;
            c_vector<double, SPACE_DIM> unit_vector_from_A_to_B2 = -r_node_A_location + p_cell_B_second_focus ;

            //set the distances between the various focii
            d_A_B1 = norm_2(unit_vector_from_A_to_B1);
            d_A_B2 = norm_2(unit_vector_from_A_to_B2);
          
            PRINT_VARIABLE(d_A_B1);
            PRINT_VARIABLE(d_A_B2);
          
            //normalise our vectors between the focii
            unit_vector_from_A_to_B1 /= d_A_B1;
            unit_vector_from_A_to_B2 /= d_A_B2;
          
            //Nissen distances given in radii
            d_A_B1 *= 2.0;
            d_A_B2 *= 2.0;
          
            double s = mS_TE_PrE;
            double number_of_active_forces = 0.0;
          
            if(d_A_B1/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A_B1/5.0)*unit_vector_from_A_to_B1/5.0;
               potential_gradient_repulsion = -exp(-d_A_B1)*unit_vector_from_A_to_B1;
               
               force_A_first_B_focus = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
            if(d_A_B2/2.0 < this->GetCutOffLength())
            {
               potential_gradient = exp(-d_A_B2/5.0)*unit_vector_from_A_to_B2/5.0;
               potential_gradient_repulsion = -exp(-d_A_B2)*unit_vector_from_A_to_B2;
               
               force_A_second_B_focus = potential_gradient*s + potential_gradient_repulsion;
               number_of_active_forces += 1.0;
            }
          
            if(number_of_active_forces == 0.0)
            {
               return zeroes;
            }
            else
            {
               force = (force_A_first_B_focus + force_A_second_B_focus)/number_of_active_forces;
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
double NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::GetS_TE_ICM()
{
    return mS_TE_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::SetS_TE_ICM(double s)
{
    mS_TE_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::GetS_TE_EPI()
{
    return mS_TE_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::SetS_TE_EPI(double s)
{
    mS_TE_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::GetS_TE_PrE()
{
    return mS_TE_PrE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::SetS_TE_PrE(double s)
{
    mS_TE_PrE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::GetS_TE_TE()
{
    return mS_TE_TE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::SetS_TE_TE(double s)
{
    mS_TE_TE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::SetGrowthDuration(double GrowthDuration)
{
    mGrowthDuration = GrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceTrophectoderm<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<S_TE_TE>" << mS_TE_TE << "</S_TE_TE>\n";
    *rParamsFile << "\t\t\t<S_TE_ICM>" << mS_TE_ICM << "</S_TE_ICM>\n";
    *rParamsFile << "\t\t\t<S_TE_EPI>" << mS_TE_EPI << "</S_TE_EPI>\n";
    *rParamsFile << "\t\t\t<S_TE_PrE>" << mS_TE_PrE << "</S_TE_PrE>\n";
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

//Explicit Instantiation of the Force
template class NissenForceTrophectoderm<1,1>;
template class NissenForceTrophectoderm<1,2>;
template class NissenForceTrophectoderm<2,2>;
template class NissenForceTrophectoderm<1,3>;
template class NissenForceTrophectoderm<2,3>;
template class NissenForceTrophectoderm<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NissenForceTrophectoderm)
