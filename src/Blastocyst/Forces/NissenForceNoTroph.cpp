
#include "NissenForceNoTroph.hpp"
#include "PolarityCellProperty.hpp"
#include "TrophectodermCellProliferativeType.hpp"
#include "EpiblastCellProliferativeType.hpp"
#include "PrECellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellPolaritySrnModel.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::NissenForceNoTroph()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mS_ICM_ICM(0.6), // ICM-ICM interaction strength - NOTE: Before TE specification all cells are considered ICM-like in their adhesion properties
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
NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::~NissenForceNoTroph()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
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
    c_vector<double, SPACE_DIM> force;
    c_vector<double, SPACE_DIM> zeroes;
    
   /*
    * We deal with interactions of other cells. For the time being cells have been considered in ordered pairs for clarity.
    * Individual attraction factors should be set at the top of the document. Listing cell pairs allows easy and precise 
    * manipulation of attractions between cell lineages. 
    *
    * NOTE: whilst all cells carry a CellPolaritySrn Model it should have zero effect unless we specify in this force law
    */
    //CASE 2: Cell A is Undertermined ICM
    if(p_cell_A->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
    {
       //CASE 2-1: Cell B is Undertermined ICM
       if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
        {
            
            double s = mS_ICM_ICM;
            
            //if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            //{
               //AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

               //std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

               //if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
               //{
                  // Spring rest length increases from a small value to the normal rest length over 1 hour
                  //s = 5.0 + (mS_ICM_ICM - 5.0) * ageA/mGrowthDuration;
               //}
               //if (ageA + SimulationTime::Instance()->GetTimeStep() >= mGrowthDuration)
               //{
                  // This spring is about to go out of scope
                  //p_static_cast_cell_population->UnmarkSpring(cell_pair);
               //}
            //}
          
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
        }
       //CASE 2-2: Cell B is Epiblast
       else if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {

            
            double s = mS_EPI_ICM;
          
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
        }
       //CASE 2-3: Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
        {        
            double s = mS_PrE_ICM;

            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
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
            double s = mS_EPI_EPI;
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
               AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

               std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

               if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
               {
                  // Spring rest length increases from a small value to the normal rest length over 1 hour
                  s = 5.0 + (mS_EPI_EPI - 5.0) * ageA/mGrowthDuration;
               }
               if (ageA + SimulationTime::Instance()->GetTimeStep() >= mGrowthDuration)
               {
                  // This spring is about to go out of scope
                  p_static_cast_cell_population->UnmarkSpring(cell_pair);
               }
            }
          
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
       }
       //CASE 3-2 Cell B is Undetermined ICM
       else if(p_cell_B->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
       {        
            double s = mS_EPI_ICM;
    
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
       }
       //CASE 3-3 Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
       {
            double s = mS_PrE_EPI;
         
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
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
            double s = mS_PrE_ICM;
          
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
       }
       //CASE 4-2 Cell B is Epiblast
       else if(p_cell_B->GetCellProliferativeType()->template IsType<EpiblastCellProliferativeType>())
       {
            double s = mS_PrE_EPI;
            
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
       }
       //CASE 4-3 Cell B is Primitive Endoderm
       else if(p_cell_B->GetCellProliferativeType()->template IsType<PrECellProliferativeType>())
       {
            double s = mS_PrE_PrE;
            
            if (ageA < mGrowthDuration && ageB < mGrowthDuration)
            {
               AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

               std::pair<CellPtr,CellPtr> cell_pair = p_static_cast_cell_population->CreateCellPair(p_cell_A, p_cell_B);

               if (p_static_cast_cell_population->IsMarkedSpring(cell_pair))
               {
                  // Spring rest length increases from a small value to the normal rest length over 1 hour
                  s = 5.0 + (mS_PrE_PrE - 5.0) * ageA/mGrowthDuration;
               }
               if (ageA + SimulationTime::Instance()->GetTimeStep() >= mGrowthDuration)
               {
                  // This spring is about to go out of scope
                  p_static_cast_cell_population->UnmarkSpring(cell_pair);
               }
            }
          
            force = potential_gradient*s + potential_gradient_repulsion;
            return force;
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
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetS_ICM_ICM()
{
    return mS_ICM_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetS_ICM_ICM(double s)
{
    mS_ICM_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_PrE()
{
    return mS_PrE_PrE;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_PrE(double s)
{
    mS_PrE_PrE = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_EPI()
{
    return mS_PrE_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_EPI(double s)
{
    mS_PrE_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetS_PrE_ICM()
{
    return mS_PrE_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetS_PrE_ICM(double s)
{
    mS_PrE_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetS_EPI_EPI()
{
    return mS_EPI_EPI;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetS_EPI_EPI(double s)
{
    mS_EPI_EPI = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetS_EPI_ICM()
{
    return mS_EPI_ICM;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetS_EPI_ICM(double s)
{
    mS_EPI_ICM = s;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::SetGrowthDuration(double GrowthDuration)
{
    mGrowthDuration = GrowthDuration;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceNoTroph<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
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
template class NissenForceNoTroph<1,1>;
template class NissenForceNoTroph<1,2>;
template class NissenForceNoTroph<2,2>;
template class NissenForceNoTroph<1,3>;
template class NissenForceNoTroph<2,3>;
template class NissenForceNoTroph<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NissenForceNoTroph)
