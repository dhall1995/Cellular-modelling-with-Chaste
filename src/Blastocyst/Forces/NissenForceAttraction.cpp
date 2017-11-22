
#include "NissenForceAttraction.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::NissenForceAttraction()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mS(0.7),
     mBeta(2.0),
     mGrowthDuration(3.0)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::~NissenForceAttraction()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                            unsigned nodeBGlobalIndex,
                                                                                            AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);
    
    //Assign labels to each node in the pair
    Node<SPACE_DIM>* p_node_A = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<SPACE_DIM>* p_node_B = rCellPopulation.GetNode(nodeBGlobalIndex);

    //Find locations of each node in the pair
    const c_vector<double, SPACE_DIM>& r_node_A_location = p_node_A->rGetLocation();
    const c_vector<double, SPACE_DIM>& r_node_B_location = p_node_B->rGetLocation();

    //Work out the vector from node A to node B and use the GetVector method from rGetMesh
    c_vector<double, SPACE_DIM> unit_vector_from_A_to_B;
    unit_vector_from_A_to_B = rCellPopulation.rGetMesh().GetVectorFromAtoB(r_node_A_location, r_node_B_location);

    //Distance between the two nodes
    double d = norm_2(unit_vector_from_A_to_B);
    
    //Normalise the vector between A and B
    unit_vector_from_A_to_B /= d;
    
    //NISSEN DISTANCES ARE GIVEN IN UNITS OF RADII
    d = 2.0*d;
    
    //Get Ages of Cells
    CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

    double ageA = p_cell_A->GetAge();
    double ageB = p_cell_B->GetAge();
    
    //Check that the cells actually have ages
    assert(!std::isnan(ageA));
    assert(!std::isnan(ageB));
    
    
    //Work out the actual force acting between A and B
    c_vector<double, SPACE_DIM> potential_gradient;
    potential_gradient = mS*exp(-d/mBeta)*unit_vector_from_A_to_B/mBeta;
    c_vector<double, SPACE_DIM> force;

    if (this->mUseCutOffLength)
    {
        if (d/2.0 >= this->GetCutOffLength())  //remember chaste distances given in DIAMETERS
        {
            return force;
        }
    }
    
    if(ageA < mGrowthDuration && ageB < mGrowthDuration)
    {
        force = potential_gradient;
        /*
         * If the cells are both newly divided, then the repulsion between the cells grows linearly
         * with the age of the cells.
         */
        force = (std::min(ageA, ageB)/mGrowthDuration)*force;
        return force;
    }
    else //if no other conditions are met then return the force from Nissen
    {
        force = potential_gradient;
        return force;
    }
}

//Defining the GetS method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::GetS()
{
    return mS;
}

//defining the SetS method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::SetS(double s)
{
    mS = s;
}

//defining the GetBeta Method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::GetBeta()
{
    return mBeta;
}

//Defining the SetBeta Method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::SetBeta(double beta)
{
    mBeta = beta;
}

//defining the GetGrowthDuration Method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

//Defining the SetGrowthDuration Method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::SetGrowthDuration(double GrowthDuration)
{
    mBeta = GrowthDuration;
}

//Saving the force parameters to rParamsFile
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForceAttraction<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Beta>" << mBeta << "</Beta>\n";
    *rParamsFile << "\t\t\t<S>" << mS << "</S>\n";
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";
    AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(rParamsFile);
}

//Explicit Instantiation of the Force
template class NissenForceAttraction<1,1>;
template class NissenForceAttraction<1,2>;
template class NissenForceAttraction<2,2>;
template class NissenForceAttraction<1,3>;
template class NissenForceAttraction<2,3>;
template class NissenForceAttraction<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(NissenForceAttraction)
