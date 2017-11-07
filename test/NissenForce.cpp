
#include "NissenForce.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NissenForce<ELEMENT_DIM,SPACE_DIM>::NissenForce()
   : AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM>(),
     mS(0.6),
     mBeta(5.0)
{
}

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

    //Work out the actual force acting between A and B
    c_vector<double, SPACE_DIM> potential_gradient;
    potential_gradient = exp(-d)*unit_vector_from_A_to_B - mS*exp(-d/mBeta)*unit_vector_from_A_to_B/mBeta;
    c_vector<double, SPACE_DIM> force;
    force -= potential_gradient;
    return force;
}

//Defining the GetS method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetS()
{
    return mS;
}

//defining the SetS method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetS(double s)
{
    mS = s;
}

//defining the GetBeta Method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NissenForce<ELEMENT_DIM,SPACE_DIM>::GetBeta()
{
    return mBeta;
}

//Defining the SetBeta Method
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::SetBeta(double beta)
{
    mBeta = beta;
}

//Saving the force parameters to rParamsFile
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NissenForce<ELEMENT_DIM,SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Beta>" << mBeta << "</Beta>\n";
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
