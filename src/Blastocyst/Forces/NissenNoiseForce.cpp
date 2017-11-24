
#include "NissenNoiseForce.hpp"
#include "NodeBasedCellPopulation.hpp"

// Constructor
template<unsigned DIM>
NissenNoiseForce<DIM>::NissenNoiseForce()
    : AbstractForce<DIM>(),
      mNoiseStandardDev(1.0e-3) // default to Value in Nissen paper
{
}

// Destructor (trivial)
template<unsigned DIM>
NissenNoiseForce<DIM>::~NissenNoiseForce()
{
}

//Tool to set the Noise Standard Deviation in the system
template<unsigned DIM>
void NissenNoiseForce<DIM>::SetNoiseStandardDev(double newValue)
{
    assert(newValue > 0.0);
    mNoiseStandardDev = newValue;
}

// Return the maximum noise standard deviation
template<unsigned DIM>
double NissenNoiseForce<DIM>::GetNoiseStandardDev()
{
    return mNoiseStandardDev;
}

template<unsigned DIM>
void NissenNoiseForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        // Define a force contribution for each node taken from a normal distribution
        c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            double xi = RandomNumberGenerator::Instance()->NormalRandomDeviate(0.0, mNoiseStandardDev);

            force_contribution[i] = xi;
        }
        node_iter->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void NissenNoiseForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class NissenNoiseForce<1>;
template class NissenNoiseForce<2>;
template class NissenNoiseForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NissenNoiseForce)
