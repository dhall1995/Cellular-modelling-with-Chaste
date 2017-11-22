/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "NissenNoiseForce.hpp"
#include "NodeBasedCellPopulation.hpp"

//Create constructor
template<unsigned DIM>
NissenNoiseForce<DIM>::NissenNoiseForce()
    : AbstractForce<DIM>(),
      mNoiseStandardDev(1.0e-3) // default to Value in Nissen paper
{
}

//Create Destructor (trivial)
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

//Return the maximum noise standard deviation
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
