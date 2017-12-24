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

#include "PolaritySecondFocusVectorWriter.hpp"
#include "TrophectodermCellProliferativeType.hpp"
#include "CellPolaritySrnModel.hpp"
#include "NodeBasedCellPopulation.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PolaritySecondFocusVectorWriter<ELEMENT_DIM, SPACE_DIM>::PolaritySecondFocusVectorWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("PolaritySecondFocusVector.dat")
{
    this->mVtkVectorCellDataName = "Polarity Second Focus Vector";
    this->mOutputScalarData = false;
    this->mOutputVectorData = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PolaritySecondFocusVectorWriter<ELEMENT_DIM, SPACE_DIM>::GetVectorCellDataForVtkOutput(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    assert(this->mOutputVectorData);

    c_vector<double, SPACE_DIM> orientation;
    
    CellPolaritySrnModel* p_srn_model_A = static_cast<CellPolaritySrnModel*>(pCell->GetSrnModel());
    
    // NOTE: Here we assert that the cell does actually have the right SRN model
	  assert(p_srn_model_A != nullptr);
    double angle_A = p_srn_model_A->GetPolarityAngle();
    if(pCell->GetCellProliferativeType()->template IsType<TrophectodermCellProliferativeType>())
    {
        orientation[0] = -sin(angle_A);
        orientation[1] = cos(angle_A);
    }
    else
    {
        orientation[0] = 0.0;
        orientation[1] = 0.0;
    }

    return orientation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolaritySecondFocusVectorWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> cell_orientation = GetVectorCellDataForVtkOutput(pCell, pCellPopulation);
    c_vector<double, SPACE_DIM> cell_location_second_focus = pCellPopulation->GetLocationOfCellCentre(pCell) - 0.25*cell_orientation;		

    *this->mpOutStream << location_index << " " << cell_id << " ";
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location_second_focus[i] << " ";    
    }

    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_orientation[i] << " ";
    }
}

// Explicit instantiation
template class PolaritySecondFocusVectorWriter<1,1>;
template class PolaritySecondFocusVectorWriter<1,2>;
template class PolaritySecondFocusVectorWriter<2,2>;
template class PolaritySecondFocusVectorWriter<1,3>;
template class PolaritySecondFocusVectorWriter<2,3>;
template class PolaritySecondFocusVectorWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PolaritySecondFocusVectorWriter)
