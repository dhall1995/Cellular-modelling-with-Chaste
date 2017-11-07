/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef TESTNODEBASEDMORULA_HPP_
#define TESTNODEBASEDMORULA_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NissenForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"


class TestNodeBasedMorula : public AbstractCellBasedTestSuite
{
public:
    void TestNodeBasedEarlyMorula() throw (Exception)
    {
    	//Node Based simulations don't work in parallel
    	EXIT_IF_PARALLEL;

    	//Instantiate a cell mesh generator of dimension 1x1
    	HoneycombMeshGenerator generator(1, 1);
    	MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

    	//Create a nodes only mesh, construct the mesh with a interaction distance
    	//of 5 (As specified in Nissen)
    	NodesOnlyMesh<2> mesh;
    	mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 2.5);

    	//Instantiate vector to store cell pointers, create pointer to the cellproliferativetype
    	std::vector<CellPtr> cells;
    	MAKE_PTR(TransitCellProliferativeType, p_transit_type);

    	//Create a cell generator with cells obeying the UniformCellCycleModel in 2 dimensions and
    	//generate the cells obeying p_transit_type
        UniformCellCycleModel EarlyMorulaModel;
        EarlyMorulaModel.SetMinCellCycleDuration(18.0);
        EarlyMorulaModel.SetMaxCellCycleDuration(20.0);
    	CellsGenerator<EarlyMorulaModel, 2> cells_generator;
    	cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

    	//Link the cells with the mesh created at the start
    	NodeBasedCellPopulation<2> cell_population(mesh, cells);

    	//Instantiate the simulator cell_population, saving results in NodeBasedMorula,
    	//simulating for 30 seconds
    	OffLatticeSimulation<2> simulator(cell_population);
    	simulator.SetOutputDirectory("NodeBasedMorula");
    	simulator.SetSamplingTimestepMultiple(12);
    	simulator.SetEndTime(60.0);

    	//Make pointer to the NissenForce and add it to the simulator
    	MAKE_PTR(NissenForce<2>, p_force);
    	simulator.AddForce(p_force);

    	//Solve the simulation
    	simulator.Solve();

    }
};

#endif // TESTNODEBASEDMORULA_HPP_
