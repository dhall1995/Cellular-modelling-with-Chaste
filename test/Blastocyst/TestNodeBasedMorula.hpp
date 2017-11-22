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
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"

//Cycle Models
#include "PreCompactionCellCycleModel.hpp"

//Cell Properties
#include "PolarityCellProperty.hpp"
#include "CellLabel.hpp"

//Proliferative Types
#include "TrophectodermCellProliferativeType.hpp"
#include "EpiblastCellProliferativeType.hpp"
#include "PrECellProliferativeType.hpp"


//Mesh Generators
#include "HoneycombMeshGenerator.hpp"

//Force Models
#include "NissenForceAttractionTest.hpp"
#include "NissenForceRepulsion.hpp"
#include "NissenNoiseForce.hpp"

//Simulation Files
#include "OffLatticeSimulation.hpp"
#include "CellPolaritySrnModel.hpp"
#include "CellPolarityTrackingModifier.hpp"

#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

//Visualising
#include "CellProliferativeTypesCountWriter.hpp"
#include "Debug.hpp"


class TestNodeBasedMorula : public AbstractCellBasedTestSuite
{
private:
    double SIMULATOR_END_TIME = 80.0;
    
    //Function to call when we wish to make trophectoderm specification at E3.5. Trophectoderm Cells are also given initial polarity vectors  
    void LabelTrophectodermCells(NodeBasedCellPopulation<2>& cell_population)
    {
        const c_vector<double, 2> morula_centre = cell_population.GetCentroidOfCellPopulation();
        
        //Make the relevant pointers to cell properties for polarity, labels
        MAKE_PTR(PolarityCellProperty, p_pol);
        MAKE_PTR(CellLabel, p_label);
        
        //Make pointer to the trophectoderm cell proliferative type
        MAKE_PTR(TrophectodermCellProliferativeType, p_troph);

	    //Make pointer to the Cell Polarity SRN Model
	    CellPolaritySrnModel* p_pol_srn = new CellPolaritySrnModel();

	    for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            //Set initial polarity vectors on every cell to zero
            cell_iter->GetCellData()->SetItem("Polarity Angle", 0.0);
            
            std::vector<double> initial_conditions;
            initial_conditions.push_back(0.0);

            if (cell_population.GetNeighbouringNodeIndices(cell_population.GetLocationIndexUsingCell(*cell_iter)).size() <= 5.0)
            {
                
                //Add cell properties for labels, and polarity
                cell_iter->AddCellProperty(p_pol);
                cell_iter->AddCellProperty(p_label);
                cell_iter->SetCellProliferativeType(p_troph);

                //Get Cell Location to work out polarity vector
                const c_vector<double, 2> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
                c_vector<double,2> unit_vector_from_centroid_to_cell = morula_centre - cell_location;
			
                //norm of the distance from the centroid to the cell
                double d = sqrt(unit_vector_from_centroid_to_cell[0]*unit_vector_from_centroid_to_cell[0] + unit_vector_from_centroid_to_cell[1]*unit_vector_from_centroid_to_cell[1]);

                //Normalise our vector from centroid to cell
                unit_vector_from_centroid_to_cell /= d;

                //Store the vector x and y values in cell data for that cell
                double cell_x_value = unit_vector_from_centroid_to_cell[0];
                double cell_y_value = unit_vector_from_centroid_to_cell[1];

                double angle = atan(cell_y_value/cell_x_value);
                cell_iter->GetCellData()->SetItem("Polarity Angle", angle);

                //Initialise a std::vector to store the calculated initial angle as the initial condition for
                //the srn model
                initial_conditions[0] = angle;
                p_pol_srn->SetInitialConditions(initial_conditions);

            }
            
            //Give the cell the polarity srn model
            cell_iter->SetSrnModel(p_pol_srn);
        }
    }

	


public:
    void TestNodeBasedEarlyMorula() throw (Exception)
    {
    	//Give the cell the polarity srn model
        cell_iter->SetSrnModel(p_pol_srn);
	
	//Node Based simulations don't work in parallel
    	EXIT_IF_PARALLEL;

        //Re-seed Random Number Generator
        RandomNumberGenerator::Instance()->Reseed(57);

    	//Instantiate a cell mesh generator of dimension 1x1
    	HoneycombMeshGenerator generator(1, 1);
    	MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

    	//Create a nodes only mesh, construct the mesh with a interaction distance
    	//of 5 (As specified in Nissen)
    	NodesOnlyMesh<2> mesh;
    	mesh.ConstructNodesWithoutMesh(*p_generating_mesh,2.5);

    	//Instantiate vector to store cell pointers, create pointer to the cellproliferativetype
    	std::vector<CellPtr> cells;
    	MAKE_PTR(TransitCellProliferativeType, p_transit_type);

    	//Create a cell generator with cells obeying the PreCompactionCycleModel/UniformCellCycleModel in 2 dimensions and
    	//generate the cells obeying p_transit_type
    	CellsGenerator<PreCompactionCellCycleModel, 2> cells_generator;
    	cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

    	//Link the cells with the mesh created at the start
    	NodeBasedCellPopulation<2> cell_population(mesh, cells);
	

        for(unsigned index=0; index < cell_population.rGetMesh().GetNumNodes(); index++)
     	{
		cell_population.rGetMesh().GetNode(index)->SetRadius(0.5);
        }

    	//Instantiate the simulator cell_population, saving results in NodeBasedMorula,
    	//simulating for 30 seconds
    	OffLatticeSimulation<2> simulator(cell_population);
    	simulator.SetOutputDirectory("NodeBasedMorula");
    	simulator.SetSamplingTimestepMultiple(24);
        double Dt = simulator.GetDt();
        Dt /= 2.0;
        simulator.SetDt(Dt);
    	simulator.SetEndTime(SIMULATOR_END_TIME);

    	//Make pointer to the NissenForce and add it to the simulator
    	
        MAKE_PTR(NissenForceRepulsion<2>, p_force_repulsion);
        MAKE_PTR(NissenForceAttractionTest<2>, p_force_attraction);

        p_force_repulsion->SetCutOffLength(2.5);
        p_force_attraction->SetCutOffLength(2.5);

    	simulator.AddForce(p_force_attraction);
        simulator.AddForce(p_force_repulsion);
        

        //Make pointer to the NissenNoiseForce and add it to the simulator
        MAKE_PTR(NissenNoiseForce<2>, p_noise_force);
        simulator.AddForce(p_noise_force);

    	//Solve the simulation
    	simulator.Solve();
        
        //Make trophectoderm specification and add a writer for cell proliferative types
        LabelTrophectodermCells(cell_population);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        //Make sure to add the simulation modifier for tracking cell polarity
        MAKE_PTR(CellPolarityTrackingModifier<2>, p_pol_tracking_modifier);
        simulator.AddSimulationModifier(p_pol_tracking_modifier);
        
        //Run Simulator for a small amount more time in order to allow trophectoderm cells to reach equilibirum
        simulator.SetEndTime(SIMULATOR_END_TIME + 10.0);
        simulator.Solve();
        

    }
};

#endif // TESTNODEBASEDMORULA_HPP_
