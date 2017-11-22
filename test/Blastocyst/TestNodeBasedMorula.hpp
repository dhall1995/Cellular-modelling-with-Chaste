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
    
    /*
    *Function to call when we wish to make trophectoderm specification at E3.5. This is done via assigining
    *cells the polarity cell property, a cell label and the trophectoderm cell proliferative type.
    *
    *Trophectoderm Cells are also given initial polarity angles which are initialised to point away from the
    *centroid of the cell population.
    */
    void LabelTrophectodermCells(NodeBasedCellPopulation<2>& cell_population)
    {
        const c_vector<double, 2> morula_centre = cell_population.GetCentroidOfCellPopulation();
        
        //Make the relevant pointers to cell properties for polarity, labels
        MAKE_PTR(PolarityCellProperty, p_pol);
        MAKE_PTR(CellLabel, p_label);
        
        //Make pointer to the trophectoderm cell proliferative type
        MAKE_PTR(TrophectodermCellProliferativeType, p_troph);
	
	//re-set all polarity angles to zero
	cell_population.SetDataOnAllCells("Polarity Angle", 0.0);

	for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {  
            std::vector<double> initial_conditions;

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
                initial_conditions.push_back = angle;
                p_srn_model->SetInitialConditions(initial_conditions);

            }
        }
    }

	


public:
    void TestNodeBasedEarlyMorula(double x) throw (Exception)
    {
	//Node Based simulations don't work in parallel
    	EXIT_IF_PARALLEL;

        //Re-seed Random Number Generator to run multiple ctest runs witout recompiling 
        RandomNumberGenerator::Instance()->Reseed(x);

    	/*
	* Set of methods to generate our initial cell with WildTypeMutationState and TransitCellProliferativeType, the 
	* cell polarity srn model, birthed at time 0.0, in 2 dimensions, with a target area of 1.0
	*/
	
	//shared pointers to the mutation state and proliferative type
	boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
	
	//initial condition for the polarity angle
	std::vector<double> initial_polarity_angle;
	initial_polarity_angle.push_back(0.0);
        
	//Initialise pointers to the cell cycle model and the srn model for the cell
	PreCompactionCellCycleModel* p_cc_model = new PreCompactionCellCycleModel();
	p_cc_model->SetDimension(2);
	
	CellPolaritySrnModel* p_srn_model = new CellPolaritySrnModel();
	p_srn_model->SetInitialConditions(initial_polarity_angle);
	
	//Create a pointer to the cell with the given srn, cell cycle and mutation state, set it's proliferative type to transit
	CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
        p_cell->SetCellProliferativeType(p_prolif_type);
	
	//Set the cell's birth time to 0.0
	double birth_time = 0.0;
        p_cell->SetBirthTime(birth_time);
	
	//Make the target area of the cells 1.0 (non-dimensionalised), the polarity angle 0.0 and add to the vector of cells
	p_cell->GetCellData()->SetItem("target area", 1.0);
	p_cell->GetCellData()->SetItem("Polarity Angle", 0.0);
	
	std::vector<CellPtr> rCells; 
        rCells.push_back(p_cell);
	
	/*
	* Methods to generate the mesh on which our initialised cell will be associated with. We want to make a Nodes only
	* Mesh. 
	*/
	//Instantiate a cell mesh generator of dimension 1x1
    	HoneycombMeshGenerator generator(1, 1);
    	MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

    	//Create a nodes only mesh, construct the mesh with a interaction distance
    	//of 5 (As specified in Nissen)
    	NodesOnlyMesh<2> mesh;
    	mesh.ConstructNodesWithoutMesh(*p_generating_mesh,2.5);

        
	/*
	 * Joining our cell population to our mesh, Instantiating a simulator and setting the time-steps, sampling steps,
	 * simulation time, etc. 
	 */
	
    	//Link the cells with the mesh created at the start
    	NodeBasedCellPopulation<2> cell_population(mesh, rCells);

    	//Instantiate the simulator cell_population, saving results in NodeBasedMorula,
    	//simulating for with SIMULATOR_END_TIME hours
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
	
	//Make sure to add the simulation modifier for tracking cell polarity
        MAKE_PTR(CellPolarityTrackingModifier<2>, p_pol_tracking_modifier);
        simulator.AddSimulationModifier(p_pol_tracking_modifier);

    	//Solve the simulation the first time round
    	simulator.Solve();
	
	/*
	* At this point we should be at an early morula stage of development and ready to specify our outer cells as
	* trophectoderm which we do via the private member function LabelTrophectodermCells to which we pass our\
	* cell_population.
	*/
        
        //Make trophectoderm specification and add a writer for cell proliferative types
        LabelTrophectodermCells(cell_population);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        
        //Run Simulator for a small amount more time in order to allow trophectoderm cells to reach equilibirum
        simulator.SetEndTime(SIMULATOR_END_TIME + 10.0);
        simulator.Solve();
        

    }
};

#endif // TESTNODEBASEDMORULA_HPP_
