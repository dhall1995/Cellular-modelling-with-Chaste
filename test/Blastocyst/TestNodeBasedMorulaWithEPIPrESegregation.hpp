#ifndef TESTNODEBASEDMORULAWITHEPIPRESEGREGATION_HPP_
#define TESTNODEBASEDMORULAWITHEPIPRESEGREGATION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"

// Cell cycle models
#include "PreCompactionCellCycleModel.hpp"

// Cell properties
#include "PolarityCellProperty.hpp"
#include "CellLabel.hpp"

// Cell proliferative types
#include "TrophectodermCellProliferativeType.hpp"
#include "EpiblastCellProliferativeType.hpp"
#include "PrECellProliferativeType.hpp"

// Mesh generators
#include "HoneycombMeshGenerator.hpp"

// Force models
#include "NissenForce.hpp"
#include "NissenNoiseForce.hpp"

// Division Rules
#include "NissenBasedDivisionRule.hpp"

// Simulation files
#include "OffLatticeSimulation.hpp"
#include "CellPolaritySrnModel.hpp"
#include "CellPolarityTrackingModifier.hpp"

#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

// Visualising
#include "CellProliferativeTypesCountWriter.hpp"
#include "PolarityVectorWriter.hpp"
#include "Debug.hpp"

class TestNodeBasedMorulaWithEPIPrESegregation : public AbstractCellBasedWithTimingsTestSuite
{
private:
    double SIMULATOR_END_TIME = 220.0;
    
    /*
     * Function to call when we wish to make trophectoderm specification at E3.5. This is done by assigning
     * cells the polarity cell property, a cell label and the trophectoderm cell proliferative type.
     *
     * Trophectoderm cells are also given polarity angles, which are initialised to point away from the
     * centroid of the cell population.
     */
    void LabelEpiblastPrECells(NodeBasedCellPopulation<2>& cell_population)
    {

        /*
         * Make the relevant pointers to cell properties for trophectoderm polarity, label and proliferative type
         *
         * NOTE: When we call this method, the cell population has already called
         *       CellPropertyRegistry::TakeOwnership(). We must therefore access
         *       any new or existing cell properties via the existing cell property
         *       registry, as implemented below.
         */
        CellPropertyRegistry* p_registry = cell_population.Begin()->rGetCellPropertyCollection().GetCellPropertyRegistry();
        boost::shared_ptr<AbstractCellProperty> p_epi_label = p_registry->Get<CellLabel>();      
        boost::shared_ptr<AbstractCellProperty> p_pre_label = p_registry->Get<CellLabel>();
        boost::shared_ptr<AbstractCellProperty> p_epi = p_registry->Get<EpiblastCellProliferativeType>();
        boost::shared_ptr<AbstractCellProperty> p_pre = p_registry->Get<PrECellProliferativeType>();

        for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
	    std::set<unsigned> neighbour_indices = cell_population.GetNeighbouringNodeIndices(node_index);
            double random = RandomNumberGenerator::Instance()->ranf();
            if (random < 0.5)
            {
                // Add cell properties for labels, and polarity
                cell_iter->AddCellProperty(p_epi_label);
                cell_iter->SetCellProliferativeType(p_epi);
		cell_iter->GetCellData()->SetItem("target area", 2.0);
            }
            else
            {
                cell_iter->AddCellProperty(p_pre_label);
                cell_iter->SetCellProliferativeType(p_pre);
		cell_iter->GetCellData()->SetItem("target area", 1.0);
            }
        }
    }

public:
    void TestNodeBasedEPIPrESegregation() throw (Exception)
    {
        // Node-based simulations don't work in parallel
    	EXIT_IF_PARALLEL;

        // Re-seed random number generator to run multiple ctest runs without recompiling
        RandomNumberGenerator::Instance()->Reseed(50);

    	/*
	     * Set of methods to generate our initial cell with WildTypeMutationState and TransitCellProliferativeType, the
	     * cell polarity SRN model, birthed at time 0.0, in 2 dimensions, with a target area of 1.0
	     */
	
        // Shared pointers to the mutation state and proliferative type
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        // Initialise pointer to the cell cycle model for the cell
        PreCompactionCellCycleModel* p_cc_model = new PreCompactionCellCycleModel();
        p_cc_model->SetDimension(2);

        // Create a pointer to the cell with the given cell cycle and mutation state, set it's proliferative type to transit
        CellPtr p_cell(new Cell(p_state, p_cc_model));
        p_cell->SetCellProliferativeType(p_prolif_type);
	
        // Set the cell's birth time to 0.0
        double birth_time = 0.0;
        p_cell->SetBirthTime(birth_time);
	
        // Make the target area of the cells 1.0 (non-dimensionalised), the polarity angle 0.0 and add to the vector of cells
        p_cell->GetCellData()->SetItem("target area", 1.0);

        std::vector<CellPtr> rCells;
        rCells.push_back(p_cell);
	
        /*
         * Methods to generate the mesh on which our initialised cell will be associated with. We want to make a Nodes only
         * Mesh.
         */
	// Instantiate a cell mesh generator of dimension 1x1
    	HoneycombMeshGenerator generator(1, 1);
    	MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

    	// Create a nodes only mesh, construct the mesh with a interaction distance of 5 (as specified in Nissen)
    	NodesOnlyMesh<2> mesh;
    	mesh.ConstructNodesWithoutMesh(*p_generating_mesh,6.0);

        /*
         * Joining our cell population to our mesh, Instantiating a simulation and setting the time-steps, sampling steps,
         * simulation time, etc.
         */
	
    	// Link the cells with the mesh created at the start
    	NodeBasedCellPopulation<2> cell_population(mesh, rCells);

	

    	// Instantiate the simulation, saving results in NodeBasedMorula, simulating for SIMULATOR_END_TIME hours
    	OffLatticeSimulation<2> simulation(cell_population);
    	simulation.SetOutputDirectory("NodeBasedMorulaWithEPIPrESegregation");
    	simulation.SetSamplingTimestepMultiple(24);
        double dt = 0.5*simulation.GetDt();
        simulation.SetDt(dt);
    	simulation.SetEndTime(SIMULATOR_END_TIME);

    	// Make pointer to the NissenForce and add it to the simulation
    	MAKE_PTR(NissenForce<2>, p_force);
	p_force->SetCutOffLength(2.5);

        simulation.AddForce(p_force);

        // Make pointer to the NissenNoiseForce and add it to the simulation
        MAKE_PTR(NissenNoiseForce<2>, p_noise_force);
        simulation.AddForce(p_noise_force);

    	// Solve the simulation the first time round
    	simulation.Solve();
        TRACE("finished first simulation up to early morula");
        
	
      /*
       * At this point we should be at an early morula stage of development and ready to specify our outer cells as
       * trophectoderm which we do via the private member function LabelTrophectodermCells to which we pass our\
       * cell_population.
       */

      // Make trophectoderm specification and add a writer for cell proliferative types
      LabelEpiblastPrECells(cell_population);
      cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        
      // Run simulation for a small amount more time in order to allow trophectoderm cells to reach equilibirum
      simulation.SetEndTime(SIMULATOR_END_TIME + 40.0);
      simulation.Solve();
      
      double number_of_PrE_cells = 0.0;
      double number_of_isolated_PrE_cells = 0.0;

      for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
      {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);
	    std::set<unsigned> neighbour_indices = cell_population.GetNeighbouringNodeIndices(node_index);
            if (cell_iter->GetCellProliferativeType()->template IsType<PrECellProliferativeType>() == true)
            {
                 number_of_PrE_cells += 1.0;
		 double number_of_neighbouring_PrE_cells = 0.0;
		 for (std::set<unsigned>::iterator cell_B_iter = neighbour_indices.begin();
		 cell_B_iter != neighbour_indices.end(); ++cell_B_iter)
		 {
		 	if(cell_iter->GetCellProliferativeType()->template IsType<PrECellProliferativeType>() == true)
			{
				number_of_neighbouring_PrE_cells += 1.0;
			}
		 }
		 if(number_of_neighbouring_PrE_cells < 2.0)
		 {
		 	number_of_isolated_PrE_cells += 1.0;
		}
            }

      }
      double efficiency = 1.0 - (number_of_isolated_PrE_cells/number_of_PrE_cells);
      PRINT_VARIABLE(efficiency);
      PRINT_VARIABLE(number_of_PrE_cells);
      PRINT_VARIABLE(number_of_isolated_PrE_cells);
    }
};

#endif //TESTNODEBASEDMORULAWITHEPIPRESEGREGATION_HPP_
