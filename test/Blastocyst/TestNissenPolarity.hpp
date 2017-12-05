#ifndef TESTNISSENPOLARITY_HPP_
#define TESTNISSENPOLARITY_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"

// Cell cycle models
#include "PreCompactionCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"

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
#include "Debug.hpp"

class TestNissenPolarity : public AbstractCellBasedWithTimingsTestSuite
{
 private: 
  void GenerateTrophectodermCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<TrophectodermCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            std::vector<double> initial_conditions;
            initial_conditions.push_back(0.0);

            NoCellCycleModel* p_cc_model = new NoCellCycleModel();
            p_cc_model->SetDimension(2);

            CellPolaritySrnModel* p_srn_model = new CellPolaritySrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_prolif_type);

            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);

            p_cell->GetCellData()->SetItem("target area", 1.0);
            rCells.push_back(p_cell);
        }
    }
  
  public:
    void TestInLine()
    {
        // Node-based simulations don't work in parallel
    	  EXIT_IF_PARALLEL;
        
        // Create a simple mesh
        unsigned num_ghosts = 0;
        HoneycombMeshGenerator generator(1, 10, num_ghosts);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 2.5);
        
        // Set up cells, one for each Node
        std::vector<CellPtr> cells;
        GenerateTrophectodermCells(mesh.GetNumNodes(),cells);

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        
        // Set population to output all data to results files
        cell_population.AddCellWriter<CellProliferativeTypesCountWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedNissenPolarityTest");
        simulator.SetOutputDivisionLocations(true);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(30.0);
        
        MAKE_PTR(NissenForce<2>, p_force);
        simulator.AddForce(p_force);
        p_force->SetCutOffLength(2.5);
        
        // Add Cell Polarity modifier
        MAKE_PTR(CellPolarityTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
        
        simulator.Solve();     
    }
}

#endif //TESTNISSENPOLARITY_HPP_
