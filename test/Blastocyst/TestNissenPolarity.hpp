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

class TestNodeBasedPolarity : public AbstractCellBasedWithTimingsTestSuite
{
  public: 
    void TestInLine()
    {
      // Node-based simulations don't work in parallel
    	EXIT_IF_PARALLEL;
      
      HoneycombMeshGenerator generator(1, 10);
      MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
      
      NodesOnlyMesh<2> mesh;
      mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 2.5);
      
      std::vector<CellPtr> cells;
      MAKE_PTR(TrophectodermCellProliferativeType, p_troph_type);
      CellsGenerator<NoCellCycleModel, 2> cells_generator;
      cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_troph_type);
      
      NodeBasedCellPopulation<2> cell_population(mesh, cells);
      
      OffLatticeSimulation<2> simulator(cell_population);
      simulator.SetOutputDirectory("NodeBasedNissenPolarityTest");
      simulator.SetSamplingTimestepMultiple(12);
      simulator.SetEndTime(30.0);
      
      MAKE_PTR(NissenForce<2>, p_force);
      simulator.AddForce(p_force);
      
      simulator.Solve();
      
    }
}

