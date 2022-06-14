#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
//#include "../src/Du2013_CalibNeur.hpp"
#include "../src/DummyCell.hpp"
#include "../src/imtiaz_2002d_noTstart_COR.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "Debug.hpp"
#include "AbstractElement.hpp"
#include "Node.hpp"
#include "MonodomainProblem.hpp"
#include "ChasteEllipsoid.hpp"
#include "AbstractElement.hpp"

using namespace std;

struct coordinateV_st
{
  double x;
  double y;
  double z;
  double V;
};

class ICCCellFactory : public AbstractCardiacCellFactory<2,3>
{
private:
  std::vector<coordinateV_st> LaplaceInfo;
public:
  ICCCellFactory() : AbstractCardiacCellFactory<2,3>()
  {
    ReadLaplaceFile();
  }

  void ReadLaplaceFile()
  {
    std::ifstream inLaplaceInfo("projects/mesh/Stomach2D/rat_scaffold_32_32_2_2D_laplace_longi.txt");
    if(!inLaplaceInfo)
    {
      EXCEPTION("Reading laplace solution error");
    }
    std::string line;
    coordinateV_st lapInfo;

    while(std::getline(inLaplaceInfo, line))
    {
      stringstream cordinateLap(line);
      cordinateLap >> lapInfo.x >> lapInfo.y >> lapInfo.z >> lapInfo.V;
      LaplaceInfo.push_back(lapInfo);
    }
  }

  AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
  {
    double x = pNode->rGetLocation()[0];
    double y = pNode->rGetLocation()[1];
    double z = pNode->rGetLocation()[2];

    coordinateV_st info;
    int counter = 0;
    double V_val = 0;
    for(std::vector<coordinateV_st>::iterator itr = LaplaceInfo.begin(); itr!=LaplaceInfo.end();itr++)
    {
      info = *itr;
      if(info.x > x-0.001 && info.x < x+0.001  && info.y > y-0.001 && info.y < y+0.001 && info.z > z-0.001 && info.z < z + 0.001)
      {
        counter++;
        V_val = info.V;
        break;
      }
    }
    if (counter != 1)
    {
      PRINT_4_VARIABLES(x,y,z,counter);
      EXCEPTION("Coordinates not found in Laplace file");
    }

    //CellDu2013_CalibNeurFromCellML* cell = new CellDu2013_CalibNeurFromCellML(mpSolver, mpZeroStimulus);
    //cell->SetParameter("eta", 0.037);

    Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(mpSolver, mpZeroStimulus);
    cell->SetParameter("eta", 0.045);
    if (V_val > 45)
    {
      return new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
    }
    double r = 0.2;
    // for rat_scaffold_32_32_2_2D
    double xv = 0.06;
    double yv = -1.42;
    double zv = -2.54;

    // // for rat_scaffold_64_64_2_2D
    // double xv = 0.03;
    // double yv = -1.09;
    // double zv = -1.88;

    if (((x-xv)*(x-xv)+(y-yv)*(y-yv)+(z-zv)*(z-zv)) < r*r)
    {
      //cell->SetParameter("Correction", 1.045);
      cell->SetParameter("eta", 0.037);
    }
    return cell;
  }
};

class TestRatStomach2D : public CxxTest::TestSuite
{
public:
  void TestSimulation()
  {
    ///// Input file
    string fname = "rat_scaffold_32_32_2_2D";
    HeartConfig::Instance()->Reset();
    HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach2D/"+fname, cp::media_type::Orthotropic);

    ///// Simulation settings
    int sim_dur = 90000; //ms
    int print_step = 1000; //ms
    HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
    HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-3);
    HeartConfig::Instance()->SetCapacitance(2.5);
    HeartConfig::Instance()->SetSimulationDuration(sim_dur);  //ms.
    HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.5,print_step);
    HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.01, 0.5, 0.01));//(0.03, 3.4, 1.0));

    ///// Output file/folder
    string out_path = "test_stomach2d_monodomain_"+fname+"_"+std::to_string(sim_dur)+"ms_"+std::to_string(print_step)+"ms";
    string out_add = "_imtiaz_test1";
    HeartConfig::Instance()->SetOutputDirectory(out_path+out_add);
    HeartConfig::Instance()->SetOutputFilenamePrefix("results");
    HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);

    ///// Output visualization options
    HeartConfig::Instance()->SetVisualizeWithCmgui(false);
    HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
    HeartConfig::Instance()->SetVisualizeWithVtk(true);

    ///// Cell factory
    ICCCellFactory cell_factory;

    ///// MonodomainProblem
    MonodomainProblem<2,3> monodomain_problem(&cell_factory);
    monodomain_problem.Initialise();
    monodomain_problem.Solve();
  }
};
