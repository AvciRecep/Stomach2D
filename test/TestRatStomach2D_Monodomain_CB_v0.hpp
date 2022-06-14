#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "../src/CellICCBioPhy.hpp"
#include "../src/DummyCell.hpp"
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
    ICCCellFactory()
        :AbstractCardiacCellFactory<2,3>()
    {
        ReadLaplaceFile();
    }

    void ReadLaplaceFile()
    {
        std::ifstream inLaplaceInfo("projects/mesh/Stomach2D/rat_scaffold_16_16_2_2D_laplace_longi.txt");
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
          PRINT_4_VARIABLES(x,y,z, counter);
          EXCEPTION("Coordinates not found in Laplace file");
      }

      CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
      cell->SetParameter("V_excitation", -55);
      cell->SetParameter("live_time", 12000);

      if (V_val > 45)
      {
          return new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
      }
      double r = 0.1;
      //cout << "x: " << x << " y: " << y << " z: " << z << "\n";
      if (((x-0.1)*(x-0.1)+(z+2.1)*(z+2.1)+(y+1.2)*(y+1.2)) < r*r)
      {
          cell->SetParameter("t_start", 0);
          cout << "I was here" << "\n";
      }
      else
      {
          cell->SetParameter("t_start", 600000);
      }

      cell->SetParameter("ode_time_step", 0.1);
      cell->SetParameter("IP3Par", 0.00069);
      return cell;
    }
};

class TestRatStomach2D : public CxxTest::TestSuite
{
public:
    void TestSimulation() //throw(Exception)
    {
        ///// Input file
        string fname = "rat_scaffold_16_16_2_2D";
        HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach2D/"+fname, cp::media_type::Orthotropic);

        ///// Simulation settings
        int sim_dur = 30000; // ms
        int write_freq = 500; //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(sim_dur);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,1,write_freq);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.01, 0.30, 0.30));//(0.03, 3.4, 1.0));

        ///// Output file/folder
        string out_path = "test_stomach2d_monodomain_CB_"+fname+"_"+std::to_string(sim_dur)+"ms_"+std::to_string(write_freq)+"ms";
        string out_add = "";
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
