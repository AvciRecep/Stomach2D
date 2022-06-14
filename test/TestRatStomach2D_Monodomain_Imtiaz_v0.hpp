#include <cxxtest/TestSuite.h>
#include <cassert>

#include <set>

#include "BidomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
//#include "../src/ICCCBDerivedCa.hpp"
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

class ICCCellFactory : public AbstractCardiacCellFactory<2,3>
{

public:
    ICCCellFactory()
        :AbstractCardiacCellFactory<2,3>()
    {
    }
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        //unsigned index =  pNode->GetIndex();
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];
        cout << x << " - " << y << " - " << z << "\n";
        ChastePoint<3> centre(0.08,-1.14,-1.95);
        ChastePoint<3> radii (0.1,0.1,0.1);
        ChasteEllipsoid<3> ellipseRegion(centre, radii);
        ChastePoint<3> myPoint(x,y,z);
        Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(mpSolver, mpZeroStimulus);
        cell->SetParameter("eta", 0.045);
        if(ellipseRegion.DoesContain(myPoint))
        {
            cell->SetParameter("eta", 0.037);
        }
        return cell;
    }
};

class TestRatStomach2D : public CxxTest::TestSuite
{
public:
    void TestSimulation() //throw(Exception)
    {
        ///// Read and contruct mesh
        DistributedTetrahedralMesh<2,3> mesh;
        std::string meshFile = "projects/mesh/Stomach2D/rat_scaffold_64_64_2_2D";
        TrianglesMeshReader<2,3> mesh_reader(meshFile);
        mesh.ConstructFromMeshReader(mesh_reader);
        //HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach2D/rat_scaffold_16_16_2_2D");

        ///// Mesh check
        /*unsigned nElements = 0;
        nElements = mesh.GetNumLocalElements();
        TRACE("Number of elements: " << nElements);
        for (DistributedTetrahedralMesh<2,3>::ElementIterator iter = mesh.GetElementIteratorBegin(); iter != mesh.GetElementIteratorEnd();
             ++iter)
        {
            Node<3>* pNode = iter->GetNode(0);
            double x = pNode->rGetLocation()[0];
            double y = pNode->rGetLocation()[1];
            double z = pNode->rGetLocation()[2];
            cout << x << " - " << y << " - " << z << "\n";
        }*/

        ///// Define monodomain
        ICCCellFactory cell_factory;
        MonodomainProblem<2,3> monodomain_problem(&cell_factory);
	      monodomain_problem.SetMesh(&mesh);

        ///// Reset heart config
        HeartConfig::Instance()->Reset();

        ///// Output visualization options, we ask for meshalyzer and cmgui
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        ///// Simulation settings
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(30000);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,0.1,500);

        ///// Output file/folder
        HeartConfig::Instance()->SetOutputDirectory("rat_scaffold_64_64_2_2D_v0");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        ///// Set conductivity values
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.3, 0.3, 0.3));//(0.03, 3.4, 1.0));
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }
};
