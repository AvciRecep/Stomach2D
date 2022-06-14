
#ifndef TestGenerateSheetFibre_HPP_
#define TestGenerateSheetFibre_HPP_

#include <cxxtest/TestSuite.h>

#include "UblasIncludes.hpp"
/* This class represents the mesh internally. */
#include "TetrahedralMesh.hpp"
/* These are used to specify boundary conditions for the PDEs. */
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
/* This class helps us deal with output files. */
#include "OutputFileHandler.hpp"
/* The following header must be included in every test that uses PETSc. Note that it
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"

// My includes
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include <../src/Eigen/Core>
#include <../src/Eigen/Geometry>

#define PI 3.14159265

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::AngleAxisd;
using Eigen::Vector3d;
using Eigen::Quaterniond;
#include "VtkMeshWriter.hpp"

std::string file_name = "rat_scaffold_32_32_2_2D";
std::string laplace_path = "/hpc/ravc486/Projects/CHASTE/chaste-hpc/chaste_ravc486/testoutput/";

class TestGenerateFibreSheetNormal : public CxxTest::TestSuite
{
private:

public:
    void TestCalculateFibreSheetNormal()
    {
        // Read mesh files
        std::string full_path = "projects/mesh/Stomach2D/" + file_name;
        TrianglesMeshReader<2,3> mesh_reader(full_path);
        // Now declare a tetrahedral mesh with the same dimensions... //
        TetrahedralMesh<2,3> mesh;
        // ... and construct the mesh using the mesh reader. //
        mesh.ConstructFromMeshReader(mesh_reader);

        full_path = laplace_path + "test_laplace_longi_" + file_name + "/" + file_name + "_laplace_longi_grad.txt";
        std::ifstream grad_longi_ijk(full_path);
        if (!grad_longi_ijk)
        {
            cout << "There was a problem opening laplace longitudinal gradient for reading " << endl;
        }

        full_path = laplace_path + "test_laplace_circum_" + file_name + "/" + file_name + "_laplace_circum_grad.txt";
        std::ifstream grad_circum_ijk(full_path);
        if (!grad_circum_ijk)
        {
            cout << "There was a problem opening laplace circumferential gradient for reading " << endl;
        }

        full_path = "test_laplace_ortho_"+file_name;
        OutputFileHandler output_file_handler(full_path);
        full_path = file_name+".ortho";
        out_stream p_file = output_file_handler.OpenOutputFile(full_path);

        // Read the first line in element file to get number of element and write it to the ortho file
        full_path = "projects/mesh/Stomach2D/" + file_name + ".ele";
        std::ifstream inElem(full_path);
        if (!inElem)
        {
            cout << "There was a problem opening faces for reading " << endl;
        }
        std::string line;
        if(!std::getline(inElem, line))
        {
            cout << "Error reading file line" << endl;
        }
        unsigned int numElem, dummy1, dummy2;
        stringstream numElemLine(line);
        numElemLine >> numElem >> dummy1 >> dummy2;
        (*p_file) << numElem << "\n";

        // Compute and write fibre, sheet, normal vectors
        double x, y, z;
        Vector3d fibre;
        Vector3d sheet;
        Vector3d normal;

        std::vector<c_vector<double, 3u> > fibre_directions;
        std::vector<c_vector<double, 3u> > sheet_directions;
        std::vector<c_vector<double, 3u> > normal_directions;
        while (std::getline(grad_longi_ijk, line))
        {
            std::stringstream fibreStream(line);
            fibreStream >> x >> y >> z;
            fibre(0) = x;
            fibre(1) = y;
            fibre(2) = z;
            fibre = fibre.normalized();
            if (((fibre.array() != fibre.array())).all())
            {
                fibre(0) = 1;
                fibre(1) = 0;
                fibre(2) = 0;
            }

            c_vector<double, 3u> fibre_direction;
            fibre_direction[0] = fibre(0);
            fibre_direction[1] = fibre(1);
            fibre_direction[2] = fibre(2);
            fibre_directions.push_back(fibre_direction);

            // sheet directions //
            std::getline(grad_circum_ijk, line);
            std::stringstream sheetStream(line);
            sheetStream >> x >> y >> z;
            sheet(0) = x;
            sheet(1) = y;
            sheet(2) = z;
            sheet = sheet.normalized();
            if (((sheet.array() != sheet.array())).all())
            {
                sheet(0) = 0;
                sheet(1) = 0;
                sheet(2) = 1;
            }
            double dot_product = fibre.dot(sheet);
            Vector3d temp;
            if( fabs(dot_product) != 0)
            {
                // Fix orthoganality //
                temp = fibre.cross(sheet);
                sheet = temp.cross(fibre).normalized();
            }
            c_vector<double, 3u> sheet_direction;
            sheet_direction[0] = sheet(0);
            sheet_direction[1] = sheet(1);
            sheet_direction[2] = sheet(2);
            sheet_directions.push_back(sheet_direction);

            // normal directions //
            normal = fibre.cross(sheet);
            c_vector<double, 3u> normal_direction;
            normal_direction[0] = normal(0);
            normal_direction[1] = normal(1);
            normal_direction[2] = normal(2);
            normal_directions.push_back(normal_direction);

            // Write CHASTE ortho file //
            if (((fibre.array() != fibre.array())).all() ||((sheet.array() != sheet.array())).all() )
            {
                fibre(0) = 1;
                fibre(1) = 0;
                fibre(2) = 0;
                sheet(0) = 0;
                sheet(1) = 1;
                sheet(2) = 0;
                normal(0) = 0;
                normal(1) = 0;
                normal(2) = 1;
            }
            (*p_file) << fibre(0) << " " << fibre(1) << " " << fibre(2) << " "
                      << sheet(0) << " " << sheet(1) << " " << sheet(2) << " "
                      << normal(0) << " " << normal(1) << " " << normal(2) << "\n";
        }
        p_file->close();

        // Write mesh with fibres as a vtk file
        full_path = "test_laplace_ortho_" + file_name;
        VtkMeshWriter<2u, 3u> mesh_writer(full_path, "mesh", false);
        mesh_writer.AddCellData("Fibre Direction", fibre_directions);
        mesh_writer.AddCellData("Sheet Direction", sheet_directions);
        mesh_writer.AddCellData("Normal Direction", normal_directions);
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
};
#endif
