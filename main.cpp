#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "Mesh.h"
#include "CalabiFlow.h"

int main()
{
	// Load mesh from file
	Mesh input_mesh;
	if (!OpenMesh::IO::read_mesh(input_mesh, "data/models/camelhead.obj"))
	{
		std::cerr << "Error loading mesh from file input.obj" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Discrete Calabi Flow
	CalabiFlow calabi_flow(input_mesh);
	calabi_flow.flow();

	// Save the result to file
	auto uv = calabi_flow.get_uv();
	Mesh output_mesh = input_mesh;
	for (const auto &v : output_mesh.vertices())
		output_mesh.set_point(v, Mesh::Point(uv(v.idx(), 0), uv(v.idx(), 1), 0.0));

	if (!OpenMesh::IO::write_mesh(output_mesh, "output.obj"))
	{
		std::cerr << "Error writing mesh to file output.obj" << std::endl;
		exit(EXIT_FAILURE);
	}

	return 0;
}