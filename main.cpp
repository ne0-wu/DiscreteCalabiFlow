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
	if (!OpenMesh::IO::read_mesh(input_mesh, "data/models/cathead.obj"))
	{
		std::cerr << "Error loading mesh from file input.obj" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Discrete Calabi Flow
	CalabiFlow calabi_flow(input_mesh);
	calabi_flow.flow(100);

	return 0;
}