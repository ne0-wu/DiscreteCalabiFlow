#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include "Mesh.h"
#include "CalabiFlow.h"

#include "DrawMesh.h"

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
	calabi_flow.flow();
	auto uv = calabi_flow.get_uv();

	// Save the result to file
	Mesh output_mesh = input_mesh;
	for (const auto &v : output_mesh.vertices())
		output_mesh.set_point(v, Mesh::Point(uv(v.idx(), 0), uv(v.idx(), 1), 0.0));

	if (!OpenMesh::IO::write_mesh(output_mesh, "output.obj"))
	{
		std::cerr << "Error writing mesh to file output.obj" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Draw the mesh
	MeshDrawer drawer(2000, 2000);

	std::vector<std::pair<double, double>> vertices;
	vertices.reserve(input_mesh.n_vertices());
	for (const auto &v : input_mesh.vertices())
		vertices.push_back({uv(v.idx(), 0), uv(v.idx(), 1)});

	std::vector<int> faces;
	faces.reserve(input_mesh.n_faces());
	for (const auto &f : input_mesh.faces())
	{
		for (const auto &v : input_mesh.fv_range(f))
			faces.push_back(v.idx());
	}

	drawer.draw_mesh(vertices, faces, std::vector<int>(input_mesh.n_faces(), 3));
	drawer.save_image("output.png");

	return 0;
}