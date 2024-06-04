#pragma once

#include <Eigen/Sparse>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include "Mesh.h"

class CalabiFlow
{
public:
	CalabiFlow(Mesh &input_mesh);

	void flow(int num_iterations = 100);

private:
	Mesh mesh;

	OpenMesh::EProp<double> l;		// edge length
	OpenMesh::EProp<double> l_orig; // eriginal edge length
	OpenMesh::HProp<double> angles; // angles opposite to halfedges

	Eigen::VectorXd u;		  // discrete conformal factor
	Eigen::VectorXd K;		  // Gaussian curvature
	Eigen::VectorXd K_target; // target Gaussian curvature

	Eigen::SparseMatrix<double> L; // Laplacian matrix

	void iterate_once();
	void compute_edge_lengths();
	void compute_angles();
	void compute_gaussian_curvature();
	void compute_laplacian_matrix();

	double energy();
};
