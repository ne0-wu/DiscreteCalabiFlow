#pragma once

#include <Eigen/Sparse>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include "Mesh.h"

class CalabiFlow
{
public:
	CalabiFlow(Mesh &input_mesh);

	void flow(int num_iterations = 100);

	Eigen::MatrixX2d &get_uv() { return uv; }

private:
	Mesh &mesh;

	OpenMesh::EProp<double> l;		// edge length
	OpenMesh::EProp<double> l_orig; // original edge length
	OpenMesh::HProp<double> angles; // angle opposite to halfedge

	Eigen::VectorXd u;		  // discrete conformal factor
	Eigen::VectorXd K;		  // Gaussian curvature
	Eigen::VectorXd K_target; // target Gaussian curvature

	Eigen::SparseMatrix<double> L; // Laplacian matrix

	Eigen::MatrixX2d uv; // uv coordinates

	void compute_edge_lengths();
	void compute_angles();
	void compute_gaussian_curvature();
	void compute_laplacian_matrix();

	void local_global(int num_iterations = 100); // generate uv from metric
};
