#include <iostream>

#include <Eigen/IterativeLinearSolvers>

#include "CalabiFlow.h"

CalabiFlow::CalabiFlow(Mesh &input_mesh)
	: mesh(input_mesh),
	  l(mesh), l_orig(mesh), angles(mesh),
	  u(mesh.n_vertices()),
	  K(mesh.n_vertices()), K_target(mesh.n_vertices()),
	  L(mesh.n_vertices(), mesh.n_vertices())
{
	for (const auto &e : mesh.edges())
	{
		l_orig[e] = mesh.calc_edge_length(e);
		l[e] = l_orig[e];
	}

	for (const auto &v : mesh.vertices())
	{
		u[v.idx()] = 0.0;
		K_target[v.idx()] = 0.0;
	}
}

void CalabiFlow::flow(int num_iterations)
{
	compute_edge_lengths();
	compute_angles();
	compute_gaussian_curvature();
	std::cout << "Initial Energy: " << energy() << std::endl;

	for (int i = 0; i < num_iterations; ++i)
	{
		iterate_once();
		std::cout << "Iteration " << i << " Energy: " << energy() << std::endl;
	}
}

void CalabiFlow::iterate_once()
{
	compute_edge_lengths();
	compute_angles();
	compute_gaussian_curvature();
	compute_laplacian_matrix();

	Eigen::VectorXd b = K_target - K;

	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	solver.compute(L);
	Eigen::VectorXd delta_u = solver.solve(b);

	u = u + 0.5 * delta_u;
}

void CalabiFlow::compute_edge_lengths()
{
	for (auto &e : mesh.edges())
	{
		auto fr = e.halfedge().from();
		auto to = e.halfedge().to();
		l[e] = l_orig[e] * exp(u[fr.idx()] + u[to.idx()]);
	}
}

void CalabiFlow::compute_angles()
{
	for (const auto &he : mesh.halfedges())
	{
		if (he.is_boundary())
			continue;

		auto edge = he.edge();
		auto prev = he.prev().edge();
		auto next = he.next().edge();
		angles[he] = acos((l[prev] * l[prev] + l[next] * l[next] - l[edge] * l[edge]) / (2.0 * l[prev] * l[next]));
	}
}

void CalabiFlow::compute_gaussian_curvature()
{
	for (auto &v : mesh.vertices())
	{
		if (v.is_boundary())
			K[v.idx()] = M_PI;
		else
			K[v.idx()] = 2.0 * M_PI;

		for (const auto &he : mesh.voh_range(v))
		{
			if (he.is_boundary())
				continue;

			K[v.idx()] -= angles[he.next()];
		}
	}
}

void CalabiFlow::compute_laplacian_matrix()
{
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(4 * mesh.n_halfedges());

	for (const auto &he : mesh.halfedges())
	{
		if (he.is_boundary())
			continue;

		auto fr = he.from();
		auto to = he.to();

		auto cot = 1.0 / tan(angles[he]);

		triplets.push_back({fr.idx(), to.idx(), -cot});
		triplets.push_back({to.idx(), fr.idx(), -cot});
		triplets.push_back({fr.idx(), fr.idx(), cot});
		triplets.push_back({to.idx(), to.idx(), cot});
	}

	L.setFromTriplets(triplets.begin(), triplets.end());
}

double CalabiFlow::energy()
{
	compute_edge_lengths();
	compute_angles();
	compute_gaussian_curvature();
	return (K - K_target).norm();
}