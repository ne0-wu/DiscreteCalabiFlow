#include <iostream>

#include <Eigen/IterativeLinearSolvers>

#include "CalabiFlow.h"

CalabiFlow::CalabiFlow(Mesh &input_mesh)
	: mesh(input_mesh),
	  l(mesh), l_orig(mesh), angles(mesh),
	  u(mesh.n_vertices()),
	  K(mesh.n_vertices()), K_target(mesh.n_vertices()),
	  L(mesh.n_vertices(), mesh.n_vertices()),
	  uv(mesh.n_vertices(), 2)
{
	for (const auto &e : mesh.edges())
		l_orig[e] = mesh.calc_edge_length(e);

	// Count the number of boundary vertices
	OpenMesh::SmartHalfedgeHandle first_boundary_he;
	for (const auto &he : mesh.halfedges())
		if (he.is_boundary())
		{
			first_boundary_he = he;
			break;
		}
	int boundary_size = 0;
	for (auto he = first_boundary_he; he.is_valid(); he = he.next())
	{
		boundary_size++;
		if (he.next() == first_boundary_he)
			break;
	}

	// Initialize the discrete conformal factor and the target Gaussian curvature
	for (const auto &v : mesh.vertices())
	{
		u[v.idx()] = 0.0;
		if (v.is_boundary())
			K_target[v.idx()] = 2 * M_PI / boundary_size; // Gaussian curvature on the boundary of a disk
		else
			K_target[v.idx()] = 0.0;
	}
}

void CalabiFlow::flow(int num_iterations)
{
	compute_edge_lengths();
	compute_angles();
	compute_gaussian_curvature();

	double energy_prev = (K - K_target).squaredNorm();
	std::cout << "Initial Energy: " << energy_prev << std::endl;

	for (int iter = 1; iter <= num_iterations; iter++)
	{
		compute_laplacian_matrix();

		Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
		solver.compute(L * L);
		if (solver.info() != Eigen::Success)
		{
			std::cerr << "Decomposition failed!" << std::endl;
			exit(EXIT_FAILURE);
		}

		auto b = K_target - K;
		Eigen::VectorXd delta_u = solver.solve(L * b);

		u = u + 0.5 * delta_u;

		compute_edge_lengths();
		compute_angles();
		compute_gaussian_curvature();

		double energy = (K - K_target).squaredNorm();

		std::cout << "Iteration: " << iter
				  << " \t Energy: " << pow((K - K_target).norm(), 2)
				  << " \t Max Error: " << (K - K_target).cwiseAbs().maxCoeff()
				  << " \t Relative Energy Change: " << abs(energy_prev - energy) / energy
				  << std::endl;

		if (abs(energy_prev - energy) / energy < 1e-3)
		{
			std::cout << "Converged!" << std::endl;
			break;
		}
		energy_prev = energy;
	}

	local_global();
}

void CalabiFlow::compute_edge_lengths()
{
	for (const auto &e : mesh.edges())
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

struct SVD22
{
	// For an input matrix A, the SVD decomposition is given by A = U * S * V^T
	Eigen::Matrix2d U;
	Eigen::Matrix2d V;
	Eigen::Matrix2d S;

	SVD22(Eigen::Matrix2d &A)
	{
		// This algorithm actually provides the signed SVD decomposition
		// where sigma 2 could be negative
		double E = (A(0, 0) + A(1, 1)) / 2,
			   F = (A(0, 0) - A(1, 1)) / 2,
			   G = (A(1, 0) + A(0, 1)) / 2,
			   H = (A(1, 0) - A(0, 1)) / 2;

		double Q = sqrt(E * E + H * H), R = sqrt(F * F + G * G);

		double S1 = Q + R, S2 = Q - R;

		double T1 = atan2(G, F), T2 = atan2(H, E);

		double theta = (T2 - T1) / 2, phi = (T2 + T1) / 2;

		U = Eigen::Rotation2Dd(phi).toRotationMatrix();
		V = Eigen::Rotation2Dd(theta).toRotationMatrix().transpose();
		S << S1, 0.0, 0.0, S2;
	}

	Eigen::Matrix2d best_fit_ARAP()
	{
		return U * V.transpose();
	}

	Eigen::Matrix2d best_fit_ASAP()
	{
		double meanS = (S(0, 0) + S(1, 1)) / 2;
		Eigen::Matrix2d S;
		S << meanS, 0.0, 0.0, meanS;
		return U * S * V.transpose();
	}
};

void CalabiFlow::local_global(int num_iterations)
{
	// Tutte parameterization as initial guess
	// =======================================

	// Compute the Laplacian matrix
	compute_laplacian_matrix();

	// Fix the boundary vertices to the unit circle
	OpenMesh::SmartHalfedgeHandle first_boundary_he;
	for (const auto &he : mesh.halfedges())
		if (he.is_boundary())
		{
			first_boundary_he = he;
			break;
		}
	std::vector<int> boundary_indices;
	for (auto he = first_boundary_he; he.is_valid(); he = he.next())
	{
		boundary_indices.push_back(he.from().idx());
		if (he.next() == first_boundary_he)
			break;
	}
	Eigen::MatrixX2d b(mesh.n_vertices(), 2);
	b.setZero();
	for (int i = 0; i < boundary_indices.size(); i++)
		b.row(boundary_indices[i]) = Eigen::Vector2d(cos(2 * M_PI * i / boundary_indices.size()),
													 sin(2 * M_PI * i / boundary_indices.size()));

	// Construct the matrix A for the linear system
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(L.nonZeros());

	for (const auto &bdr : boundary_indices)
		triplets.push_back({bdr, bdr, 1.0});

	for (int k = 0; k < L.outerSize(); k++)
	{
		if (mesh.is_boundary(mesh.vertex_handle(k)))
			continue;

		for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
			if (L.Flags & Eigen::RowMajorBit)
				triplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value())); // Row major
			else
				triplets.push_back(Eigen::Triplet<double>(it.col(), it.row(), it.value())); // Column major
	}

	Eigen::SparseMatrix<double> A(mesh.n_vertices(), mesh.n_vertices());
	A.setFromTriplets(triplets.begin(), triplets.end());

	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	if (solver.info() != Eigen::Success)
	{
		std::cerr << "Decomposition failed!" << std::endl;
		exit(EXIT_FAILURE);
	}

	uv = solver.solve(b);

	// Local-Global ARAP parameterization
	// ==================================

	// Directly flatten the triangle pieces to 2D
	std::vector<Eigen::Matrix2d> tri_x(mesh.n_faces());
	for (const auto &f : mesh.faces())
	{
		auto he = f.halfedge();
		tri_x[f.idx()] << l[he.edge()], 0.0,
			l[he.prev().edge()] * cos(angles[he.next()]), l[he.prev().edge()] * sin(angles[he.next()]);
	}

	solver.compute(L);

	for (int iter = 1; iter <= num_iterations; iter++)
	{
		// Local step
		b.setZero();

		for (const auto &f : mesh.faces())
		{
			auto he01 = f.halfedge();
			int v0 = he01.from().idx(),
				v1 = he01.to().idx(),
				v2 = he01.next().to().idx();

			Eigen::Matrix2d tri_uv;
			tri_uv.row(0) = uv.row(v1) - uv.row(v0);
			tri_uv.row(1) = uv.row(v2) - uv.row(v0);

			Eigen::Matrix2d J = tri_x[f.idx()].inverse() * tri_uv;
			Eigen::Matrix2d R = SVD22(J).best_fit_ARAP();

			Eigen::Matrix<double, 1, 2> x0, x1, x2;
			x0.setZero();
			x1 = tri_x[f.idx()].row(0);
			x2 = tri_x[f.idx()].row(1);

			auto he12 = he01.next();
			auto he20 = he12.next();

			auto cot = [](auto angle)
			{ return 1.0 / tan(angle); };

			b.row(v0) += cot(angles[he01]) * (x0 - x1) * R + cot(angles[he20]) * (x0 - x2) * R;
			b.row(v1) += cot(angles[he12]) * (x1 - x2) * R + cot(angles[he01]) * (x1 - x0) * R;
			b.row(v2) += cot(angles[he20]) * (x2 - x0) * R + cot(angles[he12]) * (x2 - x1) * R;
		}

		// Global step
		uv = solver.solve(-b);
	}
}