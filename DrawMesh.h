#pragma once

#include <vector>
#include <cmath>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

class MeshDrawer
{
public:
	MeshDrawer(int width, int height)
		: width(width), height(height)
	{
		image = std::vector<unsigned char>(width * height * 3, 255); // Initialize with white color
	}

	void draw_mesh(const std::vector<std::pair<double, double>> &vertices, const std::vector<int> &faces, const std::vector<int> &faceSizes)
	{
		auto transformed_vertices = transform_vertices(vertices);

		size_t index = 0;
		for (const auto &faceSize : faceSizes)
		{
			std::vector<int> face(faceSize);
			for (int i = 0; i < faceSize; ++i)
			{
				face[i] = faces[index++];
			}
			for (size_t i = 0; i < face.size(); ++i)
			{
				int v1 = face[i];
				int v2 = face[(i + 1) % face.size()];
				draw_line(transformed_vertices[v1], transformed_vertices[v2]);
			}
		}
	}

	void save_image(const char *filename)
	{
		stbi_write_png(filename, width, height, 3, image.data(), width * 3);
	}

private:
	int width, height;
	std::vector<unsigned char> image;

	std::vector<std::pair<double, double>> transform_vertices(const std::vector<std::pair<double, double>> &vertices)
	{
		double min_x = vertices[0].first, max_x = vertices[0].first;
		double min_y = vertices[0].second, max_y = vertices[0].second;

		for (const auto &vertex : vertices)
		{
			if (vertex.first < min_x)
				min_x = vertex.first;
			if (vertex.first > max_x)
				max_x = vertex.first;
			if (vertex.second < min_y)
				min_y = vertex.second;
			if (vertex.second > max_y)
				max_y = vertex.second;
		}

		double scale_x = width / (max_x - min_x);
		double scale_y = height / (max_y - min_y);
		double scale = std::min(scale_x, scale_y);

		double translate_x = (width - (max_x - min_x) * scale) / 2;
		double translate_y = (height - (max_y - min_y) * scale) / 2;

		std::vector<std::pair<double, double>> transformed_vertices;
		for (const auto &vertex : vertices)
		{
			double x = (vertex.first - min_x) * scale + translate_x;
			double y = (vertex.second - min_y) * scale + translate_y;
			transformed_vertices.emplace_back(x, y);
		}

		return transformed_vertices;
	}

	void draw_line(const std::pair<double, double> &start, const std::pair<double, double> &end)
	{
		int x1 = static_cast<int>(start.first);
		int y1 = static_cast<int>(start.second);
		int x2 = static_cast<int>(end.first);
		int y2 = static_cast<int>(end.second);

		int dx = std::abs(x2 - x1);
		int dy = std::abs(y2 - y1);
		int sx = x1 < x2 ? 1 : -1;
		int sy = y1 < y2 ? 1 : -1;
		int err = (dx > dy ? dx : -dy) / 2;
		int e2;

		while (true)
		{
			set_pixel(x1, y1, 0, 0, 0); // Set pixel to black
			if (x1 == x2 && y1 == y2)
				break;
			e2 = err;
			if (e2 > -dx)
			{
				err -= dy;
				x1 += sx;
			}
			if (e2 < dy)
			{
				err += dx;
				y1 += sy;
			}
		}
	}

	void set_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
	{
		if (x >= 0 && x < width && y >= 0 && y < height)
		{
			int index = (y * width + x) * 3;
			image[index] = r;
			image[index + 1] = g;
			image[index + 2] = b;
		}
	}
};
