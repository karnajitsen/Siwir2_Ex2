#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Precision.h"
#include "constants.h"

using namespace std;

class Vertex {
	private:
		real m_x;
		real m_y;
		size_t m_index;
		vector<Vertex*> m_neighbours;
		
	public:
		Vertex(size_t index, real x, real y);
		real getX() const;
		real getY() const;
		void addNeighbour(Vertex* neighbour);
		void print();
		void printNeighbours();
		size_t getIndex() const;
		VertexComp compareWith(Vertex &vert);
};
