#pragma once
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include "Precision.h"

using namespace std;

// just tris for now ;)
class Polygon {
	private:
		vector<size_t> m_indices;
		
	public:
		Polygon(size_t vert1, size_t vert2, size_t vert3);
		void print();
		void print(size_t number);
		bool containsVertex(size_t index);
		vector<size_t> getOther(size_t index);
		vector<size_t> getIndices();
		bool replaceVertex(size_t original, size_t replacement);
};
