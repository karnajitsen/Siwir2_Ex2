#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "Vertex.h"
#include "Polygon.h"

using namespace std;

class Mesh {
	private:
		vector<Vertex> m_vertice;
		vector<Polygon> m_polygons;
		size_t m_vertexCount;
		size_t m_polyCount;
		
		size_t getNumberOfElements(string line);
		Vertex extractVertex(string line);
		Polygon extractPolygon(string line);
		void buildNeighbours();
		Vertex* getVertex(size_t index);
		void removeDoubles();
		void subdivide();
		
		
	public:
		Mesh(string filename);
		void printVertice();
		void printPolys();
		void printNeighbours();
		void printMeta();
		void subdivide(size_t level);
		bool writeFile(string filename);
		bool calculate(real(*Function)(real, real, real), string filename, real delta) const;

};
