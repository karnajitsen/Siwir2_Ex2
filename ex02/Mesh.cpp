#include "Mesh.h"



Mesh::Mesh(string filename)
{
	ifstream input;
	string line;
	
	input.open(filename, ios::in);
	if (input.is_open())
	{
		// extract the number of vertex of the file
		// assume the number is first word of first line
		getline(input,line);
		m_vertexCount = getNumberOfElements(line);
		
		// ignore next line
		getline(input,line);
		
		for(size_t i = 0; i < m_vertexCount; ++i)
		{
			getline(input,line);
			m_vertice.push_back(extractVertex(line));
		}
		
		// extract the number of polygons (tris) of the file
		// assume the number is first word of first line
		getline(input,line);
		m_polyCount = getNumberOfElements(line);
		
		// ignore next line
		getline(input,line);
		
		for(size_t i = 0; i < m_polyCount; ++i)
		{
			getline(input,line);
			m_polygons.push_back(extractPolygon(line));
		}
		
		input.close();
	}
}


size_t Mesh::getNumberOfElements(string line)
{
	size_t start = line.find_first_not_of(' ', 0);
	size_t end = line.find_first_of(' ', start);
	string count = line.substr (start, end);
	return stoi(count);
}


Vertex Mesh::extractVertex(string line)
{
	// first char other than spaces is the index of the vertex
	size_t startIndex = line.find_first_not_of(' ', 0);
	size_t endIndex = line.find_first_of(' ', startIndex);
	string index = line.substr(startIndex,endIndex);
	
	// this should be the x pos of the vertex
	size_t startX = line.find_first_not_of(' ', endIndex);
	size_t endX = line.find_first_of(' ', startX);
	string posX = line.substr(startX,endX);
	
	size_t startY =  line.find_first_not_of(' ', endX);
	size_t endY = line.find_first_of(' ', startY);
	string posY = line.substr(startY,endY);
	
	return Vertex(stod(index), stod(posX), stod(posY));
}


Polygon Mesh::extractPolygon(string line)
{

	// first char other than spaces is the index of the  first vertex
	size_t start1 = line.find_first_not_of(' ', 0);
	size_t end1 = line.find_first_of(' ', start1);
	string vert1 = line.substr(start1,end1);
	
	size_t start2 = line.find_first_not_of(' ', end1);
	size_t end2 = line.find_first_of(' ', start2);
	string vert2 = line.substr(start2,end2);
	
	size_t start3 = line.find_first_not_of(' ', end2);
	size_t end3 = line.find_first_of(' ', start3);
	string vert3 = line.substr(start3,end3);
	
	return Polygon(stoi(vert1), stoi(vert2), stoi(vert3));
}


void Mesh::buildNeighbours()
{
	for(size_t i = 0; i < m_vertexCount; ++i)
	{
		for(vector<Polygon>::iterator it = m_polygons.begin(); it != m_polygons.end(); it++)
		{
			if(it->containsVertex(i))
			{
				vector<size_t> other = it->getOther(i);
				for(vector<size_t>::iterator it2 = other.begin(); it2 != other.end(); it2++)
				{
					m_vertice[i].addNeighbour(getVertex(*it2));
				}
			}
		}
	}
}


void Mesh::printVertice()
{
	for(Vertex vert : m_vertice)
	{
		vert.print();
	}
}


void Mesh::printPolys()
{
	for(Polygon poly : m_polygons)
	{
		poly.print();
	}
}

void Mesh::printNeighbours()
{
	for(Vertex vert : m_vertice)
	{
		vert.printNeighbours();
	}
}

void Mesh::printMeta()
{
	cout << "The mesh has " << m_vertexCount << " Vertice and " << m_polyCount << " Polygons" << endl;
}


void Mesh::subdivide()
{
	size_t highestIndex = m_vertice[m_vertexCount-1].getIndex();
	vector<Polygon> subdivPolys;
	
	for(Polygon poly : m_polygons)
	{
		vector<size_t> indices = poly.getIndices();
		
		Vertex* vert1 = getVertex(indices[0]);
		Vertex* vert2 = getVertex(indices[1]);
		Vertex* vert3 = getVertex(indices[2]);
		
		Vertex newVert1 = Vertex(++highestIndex, (vert1->getX()+vert2->getX())/2.0, (vert1->getY()+vert2->getY())/2.0);
		Vertex newVert2 = Vertex(++highestIndex, (vert1->getX()+vert3->getX())/2.0, (vert1->getY()+vert3->getY())/2.0);
		Vertex newVert3 = Vertex(++highestIndex, (vert2->getX()+vert3->getX())/2.0, (vert2->getY()+vert3->getY())/2.0);
		
		subdivPolys.push_back(Polygon(vert1->getIndex(), newVert1.getIndex(), newVert2.getIndex()));
		subdivPolys.push_back(Polygon(vert2->getIndex(), newVert1.getIndex(), newVert3.getIndex()));
		subdivPolys.push_back(Polygon(vert3->getIndex(), newVert2.getIndex(), newVert3.getIndex()));
		subdivPolys.push_back(Polygon(newVert1.getIndex(), newVert2.getIndex(), newVert3.getIndex()));
		
		m_vertice.push_back(newVert1);
		m_vertice.push_back(newVert2);
		m_vertice.push_back(newVert3);
	}
	
	m_polyCount = subdivPolys.size();
	m_vertexCount = m_vertice.size();
	m_polygons = subdivPolys;
	removeDoubles();
}


Vertex* Mesh::getVertex(size_t index)
{
	for(size_t i = 0; i < m_vertexCount; ++i)
	{
		if(m_vertice[i].getIndex() == index)
		{
			return &m_vertice[i];
		}
	}
	
	return NULL;
}


void Mesh::removeDoubles()
{
	for(vector<Vertex>::iterator it0 = m_vertice.begin(); it0 != m_vertice.end(); it0++)
	{
		for(vector<Vertex>::iterator it = m_vertice.begin(); it != m_vertice.end(); it++)
		{
			if((it0->compareWith(*it)) == SAME_POSITION)
			{
				size_t index = it->getIndex();
				
				// remoe redundant vertex from list
				m_vertice.erase(it);
				
				// replace all occourences of the redundant vertex in polygons
				for(vector<Polygon>::iterator it2 = m_polygons.begin(); it2 != m_polygons.end(); it2++)
				{
					if(it2->containsVertex(index))
					{
						it2->replaceVertex(index, it0->getIndex());
					}
				}
				
				break;
			}
		}
	}
	m_vertexCount = m_vertice.size();
}


bool Mesh::writeFile(string filename)
{
	ofstream file(filename, ofstream::out);
	if (file.is_open())
	{
		file << " " << m_vertexCount << " vertices in the domain" << endl;
		file << " Index of vertex | x coordinate | y coordinate" << endl;
		
		for(Vertex vert : m_vertice)
		{
			file << setw(15) << vert.getIndex() << setw(14) << vert.getX() << setw(14) << vert.getY() << endl;
		}
		
		file << " " << m_polyCount << " faces in the domain" << endl;
		file << " index of vertex 0 | index of vertex 1 | index of vertex 2 " << endl;
		
		for(Polygon poly : m_polygons)
		{
			vector<size_t> indices = poly.getIndices();
			
			for(size_t index : indices)
			{
				file << setw(15) << index;
			}
			
			file << endl;
		}
		
		file.close();
		return true;
	}

	return false;
}


bool Mesh::calculate(real(*Function)(real, real, real), string filename, real delta) const
{
	ofstream file(filename, ofstream::out);
	if (file.is_open())
	{
		for(Vertex vert : m_vertice)
		{
			file << " " << vert.getX() << " " << vert.getY() << " " << Function(vert.getX(), vert.getY(), delta) << " " << endl;
		}
		
		file.close();
		return true;
	}
	
	return false;
}


void Mesh::subdivide(size_t level)
{
	for(size_t i = 0; i < level; ++i)
	{
		subdivide();
		printMeta();
	}
	
	buildNeighbours();
}



