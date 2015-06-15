#include "Vertex.h"

Vertex::Vertex(size_t index, real x, real y)
{
	m_index = index;
	m_x = x;
	m_y = y;
}

real Vertex::getX() const
{
	return m_x;
}

real Vertex::getY() const
{
	return m_y;
}

void Vertex::addNeighbour(Vertex* neighbour)
{
	for(Vertex* existingNeighbour : m_neighbours)
	{
		if(existingNeighbour->getIndex() == neighbour->getIndex())
		{
			return;
		}
	}
	
	m_neighbours.push_back(neighbour);
}

void Vertex::print()
{
	size_t tabwidth = 14;
	cout << "Vertex" << setw(tabwidth/2) << m_index << ":" << setw(tabwidth) << m_x << setw(tabwidth) << m_y << endl;
}

void Vertex::printNeighbours()
{
	size_t tabwidth = 10;
	cout << "Vertex" << setw(tabwidth-5) << m_index << " Neighbours:";
	
	for(Vertex* neighbour : m_neighbours)
	{
		cout << setw(tabwidth) << neighbour->getIndex();
	}
	
	cout << endl;
}

size_t Vertex::getIndex() const
{
	return m_index;
}

VertexComp Vertex::compareWith(Vertex &vert)
{
	if(m_index == vert.getIndex())
	{
		return IDENTICAL;
	}
	
	//real precision = 0.0001;
    //if(((abs(m_x - vert.getX())) < precision)
    //&& ((abs(m_y -  vert.getY())) < precision))
	if(m_x == vert.getX() && m_y == vert.getY())
	{
		return SAME_POSITION;
	}
	
	return DIFFERENT;
}
