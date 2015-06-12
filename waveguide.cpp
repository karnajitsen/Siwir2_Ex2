#include<iostream>
#include<fstream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#include <map>
#include "myColsamm\Source\Colsamm.h"

typedef double Real;
typedef int unint;
#define ALLIGNMENT 64

using namespace std;
using namespace ::_COLSAMM_;

struct node{
	double xcord, ycord, value;
	size_t vertno;
	  bool boundary;
};

struct graph{
	 
     std::map<size_t,node> nodes;
     size_t length;
};

struct triang{
	size_t vertex[3];
};

Real delta, eps;

graph* __restrict ugraphs = nullptr;
node* __restrict knodes = nullptr;
node* __restrict unodes = nullptr;
triang * __restrict tri = nullptr;
size_t novert, notriangle;

inline double kxy2(const double x, const double y)
{
	return (((100.0 + delta) * exp(-50.0*(x*x + y*y))) - 100.0);
}

inline void init()
{	
	string tmp;
	ifstream ucircle;
	Real a, b, c;
	size_t d, e, f;
	node * cn;

	ucircle.open("./inputs/unit_circle.txt");
	getline(ucircle, tmp);

	novert = stoi(tmp.substr(0, tmp.find(" ") - 1));

	ugraphs = (graph*)memalign(ALLIGNMENT, novert*sizeof(graph));
	unodes = (node*)memalign(ALLIGNMENT, novert*sizeof(node));

	getline(ucircle, tmp);

	for (size_t i = 0; i<novert; i++)
	{
		ucircle >> a;
		ucircle >> b;
		ucircle >> c;
		unodes[i].xcord = b;
		unodes[i].ycord = c;
		unodes[i].value = 0.0;
		unodes[i].vertno = i;
		if ((b*b + c*c) >= 1.0)  // Check for precision error ???
			unodes[i].boundary = true;
		else
			unodes[i].boundary = false;
		
		//cout << i << " " << nodes[i].xcord << " " << nodes[i].ycord << '\n';
	}

	getline(ucircle, tmp);
	getline(ucircle, tmp);
	
	for (notriangle = 1; ucircle >> d && ucircle >> e && ucircle >> f; notriangle++)
	{
		
		ugraphs[d].nodes.emplace(d, unodes[d]);
		ugraphs[d].nodes.emplace(e, unodes[e]);
		ugraphs[d].nodes.emplace(f, unodes[f]);
		ugraphs[e].nodes.emplace(e, unodes[e]);
		ugraphs[e].nodes.emplace(d, unodes[d]);
		ugraphs[e].nodes.emplace(f, unodes[f]);
		ugraphs[f].nodes.emplace(f, unodes[f]);
		ugraphs[f].nodes.emplace(e, unodes[e]);
		ugraphs[f].nodes.emplace(d, unodes[d]);
		
		//cout << d << " " << e << " " << f << '\n';
		tri[notriangle].vertex[0] = d;
		tri[notriangle].vertex[1] = e;
		tri[notriangle].vertex[2] = f;
	}

	knodes = (node*)memalign(ALLIGNMENT, novert*sizeof(node));

	for (size_t i; i < novert; i++)
	{
		knodes[i].xcord = unodes[i].xcord;
		knodes[i].ycord = unodes[i].ycord;
		knodes[i].value = kxy2(unodes[i].xcord, unodes[i].ycord);
	}


}

inline void createLocalMatrix(size_t a, size_t b, size_t c, std::vector<std::vector<double>>& localstiff, std::vector<std::vector<double>>& localmass)
{
	ELEMENTS::Triangle my_element;
	//std::vector< std::vector< double > > loc_stiff, loc_mass;
	std::vector<double> corners(6, 0.0);
	
	corners[0] = ugraphs[a].nodes.at(a).xcord; 
	corners[1] = ugraphs[a].nodes.at(a).ycord;
	corners[2] = ugraphs[b].nodes.at(b).xcord;
	corners[3] = ugraphs[b].nodes.at(b).ycord;
	corners[4] = ugraphs[c].nodes.at(c).xcord;
	corners[5] = ugraphs[c].nodes.at(c).ycord;
	// pass the corners to the finite element
	my_element(corners);

	localstiff =
		my_element.integrate(grad(v_()) * grad(w_()));

	localmass =
		my_element.integrate(grad(v_()) * grad(w_()));
}

inline void createGlobalMatrix()
{
	std::vector<std::vector<double>> localstiff, localmass;
	size_t a, b, c;
	for (size_t i = 0; i < notriangle; i++)
	{
		a = tri[notriangle].vertex[0];
		b = tri[notriangle].vertex[1];
		c = tri[notriangle].vertex[2];
		
		createLocalMatrix(a, b, c, localstiff, localmass);

		ugraphs[a].nodes.at(a).value += localstiff[0][0];
		ugraphs[a].nodes.at(b).value += localstiff[0][1];
		ugraphs[a].nodes.at(c).value += localstiff[0][2];

		ugraphs[b].nodes.at(a).value += localstiff[1][0];
		ugraphs[b].nodes.at(b).value += localstiff[1][1];
		ugraphs[b].nodes.at(c).value += localstiff[1][2];

		ugraphs[c].nodes.at(a).value += localstiff[2][0];
		ugraphs[c].nodes.at(b).value += localstiff[2][1];
		ugraphs[c].nodes.at(c).value += localstiff[2][2];

	}
}


int main(int argc, char** argv)
{

	cout << "In the programm" << '\n';
	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	delta = atof(argv[1]);
	eps = atof(argv[2]);
    	
	init();
	createGlobalMatrix();

	return 0;

}
