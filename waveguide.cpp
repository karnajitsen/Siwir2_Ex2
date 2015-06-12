#include<iostream>
#include<fstream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#include "Colsamm.h"

typedef double Real;
typedef int unint;
#define ALLIGNMENT 64

using namespace std;
using namespace ::_COLSAMM_;

struct node{
      double xcord, ycord,value;
	  bool boundary;
};

struct graph{
     node *neighbours, center;
     size_t length;
};

struct triangle{
	node vertex[3];
};

Real delta, eps;

graph* __restrict ugraphs = nullptr;
node* __restrict knodes = nullptr;
node* __restrict unodes = nullptr;
triangle * __restrict tri = nullptr;

inline double kxy2(const double x, const double y)
{
	return (((100.0 + delta) * exp(-50.0*(x*x + y*y))) - 100.0);
}

inline void init()
{
	size_t novert, notriangle;
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
		if ((b*b + c*c) >= 1.0)  // Check for precision error ???
			unodes[i].boundary = true;
		else
			unodes[i].boundary = false;
		//unodes[i].length = 0;
		
		//cout << i << " " << nodes[i].xcord << " " << nodes[i].ycord << '\n';
	}

	getline(ucircle, tmp);
	//getline(ucircle,tmp1);
	//cout << tmp1 << '\n';

	getline(ucircle, tmp);
	//  cout << tmp1 << '\n';

	for (notriangle = 1; ucircle >> d && ucircle >> e && ucircle >> f; notriangle++)
	{
		node c1, c2, c3, c4, c5, c6;
		c1.xcord = unodes[e].xcord;
		c1.ycord = unodes[e].ycord;
		c1.value = 0.0;
		c2.xcord = unodes[f].xcord;
		c2.ycord = unodes[f].ycord;
		c2.value = 0.0;	
		ugraphs[d].center = unodes[d];
		ugraphs[d].length += 2;
		cn = (node*)memalign(ALLIGNMENT, ugraphs[d].length*sizeof(node));
		memcpy(cn, ugraphs[d].neighbours, (ugraphs[d].length - 2)*sizeof(node));
		cn[length - 2] = c1;
		cn[length - 1] = c2;
		memcpy(ugraphs[d].neighbours, cn, (ugraphs[d].length)*sizeof(node));
		free(cn);

		c3.xcord = unodes[d].xcord;
		c3.ycord = unodes[d].ycord;
		c3.value = 0.0;
		c4.xcord = unodes[f].xcord;
		c4.ycord = unodes[f].ycord;
		c4.value = 0.0;
		ugraphs[e].center = unodes[e];
		ugraphs[e].length += 2;
		cn = (node*)memalign(ALLIGNMENT, ugraphs[e].length*sizeof(node));
		memcpy(cn, ugraphs[e].neighbours, (ugraphs[e].length - 2)*sizeof(node));
		cn[length - 2] = c3;
		cn[length - 1] = c4;
		memcpy(ugraphs[e].neighbours, cn, (ugraphs[e].length)*sizeof(node));
		free(cn);

		c3.xcord = unodes[d].xcord;
		c3.ycord = unodes[d].ycord;
		c3.value = 0.0;
		c4.xcord = unodes[e].xcord;
		c4.ycord = unodes[e].ycord;
		c4.value = 0.0;
		ugraphs[f].center = unodes[f];
		ugraphs[f].length += 2;
		cn = (node*)memalign(ALLIGNMENT, ugraphs[f].length*sizeof(node));
		memcpy(cn, ugraphs[f].neighbours, (ugraphs[f].length - 2)*sizeof(node));
		cn[length - 2] = c3;
		cn[length - 1] = c4;
		memcpy(ugraphs[f].neighbours, cn, (ugraphs[f].length)*sizeof(node));
		free(cn);
		//cout << d << " " << e << " " << f << '\n';
		tri[notriangle].vertex[0] = unodes[d];
		tri[notriangle].vertex[1] = unodes[e];
		tri[notriangle].vertex[2] = unodes[f];
	}

	knodes = (childnode*)memalign(ALLIGNMENT, novert*sizeof(childnode));

	for (size_t i; i < novert, i++)
	{
		knodes[i].xcord = unodes[i].xcord;
		knodes[i].ycord = unodes[i].ycord;
		knodes[i].value = kxy2(unodes[i].xcord, unodes[i].ycord);
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
    

	return 0;

}
