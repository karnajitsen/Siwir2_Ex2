#include<iostream>
#include<fstream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
//#include "Colsamm.h"

typedef double Real;
typedef int unint;
#define ALLIGNMENT 64

using namespace std;
//using namespace ::_COLSAMM_;

struct childnode{
      double xcord, ycord;
};


struct node{
     childnode * neighbours;
     double xcord, ycord;
};

int main(int argc, char** argv)
{

	cout << "In the programm" << '\n';
	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}

	Real delta = atof(argv[1]);
	Real eps = atof(argv[2]);
    node* nodes = nullptr;

	cout << delta  << '\n';
	ifstream ucircle; 
	Real a, b, c;
    int d,e,f;
    size_t novert;
    string tmp,tmp1;
     //ucircle >> c;

    ucircle.open("./inputs/unit_circle.txt");
    getline(ucircle,tmp);

    novert = stoi(tmp.substr(0,tmp.find(" ")-1));

    nodes = (node*) memalign(ALLIGNMENT, novert*sizeof(node));

    getline(ucircle,tmp);

    for (size_t i = 0; i<novert ; i++)
	{
        ucircle >> a;
        ucircle >> b;
        ucircle >> c;
        nodes[i].xcord = b;
        nodes[i].ycord = c;
        //cout << i << " " << nodes[i].xcord << " " << nodes[i].ycord << '\n';
	}

    getline(ucircle,tmp1);
    //getline(ucircle,tmp1);
    //cout << tmp1 << '\n';

    getline(ucircle,tmp1);
     //  cout << tmp1 << '\n';

    while (ucircle >> d && ucircle >> e && ucircle >> f)
    {
        //cout << d << " " << e << " " << f << '\n';
    }

	return 0;

}
