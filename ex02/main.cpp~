#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Mesh.h"

using namespace std;

real HereinCoefficient(real x, real y, real delta)
{
	return (100.0 + delta) * exp(-50*(x*x+y*y)) - 100;
}

int main(int args, char **argv)
{
	if(args == 1)
	{
		cout << "are you sure you know what you're doing?" << endl;
		return 0;
	}
	else if(args < 3)
	{
		cout << "still something missing =/" << endl;
		return 0;
	}
	
	size_t delta = atoi(argv[1]);
	//size_t eta = atoi(argv[2]);
	size_t subdiv = 0;
	
	if(args == 4)
	{
		subdiv = atoi(argv[3]);
	}
	
	
	Mesh circle = Mesh("./inputs/unit_circle.txt");
	circle.printMeta();
	//circle.printNeighbours();
	
	circle.subdivide(subdiv);
	circle.writeFile("circle");
	circle.calculate(HereinCoefficient, "output", delta);
	
	
	return 0;
}
