#include <iostream>
#include <stdlib.h> 
#include <fstream>
//#include "Colsamm.h"
typedef double Real;

using namespace std;
//using namespace ::_COLSAMM_;

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

	cout << delta  << '\n';
	ifstream ucircle;
	ucircle.open("/inputs/unit_circle.txt");
	Real a, b, c;
	string tmp;
	ucircle >> tmp;
	cout << tmp << '\n';
	ucircle >> tmp;
	cout << tmp << '\n';
	while (ucircle >> a && ucircle >> b && ucircle >> c)
	{
		cout << a << " " << b << " " << c << '\n';
	}

	return 0;

}