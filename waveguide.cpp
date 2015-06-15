#include<iostream>
#include<fstream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#include<algorithm>
#include<vector>
#include <map>
#include "./myColsamm/Source/Colsamm.h"

#define ERRLIMIT 0.000000001

typedef double Real;
typedef int unint;
#define ALLIGNMENT 64
  
using namespace std;
using namespace ::_COLSAMM_;

struct node{
	Real xcord, ycord, stiffval, massval, uval, fval;
	size_t vertno;
	bool boundary;
};

struct graph{
	 
     std::map<size_t,node> nodes;
	 std::vector<int> index;
     size_t length;
};

struct triang{
	size_t vertex[3];
};

Real delta, eps;

graph* __restrict ugraphs = nullptr;
node* __restrict knodes = nullptr;
node* __restrict unodes = nullptr;
//node* __restrict fnodes = nullptr;
triang * __restrict tri = nullptr;
size_t novert, notriangle;

inline double kxy2(const Real x, const Real y)
{
	return (((100.0 + delta) * exp(-50.0*(x*x + y*y))) - 100.0);
}

inline void init()
{	
	string tmp;
	ifstream ucircle;
	Real a, b, c;
	int d, e, f;
	//node * cn;

	ucircle.open("./inputs/unit_circle.txt");
	getline(ucircle, tmp);

	novert = stoi(tmp.substr(0, tmp.find(" ") - 1));
	cout << "no of vertex = " << novert << '\n';

	
	unodes = new node[novert];
	getline(ucircle, tmp);
	//cout << tmp << '\n';
	for (size_t i = 0; i<novert; i++)
	{
		ucircle >> a;
		ucircle >> b;
		ucircle >> c;
		
		unodes[i].xcord = b;
		unodes[i].ycord = c;
		unodes[i].fval = 0.0;
		unodes[i].uval = 1.0;
		unodes[i].massval = 0.0;
		unodes[i].stiffval = 0.0;
		unodes[i].vertno = i;
		//cout << "2 " << a << " " << unodes[i].xcord << " " << unodes[i].ycord << '\n';
		if ((b*b + c*c) >= 1.0)  // Check for precision error ???
			unodes[i].boundary = true;
		else
			unodes[i].boundary = false;
		
		//cout << i << " " << nodes[i].xcord << " " << nodes[i].ycord << '\n';
	}

	getline(ucircle, tmp);	

	getline(ucircle, tmp);
	
	notriangle = stoi(tmp.substr(0, tmp.find(" ") - 1));;
	//ugraphs = (graph*)memalign(ALLIGNMENT, novert*sizeof(graph));
	ugraphs = new graph[novert];
	tri = new triang[notriangle];
	cout << "no of triangle = " << notriangle << '\n';
	getline(ucircle, tmp);
	
	for (size_t i = 0; ucircle >> d && ucircle >> e && ucircle >> f; i++)
	{
		ugraphs[d].nodes.emplace(d, unodes[d]);
		ugraphs[d].nodes.emplace(e, unodes[e]);
		ugraphs[d].nodes.emplace(f, unodes[f]);
		//cout << "77777" << '\n';
		if (std::find(ugraphs[d].index.begin(), ugraphs[d].index.end(), d) == ugraphs[d].index.end())
					ugraphs[d].index.emplace_back(d);
		
		if (std::find(ugraphs[d].index.begin(), ugraphs[d].index.end(), e) == ugraphs[d].index.end())
					ugraphs[d].index.emplace_back(e);
		
		if (std::find(ugraphs[d].index.begin(), ugraphs[d].index.end(), f) == ugraphs[d].index.end())
					ugraphs[d].index.emplace_back(f);
		
		ugraphs[e].nodes.emplace(e, unodes[e]);
		ugraphs[e].nodes.emplace(d, unodes[d]);
		ugraphs[e].nodes.emplace(f, unodes[f]);

		if (std::find(ugraphs[e].index.begin(), ugraphs[e].index.end(), e) == ugraphs[e].index.end())
			ugraphs[e].index.emplace_back(e);

		if (std::find(ugraphs[e].index.begin(), ugraphs[e].index.end(), d) == ugraphs[e].index.end())
			ugraphs[e].index.emplace_back(d);		

		if (std::find(ugraphs[e].index.begin(), ugraphs[e].index.end(), f) == ugraphs[e].index.end())
			ugraphs[e].index.emplace_back(f);

		ugraphs[f].nodes.emplace(f, unodes[f]);
		ugraphs[f].nodes.emplace(e, unodes[e]);
		ugraphs[f].nodes.emplace(d, unodes[d]);

		if (std::find(ugraphs[f].index.begin(), ugraphs[f].index.end(), f) == ugraphs[f].index.end())
			ugraphs[f].index.emplace_back(f);

		if (std::find(ugraphs[f].index.begin(), ugraphs[f].index.end(), d) == ugraphs[f].index.end())
			ugraphs[f].index.emplace_back(d);

		if (std::find(ugraphs[f].index.begin(), ugraphs[f].index.end(), e) == ugraphs[f].index.end())
			ugraphs[f].index.emplace_back(e);	
		
		//cout << d << " " << e << " " << f << '\n';
		tri[i].vertex[0] = d;
		tri[i].vertex[1] = e;
		tri[i].vertex[2] = f;
		
	}

	knodes = new node[novert];
	for (size_t i=0; i < novert; i++)
	{
		knodes[i].xcord = unodes[i].xcord;
		knodes[i].ycord = unodes[i].ycord;
		knodes[i].uval = kxy2(unodes[i].xcord, unodes[i].ycord);
		std::sort(ugraphs[i].index.begin(), ugraphs[i].index.end());

		/*for (int j = 0; j < ugraphs[i].index.size(); j++)
		{
			cout << " " << ugraphs[i].index[j];
		}
		cout << '\n';*/
	}
	//cout << "222" << '\n';

}

inline void createLocalMatrix(size_t a, size_t b, size_t c, std::vector<std::vector<double>>& localstiff, std::vector<std::vector<double>>& localmass)
{
	ELEMENTS::Triangle my_element;

	//std::vector< std::vector< double > > loc_stiff, loc_mass;
	std::vector<double> corners(6, 0.0);
	//cout << "1 " << " " << b << " " << c << '\n';
	corners[0] = ugraphs[a].nodes.at(a).xcord; 
	corners[1] = ugraphs[a].nodes.at(a).ycord;
	//cout << "2 " << a << " " << b << " " << c << '\n';
	corners[2] = ugraphs[b].nodes.at(b).xcord;
	corners[3] = ugraphs[b].nodes.at(b).ycord;
	//cout << "3 " <<  a << " " << b << " " << c << '\n';
	corners[4] = ugraphs[c].nodes.at(c).xcord;
	corners[5] = ugraphs[c].nodes.at(c).ycord;
	// pass the corners to the finite element
	//cout << a << " " << b << " " << c << '\n';
	my_element(corners);

	localstiff =
		my_element.integrate(grad(v_()) * grad(w_()) - func<double>(kxy2) * v_() * w_());// -my_element.integrate(func<double>(kxy2) * v_() * w_());

	localmass =
		my_element.integrate(v_() * w_());
}

inline void createGlobalMatrix()
{
	std::vector<std::vector<double>> localstiff, localmass;
	size_t a, b, c;
	
	for (size_t i = 0; i < notriangle; i++)
	{
		//cout << "333" << '\n';
		a = tri[i].vertex[0];
		b = tri[i].vertex[1];
		c = tri[i].vertex[2];
		//cout << "1 " << " " << b << " " << c << '\n';
		createLocalMatrix(a, b, c, localstiff, localmass);
		//cout << "555" << '\n';
		ugraphs[a].nodes.at(a).stiffval += localstiff[0][0];
		ugraphs[a].nodes.at(b).stiffval += localstiff[0][1];
		ugraphs[a].nodes.at(c).stiffval += localstiff[0][2];
		//cout << "666" << '\n';
		ugraphs[b].nodes.at(a).stiffval += localstiff[1][0];
		ugraphs[b].nodes.at(b).stiffval += localstiff[1][1];
		ugraphs[b].nodes.at(c).stiffval += localstiff[1][2];
		//cout << "777" << '\n';
		ugraphs[c].nodes.at(a).stiffval += localstiff[2][0];
		ugraphs[c].nodes.at(b).stiffval += localstiff[2][1];
		ugraphs[c].nodes.at(c).stiffval += localstiff[2][2];
		//cout << "88" << '\n';
		ugraphs[a].nodes.at(a).massval += localmass[0][0];
		ugraphs[a].nodes.at(b).massval += localmass[0][1];
		ugraphs[a].nodes.at(c).massval += localmass[0][2];

		ugraphs[b].nodes.at(a).massval += localmass[1][0];
		ugraphs[b].nodes.at(b).massval += localmass[1][1];
		ugraphs[b].nodes.at(c).massval += localmass[1][2];

		ugraphs[c].nodes.at(a).massval += localmass[2][0];
		ugraphs[c].nodes.at(b).massval += localmass[2][1];
		ugraphs[c].nodes.at(c).massval += localmass[2][2];
	}
}

inline void populateFval()
{
	Real temp = 0.0;
	for (size_t i =0; i < novert; i++)
	{
		for (size_t k = 0; k < ugraphs[i].nodes.size(); k++)
		{
			//cout << "## 7777 ###" << ugraphs[i].index.size() << '\n';
			size_t id = ugraphs[i].index[k];
			//cout << "## 8888 ###" << '\n';
			temp += ugraphs[i].nodes.at(id).massval * unodes[id].uval;
		}
		unodes[i].fval = temp;
		temp = 0.0;
	}
}

inline void solveGS()
{
	Real r = 0.0,norm = 0.0, temp=0.0;
	vector<double> preu;
	for (size_t i = 0; i < novert; i++)
	{
		preu.emplace_back(0.0);
	}
	
	do{
		for (size_t i = 0; i < novert; i++)
		{
			for (size_t k = 0; k < ugraphs[i].nodes.size(); k++)
			{
				
				size_t id = ugraphs[i].index[k];
				if (i == id) continue;
				temp += ugraphs[i].nodes.at(id).stiffval * unodes[id].uval;

			}
			preu[i] = unodes[i].uval;
			unodes[i].uval = (unodes[i].fval - temp) / ugraphs[i].nodes.at(i).stiffval;
			temp = 0.0;
		}
		norm = 0.0;
		r = 0.0;
		for (size_t i = 0; i < novert; i++)
		{
			r = unodes[i].uval-preu[i];
			norm += r*r;
		}
		norm = sqrt(norm/novert);
		cout << norm << '\n';
	} while (norm > eps);
}

inline void solveCG()
{
	vector<Real> res, dirc,z;
	Real temp = 0.0, del0 = 0.0, del1 = 0.0, denom = 0.0, alpha = 0.0, beta = 0.0;
	size_t id = 0;
	//cout << "555" << '\n';
	for (size_t i = 0; i < novert; i++)
	{
		for (size_t k = 0; k < ugraphs[i].nodes.size(); k++)
		{
			id = ugraphs[i].index[k];
			temp += ugraphs[i].nodes.at(id).stiffval * unodes[id].uval;
		}
		res.emplace_back(unodes[i].fval - temp);
		dirc.emplace_back(unodes[i].fval - temp);
		z.emplace_back(0.0);
		temp = 0.0;
	}

	//for (size_t i = 0; i < novert; i++)
	//{
	//	//unodes[i].uval /= normu;
	//	cout << unodes[i].fval << " ";
	//}

	//for (size_t i = 0; i < 10; i++)
	//{
	//	//unodes[i].uval /= normu;
	//	cout << dirc[i] << " ";
	//}
	
	for (size_t i = 0; i < res.size(); i++)
		del0 += res[i] * res[i];
	//dirc = res;
	while (sqrt(del0) > eps)
	{
		
		for (size_t i = 0; i < novert; i++)
		{
			for (size_t k = 0; k < ugraphs[i].nodes.size(); k++)
			{
				id = ugraphs[i].index[k];
				temp += ugraphs[i].nodes.at(id).stiffval * dirc[id];
			}
			z[i] = temp;
			temp = 0.0;
		}

		//cout << "\n res = ";
		//for (size_t i = 0; i < 10; i++)
		//{
		//	//unodes[i].uval /= normu;
		//	cout << res[i] << " ";
		//}
		//cout << '\n';

		//cout << "\n zzzz  = ";
		//for (size_t i = 0; i < 10; i++)
		//{
		//	//unodes[i].uval /= normu;
		//	cout << z[i] << " ";
		//}
		//cout << "\n dirc = ";
		//for (size_t i = 0; i < 10; i++)
		//{
		//	//unodes[i].uval /= normu;
		//	cout << dirc[i] << " ";
		//}

		denom = 0.0;
	for (size_t i = 0; i < dirc.size(); i++)
			denom += dirc[i] * z[i];

	//cout << "denom = " << denom;
	
	alpha = del0 / denom;
	del1 = 0.0;
	//cout << "  alpha = " << alpha;
	for (size_t i = 0; i < novert; i++)
	{
		unodes[i].uval += alpha * dirc[i];	
		res[i] -= alpha * z[i];
		del1 += res[i] * res[i];
	}
	
	//cout << "delta = " << del1 << " " << del0;
	
	if (sqrt(del1) <= eps)
		return;
	beta = del1 / del0;
	//cout << "beta = " << beta;
	//cout << '\n';
	for (size_t i = 0; i < novert; i++)
	{
		dirc[i] = res[i] + beta * dirc[i];
	}

	del0 = del1;			
 }
}

inline void invPower(Real& lambda)
{
	Real lambdaold = 0.0, normu = 0.0;
	vector<Real> num, denom;
	Real t1 = 0.0, t2 = 0.0;
	size_t id = 0;
	for (size_t i = 0; i < novert; i++)
	{
		num.emplace_back(0.0);
		denom.emplace_back(0.0);
	}
	do{		
		lambdaold = lambda;
		populateFval();
		solveCG();
		//solveGS();
		//cout << "555" << '\n';
		normu = 0.0;
		for (size_t i = 0; i < novert; i++)
		{
			normu += unodes[i].uval * unodes[i].uval;
		}
		
		normu = sqrt(normu);
		//cout << normu;
		for (size_t i = 0; i < novert; i++)
		{
			unodes[i].uval = unodes[i].uval / normu;
			//cout << unodes[i].uval << " ";
		} 
		//cout << '\n';
		

		for (size_t i = 0; i < novert; i++)
		{
			for (size_t k = 0; k < ugraphs[i].nodes.size(); k++)
			{
				id = ugraphs[i].index[k];
				t1 += ugraphs[i].nodes.at(id).stiffval * unodes[id].uval;
				t2 += ugraphs[i].nodes.at(id).massval * unodes[id].uval;
			}
			
			num[i] = t1;
			denom[i] = t2;
			t1 = 0.0;
			t2 = 0.0;
		}
		Real n = 0.0, d = 0.0;
		for (size_t i = 0; i < novert; i++)
		{
			n += unodes[i].uval * num[i];
			d += unodes[i].uval * denom[i];
		}

		lambda = n / d;	
		cout << '\n' << lambda;
	} while ((abs(lambda - lambdaold)/lambdaold) > ERRLIMIT);
}

inline bool compareFiles(string sfile, string tfile)
{
	string sline, tline;
	bool flag = true;
	ifstream srcfile, tgtfile;

	srcfile.open(sfile);
	tgtfile.open(tfile);

	while (getline(srcfile, sline) && getline(tgtfile, tline))
	{
		if (sline != tline)
			flag = false;
	}

	return flag;
}

inline bool compareEigenFiles(string sfile, string tfile)
{
	string sline, tline;
	bool flag = true;
	ifstream srcfile, tgtfile;
	Real a, b, c, g, e, f;

	srcfile.open(sfile);
	tgtfile.open(tfile);

	while (srcfile >> a && srcfile >> b && srcfile >> c && tgtfile >> g && tgtfile >> e && tgtfile >> f)
	{
		if (abs(a-g) > 0.00005 || abs(b- e) > 0.00005  || abs(c-f)>0.00005)
			flag = false;
	}

	return flag;
}

int main(int argc, char** argv)
{

	if (argc < 3)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}
	Real lambda = 2.0;

	delta = atof(argv[1]);
	eps = atof(argv[2]);
    	
	init();
	createGlobalMatrix();
	invPower(lambda);
	
	cout << "\n Eigenvalue = " << lambda;

	cout << "\n Writing solution to files... \n\n";

	std::string fname1 = std::string("ksq.txt");
	std::ofstream	fOut1(fname1);
	
	for (size_t i = 0; i < novert; ++i) 
	{
		fOut1 << knodes[i].xcord << " " << knodes[i].ycord << " " << knodes[i].uval << std::endl;			
	}
		fOut1 << std::endl;
	
	fOut1.close();

	std::string fname2 = std::string("A.txt");
	std::ofstream	fOut2(fname2);
	std::string fname3 = std::string("M.txt");
	std::ofstream	fOut3(fname3);
	size_t id;
	for (size_t i = 0; i < novert; ++i)
	{
		for (size_t k = 0; k < ugraphs[i].nodes.size(); k++)
		{
			id = ugraphs[i].index[k];
			fOut2 << i << " " << id << " " << ugraphs[i].nodes[id].stiffval << std::endl;
			fOut3 << i << " " << id << " " << ugraphs[i].nodes[id].massval << std::endl;
		}
		
	}
	fOut2 << std::endl;
	fOut3 << std::endl;

	fOut2.close();
	fOut3.close();

	std::string fname4 = std::string("eigenmode.txt");
	std::ofstream	fOut4(fname4);

	for (size_t i = 0; i < novert; ++i)
	{
		fOut4 << unodes[i].xcord << " " << unodes[i].ycord << " " << unodes[i].uval << std::endl;
	}
//	fOut4 << std::endl;

	fOut4.close();
	

	string srcfile,tgtfile;
	srcfile = "./reference-outputs/ksq-ref.txt";
	tgtfile = "./ksq.txt";
	
	if (compareFiles(srcfile, tgtfile))
		cout << "ksq.txt is correct with reference file \n\n";
	else
		cout << "ksq.txt is not correct with reference file \n\n";

	srcfile = "./reference-outputs/A-ref.txt";
	tgtfile = "./A.txt";

	if (compareFiles(srcfile, tgtfile))
		cout << "A.txt is correct with reference file \n\n";
	else
		cout << "A.txt is not correct with reference file \n\n";

	srcfile = "./reference-outputs/M-ref.txt";
	tgtfile = "./M.txt";

	if (compareFiles(srcfile, tgtfile))
		cout << "M.txt is correct with reference file \n\n";
	else
		cout << "M.txt is not correct with reference file \n\n";

	srcfile = "./reference-outputs/eigenmode-ref.txt";
	tgtfile = "./eigenmode.txt";

	if (compareEigenFiles(srcfile, tgtfile))
		cout << "eigenmode.txt is correct with reference file \n\n";
	else
		cout << "eigenmode.txt is not correct with reference file \n\n";
	
	/*free(ugraphs);
	free(tri);
	free(unodes);
	free(knodes);*/
	
	return 0;

}
