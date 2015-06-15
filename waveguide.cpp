#include<iostream>
#include<fstream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#include<algorithm>
#include <stdio.h>
#include <string.h>
#include<vector>
#include <map>
#include "./myColsamm/Source/Colsamm.h"
#include <iomanip>
#define ERRLIMIT 0.000000001

typedef double Real;
typedef int unint;
#define ALLIGNMENT 64
  
using namespace std;
using namespace ::_COLSAMM_;

struct node{
	Real xcord, ycord, stiffval, massval, uval, fval;
	size_t vertno;
	};

struct graph{
	 
     std::map<size_t,node> nodes;
	 std::vector<int> index;
     size_t length;
};

struct triang{
	size_t vertex[3];
};

Real delta, eps, reflevel;

graph* __restrict ugraphs = nullptr;
node* __restrict knodes = nullptr;
node* __restrict unodes = nullptr;
triang * __restrict tri = nullptr , * __restrict otri = nullptr;
size_t novert, notriangle, tricnt = 0,b4ref;

inline double kxy2(const Real x, const Real y)
{
	return (((100.0 + delta) * exp(-50.0*(x*x + y*y))) - 100.0);
}

inline size_t searchNode(node nd)
{

    for(size_t k = b4ref ; k < novert; k++)
    {
        if(abs(unodes[k].xcord - nd.xcord) < 0.000005 && abs(unodes[k].ycord - nd.ycord) < 0.000005)
            return k;
    }

    ++novert;
  /*  node* newnode = (node *)memalign(ALLIGNMENT, novert * sizeof(node));
    memcpy(newnode,unodes ,(novert-1)* sizeof(node));

    newnode[novert-1].xcord = nd.xcord;
    newnode[novert-1].ycord = nd.ycord;
    newnode[novert-1].fval = 0.0;
    newnode[novert-1].uval = 1.0;
    newnode[novert-1].massval = 0.0;
    newnode[novert-1].stiffval = 0.0;
    newnode[novert-1].vertno = novert;
    free(unodes);
    unodes = newnode;*/
    unodes = (node *)realloc( unodes, novert * sizeof(node));
  //  memcpy(newnode,unodes ,(novert-1)* sizeof(node));

    unodes[novert-1].xcord = nd.xcord;
    unodes[novert-1].ycord = nd.ycord;
    unodes[novert-1].fval = 0.0;
    unodes[novert-1].uval = 1.0;
    unodes[novert-1].massval = 0.0;
    unodes[novert-1].stiffval = 0.0;
    unodes[novert-1].vertno = novert;
    //free(unodes);
    //unodes = newnode;

    return (novert-1);
}

inline void refinement(triang trin, size_t rfl)
{
      triang tmptri;
      size_t a,b,c;

        node newVert1 , newVert2, newVert3;
        size_t vert1 = trin.vertex[0];
        size_t vert2 = trin.vertex[1];
        size_t vert3 = trin.vertex[2];

         newVert1.xcord = (unodes[vert1].xcord + unodes[vert2].xcord) * 0.5;
	 newVert1.ycord = (unodes[vert1].ycord + unodes[vert2].ycord) * 0.5;
	 newVert2.xcord = (unodes[vert2].xcord + unodes[vert3].xcord) * 0.5;
	 newVert2.ycord = (unodes[vert2].ycord + unodes[vert3].ycord) * 0.5;
	 newVert3.xcord = (unodes[vert1].xcord + unodes[vert3].xcord) * 0.5;
	 newVert3.ycord = (unodes[vert1].ycord + unodes[vert3].ycord) * 0.5;

         a = searchNode(newVert1);
         b = searchNode(newVert2);
         c = searchNode(newVert3);

         tmptri.vertex[0] = a;
         tmptri.vertex[1] = b;
         tmptri.vertex[2] = c;
         
	 if(rfl > 1)
             refinement(tmptri,rfl - 1);
         else{
             tri[tricnt++] = tmptri;
         }

         tmptri.vertex[0] = vert1;
         tmptri.vertex[1] = a;
         tmptri.vertex[2] = b;
         if(rfl > 1)
             refinement(tmptri,rfl - 1);
         else{
             tri[tricnt++] = tmptri;
         }

	 tmptri.vertex[0] = vert2;
         tmptri.vertex[1] = a;
         tmptri.vertex[2] = c;
         if(rfl > 1)
             refinement(tmptri,rfl - 1);
         else{
             tri[tricnt++] = tmptri;
         }

         tmptri.vertex[0] = vert3;
         tmptri.vertex[1] = b;
         tmptri.vertex[2] = c;
         if(rfl > 1)
             refinement(tmptri,rfl - 1);
         else{
             tri[tricnt++] = tmptri;
         }
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
  	b4ref = novert;
        unodes = (node *)memalign(ALLIGNMENT,novert*sizeof(node));
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
	}

	getline(ucircle, tmp);	

	getline(ucircle, tmp);
	
	notriangle = stoi(tmp.substr(0, tmp.find(" ") - 1));;
    otri = new triang[notriangle];
   size_t ntr = notriangle;
    notriangle = pow(2,2*reflevel) * notriangle;
	//ugraphs = new graph[novert];
	tri = new triang[notriangle];
	getline(ucircle, tmp);

    for (size_t i = 0; ucircle >> d && ucircle >> e && ucircle >> f; i++)
    {
        otri[i].vertex[0] = d;
        otri[i].vertex[1] = e;
        otri[i].vertex[2] = f;
    }

    if(reflevel > 0)
    {
     for(size_t i = 0;i<ntr; i++)
      refinement(otri[i],reflevel);
	delete otri;
    }
    else
        tri = otri;
    cout << "no of vertex = " << novert << '\n';
    cout << "no of triangle = " << notriangle << '\n';
    ugraphs = new graph[novert];
    for (size_t i = 0;i<notriangle; i++)
	{
        d = tri[i].vertex[0];
        e = tri[i].vertex[1];
        f = tri[i].vertex[2];

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
	}

	knodes = new node[novert];
	for (size_t i=0; i < novert; i++)
	{
		knodes[i].xcord = unodes[i].xcord;
		knodes[i].ycord = unodes[i].ycord;
		knodes[i].uval = kxy2(unodes[i].xcord, unodes[i].ycord);
		std::sort(ugraphs[i].index.begin(), ugraphs[i].index.end());
	}
}

inline void createLocalMatrix(size_t a, size_t b, size_t c, std::vector<std::vector<double>>& localstiff, std::vector<std::vector<double>>& localmass)
{
	ELEMENTS::Triangle my_element;

		std::vector<double> corners(6, 0.0);
	corners[0] = ugraphs[a].nodes.at(a).xcord; 
	corners[1] = ugraphs[a].nodes.at(a).ycord;
	corners[2] = ugraphs[b].nodes.at(b).xcord;
	corners[3] = ugraphs[b].nodes.at(b).ycord;
	corners[4] = ugraphs[c].nodes.at(c).xcord;
	corners[5] = ugraphs[c].nodes.at(c).ycord;

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
		
		createLocalMatrix(a, b, c, localstiff, localmass);
		
		ugraphs[a].nodes.at(a).stiffval += localstiff[0][0];
		ugraphs[a].nodes.at(b).stiffval += localstiff[0][1];
		ugraphs[a].nodes.at(c).stiffval += localstiff[0][2];
		
		ugraphs[b].nodes.at(a).stiffval += localstiff[1][0];
		ugraphs[b].nodes.at(b).stiffval += localstiff[1][1];
		ugraphs[b].nodes.at(c).stiffval += localstiff[1][2];
		
		ugraphs[c].nodes.at(a).stiffval += localstiff[2][0];
		ugraphs[c].nodes.at(b).stiffval += localstiff[2][1];
		ugraphs[c].nodes.at(c).stiffval += localstiff[2][2];
		
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
				size_t id = ugraphs[i].index[k];
				temp += ugraphs[i].nodes.at(id).massval * unodes[id].uval;
		}
		unodes[i].fval = temp;
		temp = 0.0;
	}
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

		
		denom = 0.0;
	for (size_t i = 0; i < dirc.size(); i++)
			denom += dirc[i] * z[i];

	alpha = del0 / denom;
	del1 = 0.0;
	
	for (size_t i = 0; i < novert; i++)
	{
		unodes[i].uval += alpha * dirc[i];	
		res[i] -= alpha * z[i];
		del1 += res[i] * res[i];
	}
	
	if (sqrt(del1) <= eps)
		return;
	beta = del1 / del0;
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
	size_t id = 0 ,k = 0;
	for (size_t i = 0; i < novert; i++)
	{
		num.emplace_back(0.0);
		denom.emplace_back(0.0);
	}
	do{	
		k++;
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
		cout << setprecision(10) << "\nEigenvalue after step: " << k << " = " << lambda;
		
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

	if (argc < 4)
	{
		std::cout << "Invalid number of argument";
		exit(0);
	}
	Real lambda = 2.0;

	delta = atof(argv[1]);
	eps = atof(argv[2]);
    reflevel = atoi(argv[3]);
    	
	init();
	createGlobalMatrix();
	invPower(lambda);
	
	cout << "\n Smallest Eigenvalue = " << lambda;

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
	fOut4 << std::endl;

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
