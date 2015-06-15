#include <iostream>
#include <map>
#include<memory>
#include<malloc.h>
#include<stdlib.h>
using namespace std;

struct node{
std::map<char,int> mymap;
};
int main ()
{
 // std::map<char,int> mymap;
    node* n;
  // first insert function version (single parameter):
  n = (node*) memalign(64,4*sizeof(node));
  n[0].mymap.insert ( std::pair<char,int>('a',100) );
  n[0].mymap.insert ( std::pair<char,int>('z',200) );
  std::cout << "aaa" << n[0].mymap.at('z') << '\n'; 

/*
  
mymap.insert (it, std::pair<char,int>('b',300));  // max efficiency inserting
  mymap.insert (it, std::pair<char,int>('c',400));  // no max efficiency inserting

  // third insert function version (range insertion):
  std::map<char,int> anothermap;
  anothermap.insert(mymap.begin(),mymap.find('c'));

  // showing contents:
  std::cout << "mymap contains:\n";
  for (it=mymap.begin(); it!=mymap.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';

  std::cout << "anothermap contains:\n";
  for (it=anothermap.begin(); it!=anothermap.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';
*/
  return 0;
}
