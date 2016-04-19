
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <string.h>
#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */
#include <math.h>  
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>

using namespace std;

class Edge {
 public:
  string n1;
  string n2;
  double length;

  Edge( string n1 , string n2 , double length) { 
    init( n1 , n2 , length);
  }

  void init (  string n1 , string n2 , double length) {
   this->n1 = n1 ;
   this->n2 = n2 ; 
   this->length = length;
  }

  Edge ( string line ) {
   char* cstr = new char [line.size()+1];
   strcpy (cstr, line.c_str());
   char* p1=strtok (cstr,",");
   char* p2=strtok (NULL,",");
   char* p3=strtok (NULL,",");
   double len = atof (p3);
   init( p1, p2, len);
 }

 void print() {
  cout << " { " << n1 << " , " << n2 << " , " << length << " } " << endl;
 }
   
};

struct CompareEdge {
public:
    bool operator()(Edge t1, Edge t2)
    {
       if (t1.length < t2.length) return true;
       return false;
    }
} myobj;

void addEdges( vector<Edge> edges, priority_queue <Edge , vector<Edge> , CompareEdge >&   sortedEdges ) {
  for ( int i=0;i<edges.size();i++) {
	sortedEdges.push(edges[i]);
  }
}

string getParentSet(  map<string, string >& nodes , string node) {

  string child = node;
  string parent = nodes[child];
  int levelsTraversed = 0;
  while( child.compare(parent) != 0 ) {
    string tmp = parent ;
    string prevChild = child;
    child = parent;
    parent = nodes[tmp];  
    nodes[prevChild] = parent; // link to parent directly
    levelsTraversed++;

  } 
  if ( levelsTraversed > 5 ) 
    cout << " levels traversed is " << levelsTraversed << endl;
  return parent;
}

int kruskal( map<string, string > nodes , vector<Edge> edges ) {
  if( nodes.empty() ) 
     return -1;
  
  sort(edges.begin() , edges.end() , myobj);

  double totalCost =0.0;
  //initialize e_new
  vector<Edge> e_new ;

  int totalEdgesConsidered=0;
  
  for( int i=0;i<edges.size() ; i++ ) {
     Edge e1 = edges[i];
     totalEdgesConsidered++;
     string n1_parent = getParentSet( nodes, e1.n1) ;
     string n2_parent = getParentSet( nodes, e1.n2) ;
     int smaller = n1_parent.compare(n2_parent);
     if ( smaller == 0 ) {
       //same set ignore 
     } else if ( smaller < 0 ) {
	 e1.print();
	 e_new.push_back(e1);
	 totalCost += e1.length;
         nodes[n1_parent] = n1_parent ;
         nodes[n2_parent] = n1_parent ;
     } else if ( smaller > 0 ) {
	 e1.print();
	 e_new.push_back(e1);
	 totalCost += e1.length;
         nodes[n1_parent] = n2_parent ;
         nodes[n2_parent] = n2_parent ;
     } else {
       cout << " unexpected string comparison result " << endl ;
     } 

     if( e_new.size() %10 == 0 ) {
         cout << " Current MS edge set size : " << e_new.size() << endl;
     }
     if( totalEdgesConsidered % 100 == 0 ) {
         cout << " Total edges considered so far :" << totalEdgesConsidered << endl;
     }
  }

  cout << " # edges in e_new is " << e_new.size() << "and cost is " << totalCost << endl;
   set<string> connectedComponents;
   for( map<string , string >::iterator iter = nodes.begin() ; iter!= nodes.end(); iter++) {
	string parent = getParentSet(nodes , iter->first);
	connectedComponents.insert(parent);
  }
  cout << " # connected components is " << connectedComponents.size() << endl ;

  return 0;

}

void readFile(string str) {
  std::ifstream infile(str.c_str() );
  std::string line;
  vector<Edge> edges;
  while (std::getline(infile, line))
  {
    Edge e1(line);
    edges.push_back(e1);
  }
  map<string , string > nodes;
  cout << "Finished reading edges " <<endl;
  for( int i=0;i<edges.size();i++) {
      nodes[edges[i].n1] = edges[i].n1 ;
      nodes[edges[i].n2] = edges[i].n2 ;
  }

  cout << "Read " << edges.size() << " edges  and found " << nodes.size() << " nodes " << endl; 
  cout << " Running kruskal now " <<endl ;
  kruskal(nodes,edges);
}

int main() {

cout << "hello world " ;
//test1();

//readFile("tree2.csv");
readFile("tree_x.txt");
}
