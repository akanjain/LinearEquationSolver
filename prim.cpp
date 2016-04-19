
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

class CompareEdge {
public:
    bool operator()(Edge& t1, Edge& t2)
    {
       if (t1.length > t2.length) return true;
       return false;
    }
};

void addEdges( vector<Edge> edges, priority_queue <Edge , vector<Edge> , CompareEdge >&   sortedEdges ) {
  for ( int i=0;i<edges.size();i++) {
	sortedEdges.push(edges[i]);
  }
}


int prim( map<string, vector<Edge> > nodes , vector<Edge> edges ) {
  if( nodes.empty() ) 
     return -1;
  
  // initialize v_new
  std::set <string> v_new ;
  v_new.insert( nodes.begin()->first);
  cout << "First node is " <<  nodes.begin()->first ;

  priority_queue <Edge , vector<Edge> , CompareEdge >  cEdges ;
  cout << (nodes.begin()->second).size() <<endl; 
  addEdges( nodes.begin()->second ,  cEdges );

  double totalCost =0.0;
  //initialize e_new
  vector<Edge> e_new ;

  int totalEdgesConsidered=0;
  
  while ( v_new.size() != nodes.size() ) { 
     if ( cEdges.size() == 0 ) {
        cout << " No Candidate Edge...graph is not connected..total edges considered " << totalEdgesConsidered  <<endl ;
        cout << " # nodes in v_new is " << v_new.size() << " and # edges in e_new is " << e_new.size() << "and cost is " << totalCost << endl;
        return 0;
     }

     Edge e1 =  cEdges.top(); 
     cEdges.pop();
     totalEdgesConsidered++;
     bool v_new_contains_e1_n1 = ( v_new.find(e1.n1) != v_new.end());
     bool v_new_contains_e1_n2 = ( v_new.find(e1.n2) != v_new.end());
     if ( v_new_contains_e1_n1 && ! v_new_contains_e1_n2 ) {
	   e1.print();
           v_new.insert(e1.n2);
           e_new.push_back(e1);
           totalCost += e1.length;
           addEdges(nodes[e1.n2] , cEdges);
     } else if ( !v_new_contains_e1_n1 && v_new_contains_e1_n2 ) {
	   e1.print();
           v_new.insert(e1.n1);
	   e_new.push_back(e1);
	   totalCost += e1.length;
           addEdges(nodes[e1.n1] , cEdges);
     } else if ( ! v_new_contains_e1_n1 && ! v_new_contains_e1_n2 ) {
        cout << " unexpected edge found in cEdges " << endl;
        e1.print();
     }

     if( v_new.size() %1000 == 0 ) {
         cout << " Current MST size : " << v_new.size() << endl;
     }
  }

  cout << " # nodes in v_new is " << v_new.size() << " and # edges in e_new is " << e_new.size() << "and cost is " << totalCost << endl;
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
  map<string , vector<Edge> > nodes;
  cout << "Finished reading edges " <<endl;
  for( int i=0;i<edges.size();i++) {
    if ( nodes.find( edges[i].n1) == nodes.end() ) {
      vector<Edge> edge_list;
      edge_list.push_back(edges[i]);
      nodes[edges[i].n1] = edge_list;
   } else {
     vector<Edge> edge_list =  nodes[edges[i].n1] ;
     edge_list.push_back(edges[i]);
     nodes[edges[i].n1] = edge_list;
   }
 
    if ( nodes.find( edges[i].n2) == nodes.end() ) {
      vector<Edge> edge_list;
      edge_list.push_back(edges[i]);
      nodes[edges[i].n2] = edge_list;
   } else {
     vector<Edge> edge_list =  nodes[edges[i].n2] ;
     edge_list.push_back(edges[i]);
     nodes[edges[i].n2] = edge_list;
   }
  }

  cout << "Read " << edges.size() << " edges  and found " << nodes.size() << " nodes " << endl; 
  //maybe i want to print edge count for individual node to do some spot checks
  for( map<string , vector<Edge> >::iterator iter = nodes.begin() ; iter!= nodes.end(); iter++) {
    //cout << iter->first << iter->second.size() << endl;
  }
  cout << " Running prim now " <<endl ;
  prim(nodes,edges);
}

int main() {

cout << "hello world " ;
//test1();

//readFile("tree2.csv");
readFile("tree_y.txt");
}
