#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<map>
#include<list>
#include <algorithm>
#include <vector>
#include <math.h>
#include<sys/time.h>
using namespace std;

class parse
{
	public:
	void parseinputnodes (const char *);
	void parseinputpl (const char *);
	void parsematrix (const char *);
	void initsolution_x ();
	void initsolution_y ();
	void HPWL_x();
	void HPWL_y();
};	

double cpuTime(void)
{
    		struct timeval cputime;
    		gettimeofday(&cputime, NULL);
     		return static_cast<double>(cputime.tv_sec) + (1.e-6)*static_cast<double>(cputime.tv_usec);
}

struct node {
string name;
unsigned id;
float x, y;
bool fixed;
map <string, float> adjmap_x;
map <string, float> adjmap_y;
};


struct edge {
string source, target;
};

struct matrix_x {
unsigned row;
float data;
list<pair<int,float> > movelist;
list<pair<int,float> > pinlist;
};

struct matrix_y {
unsigned row;
float data;
list<pair<int,float> > movelist;
list<pair<int,float> > pinlist;
};


struct param {
float r;
float p;
};

struct net {
string net;
map<string,float> lent_x;
map<string,float> lent_y;
};

static map <int, float> sol_map;
static map <int, float> temp_map;
static map <int, param*> par_map;
static map <string, node*> name_to_node_map;
static map <int, node*> id_to_node_map;
static map <int, matrix_x*> matrix_mapx;
static map <int, matrix_y*> matrix_mapy;
static list<pair<edge*, float> >alist;
static list<pair<edge*, float> >fixlist;
//static vector <pair<int, matrix_x*> > matrix_vecx;
//static vector <pair<int, matrix_y*> > matrix_vecy;
//static vector <pair<int, float> > sol_vec;
//static vector <pair<int, float> > temp_vec;
//static vector <pair<int, param*> > par_vec;
static map <int, net*> hpl;
static vector<pair<int, float> >res_vec;

unsigned num_nodes = 0;
int countmove = 0;
unsigned num = 0;	

void parse :: parseinputpl (const char *filename)
{
	int count1 = 0;
	string name, line, prop;
	float w, h;
	int g = 0;
	int k = 0;
	int number = 0;
	int min, max = 0;
        ifstream myfile (filename);
	ofstream outf ("sort.txt");
        if (myfile.is_open())
        {
                while (myfile.good())
                {
			
			while (getline(myfile,line))
			{		
				if (myfile.peek() == '#' || myfile.peek() == '\n')
                        	{
                                getline (myfile,line);
                        	}
				if (line.find("")!= string :: npos)
                                {
                                }

				if (line.find("UCLA")!= string :: npos)
				{
				}
				else
				{
				istringstream ss(line);
				ss >> name >> w >> h >> prop >> prop >> prop;
//                                cout << name << "\t" << w << "\t" << h << "\t" << prop << endl;
				number = num_nodes;
				struct node* n = new node;
				n->name = name;
				n->id = num_nodes++;
				n->x = w;
				n->y = h;
		              	if(prop == "/FIXED" || prop == "/FIXED_NI") 
				{
					n->fixed = true;
					count1 = count1++;
				}
				else
				{
					n->fixed = false;
					countmove = countmove++;
				}
				name_to_node_map[name] = n;
				//name_vec.push_back(pair <string, node*>  (name, n));
				//id_vec.push_back(pair <int, node*>  (number, n));
				id_to_node_map[number] = n;
				}	
			}
		
		}	
		myfile.close();
/* 
	for (map<string,node*>::iterator it=name_to_node_map.begin(); it!= name_to_node_map.end(); ++it)
	{
		cout << it->first << "\t" << it->second->x << "\t" << it->second->y << "\t" << it->second->fixed << "\t" << it->second->id << endl;
		if (it->second->fixed == 1)
		{	
        	k = k + 1;
		if (k == 1)
		{
			min = it->second->x;
			max = it->second->x;
		}
		else
		{
	        	if ((it->second->x) < min)
			{
				min = it->second->x;
			}
			if ((it->second->x) > max)
			{
				max = it->second->x;
			}
        	}
		}
	}
	cout << "number of terminal nodes: "<< count1 << endl;
	cout << "number of movable nodes: " << countmove << endl;
	cout << "min = " << min << endl;
	cout << "max = " << max << endl;

	for (vector<pair<string, node*> >:: iterator it = name_vec.begin(); it!= name_vec.end(); it++)
	{
		cout << "name: " << it->first << "name : " << it->second->name << "x = " << it->second->x << "y = " << it->second->y << " status = " << it->second->fixed << " id = " << it->second->id << endl; 
	}
*/
	}

	else
	cout << "Unable to open file";
}
	

void parse :: parsematrix (const char *filename)
{
	int count1 = 0;
	int id1;
        string name, line, prop, maxnode_x, minnode_x, maxnode_y, minnode_y, netname;
        float min_x, max_x, min_y, max_y = 0;
	int degree = 0;
        int g = 0;
        int k = 0;
	float e = 1;
	float xc, yc = 0;
	float wt, wt1, wt2, a;
	int num = 0;
	int iter = 0;
        ifstream myfile (filename);
	ofstream outtree_x ("tree_x.txt");
	ofstream outtree_y ("tree_y.txt");


	if (myfile.is_open())
	{
		while (myfile.good())
		{
		
			while (myfile.peek() == '#' || myfile.peek() == '\n')
			{
				getline (myfile,line);
			}

			getline (myfile,line);
			
			if (line.find("UCLA") != string :: npos)
			{
			}
			if (line.find("NetDegree") != string :: npos)
			{
				iter = iter++ ;
				max_x = 0;
				max_y = 0;
				istringstream ss(line);
				ss >> name >> prop >> degree;
				struct net* n = new net;
				n->net = netname;
				hpl[iter] = n;

				string node_name[degree];
				for (k = 0; k < degree; k++)
				{
					node_name[k] = "zero";
				}
				for (g = 0; g < degree; g++)
				{
					getline (myfile,line);
					istringstream ss(line);
					ss >> name;
					node_name[g] = name;

					if (g == 0)
					{
						min_x = name_to_node_map.find(name)->second->x;
						minnode_x = name;
						max_x = name_to_node_map.find(name)->second->x;
                                                maxnode_x = name;
                                                min_y = name_to_node_map.find(name)->second->y;
                                                minnode_y = name;
                                                max_y = name_to_node_map.find(name)->second->y;
                                                maxnode_y = name;

	
					}
					else
					{
						if ((name_to_node_map.find(name)->second->x) < min_x)
						{
							min_x =  name_to_node_map.find(name)->second->x;
							minnode_x = name;
						}

						if ((name_to_node_map.find(name)->second->x) > max_x)
						{
							max_x = name_to_node_map.find(name)->second->x;
							maxnode_x = name;
						}
								
						if ((name_to_node_map.find(name)->second->y) < min_y)
						{
							min_y =  name_to_node_map.find(name)->second->y;
							minnode_y = name;
						}

						if ((name_to_node_map.find(name)->second->y) > max_y)
						{
							max_y = name_to_node_map.find(name)->second->y;
							maxnode_y = name;
						}
					}
					if (max_x == min_x && g == degree-1)
					{
						maxnode_x = name;
					}
					if (max_y == min_y && g == degree-1)
                                        {
                                                maxnode_y = name;
                                        }
				
	
				}
			
		                for (int i = 0; i<degree; i++)
		                {
      
					xc = name_to_node_map.find(node_name[i])->second->x;
                                        yc = name_to_node_map.find(node_name[i])->second->y;

					hpl[iter]->lent_x[node_name[i]] = xc;
					hpl[iter]->lent_y[node_name[i]] = yc;
					
					if (node_name[i] == minnode_x || node_name[i] == maxnode_x)
					{
						for (int j=0; j<i; j++)
						{
							wt = float (1)/((degree-1)*(fabs(xc-(name_to_node_map.find(node_name[j])->second->x))+e));
							name_to_node_map[node_name[i]]->adjmap_x[node_name[j]] =  name_to_node_map[node_name[i]]->adjmap_x[node_name[j]] + wt;

						}

						for (int j= i+1; j<degree; j++)
						{
							wt = float (1)/((degree-1)*(fabs(xc-(name_to_node_map.find(node_name[j])->second->x))+e));
						 	name_to_node_map[node_name[i]]->adjmap_x[node_name[j]] =  name_to_node_map[node_name[i]]->adjmap_x[node_name[j]] + wt;


						}

					}

					else
					{
						wt = float (1)/((degree-1)*((fabs(xc-min_x))+e));
						name_to_node_map[node_name[i]]->adjmap_x[minnode_x] =  name_to_node_map[node_name[i]]->adjmap_x[minnode_x] + wt;

						wt = float (1)/((degree-1)*(fabs(xc-max_x)+e));
						name_to_node_map[node_name[i]]->adjmap_x[maxnode_x] =  name_to_node_map[node_name[i]]->adjmap_x[maxnode_x] + wt;

					}

					if (node_name[i] == minnode_y || node_name[i] == maxnode_y)
                                        {
                                                for (int j=0; j<i; j++)
                                                {
                                                        wt = float (1)/((degree-1)*(fabs(yc-(name_to_node_map.find(node_name[j])->second->y))+e));
                                                        name_to_node_map[node_name[i]]->adjmap_y[node_name[j]] =  name_to_node_map[node_name[i]]->adjmap_y[node_name[j]] + wt;

                                                }

                                                for (int j= i+1; j<degree; j++)
                                                {
                                                        wt = float (1)/((degree-1)*(fabs(yc-(name_to_node_map.find(node_name[j])->second->y))+e));
                                                        name_to_node_map[node_name[i]]->adjmap_y[node_name[j]] =  name_to_node_map[node_name[i]]->adjmap_y[node_name[j]] + wt;


                                                }

                                        }

                                        else
                                        {
                                                wt = float (1)/((degree-1)*((fabs(yc-min_y))+e));
                                                name_to_node_map[node_name[i]]->adjmap_y[minnode_y] =  name_to_node_map[node_name[i]]->adjmap_y[minnode_y] + wt;

                                                wt = float (1)/((degree-1)*(fabs(yc-max_y)+e));
                                                name_to_node_map[node_name[i]]->adjmap_y[maxnode_y] =  name_to_node_map[node_name[i]]->adjmap_y[maxnode_y] + wt;

                                        }

	

				}	

			}
		}
	}
			


     
	for (map<string,node*>::iterator it=name_to_node_map.begin(); it!= name_to_node_map.end(); ++it)
        {
		for (map<string, float> :: iterator it1 = it->second->adjmap_x.begin(); it1!= it->second->adjmap_x.end(); ++it1)
		{
			outtree_x << it->first << "," << it1->first << "," << it1->second << endl;
		}
	
		for (map<string, float> :: iterator it1 = it->second->adjmap_y.begin(); it1!= it->second->adjmap_y.end(); ++it1)
                {
                        outtree_y << it->first << "," << it1->first << "," << it1->second << endl;
                }



	}
	
	for (map<string,node*>::iterator it=name_to_node_map.begin(); it!= name_to_node_map.end(); ++it)
        {
		a = 0;
		float xco = 0;
		if (it->second->fixed == 0)
		{
		struct matrix_x* m = new matrix_x;
                id1 = it->second->id;
                m->row = it->second->id;

		for (map<string, float> :: iterator it1 = it->second->adjmap_x.begin(); it1!= it->second->adjmap_x.end(); ++it1)
		{
			if ((it->second->fixed == 0) && (name_to_node_map[it1->first]->fixed==0))
			{
				struct edge* p = new edge;
				p->source = it->first ;
				p->target = it1->first;
		//		alist.push_back(pair<edge*, float> (p,it1->second));
				m->movelist.push_back(pair<int,float> ((name_to_node_map[it1->first]->id),-(it1->second)));
		
			}
			if ((it->second->fixed == 0) && (name_to_node_map[it1->first]->fixed==1))
                        {
                                struct edge* k = new edge;
                                k->source = it->first ;
                                k->target = it1->first;
                                //fixlist.push_back(pair<edge*, float> (k,it1->second));
                		m->pinlist.push_back(pair<int,float> ((name_to_node_map[it1->first]->id),it1->second));
				xco = (xco + ((it1->second) * (name_to_node_map[it1->first]->x))); 
		        }
			a = a + it1->second;
			m->data = xco;
		}
		m->movelist.push_back(pair<int,float> ((name_to_node_map[it->first]->id),a));
		
		//matrix_vecx.push_back(pair<int, matrix_x*> (id1, m));
		matrix_mapx[id1] = m;
		}
	}

	
	for (map<string,node*>::iterator it=name_to_node_map.begin(); it!= name_to_node_map.end(); ++it)
        {
		a = 0;
		float yco = 0;
		if (it->second->fixed == 0)
		{
		struct matrix_y* m = new matrix_y;
                id1 = it->second->id;
                m->row = it->second->id;

		for (map<string, float> :: iterator it1 = it->second->adjmap_y.begin(); it1!= it->second->adjmap_y.end(); ++it1)
		{
			if ((it->second->fixed == 0) && (name_to_node_map[it1->first]->fixed==0))
			{
				struct edge* p = new edge;
				p->source = it->first ;
				p->target = it1->first;
		//		alist.push_back(pair<edge*, float> (p,it1->second));
				m->movelist.push_back(pair<int,float> ((name_to_node_map[it1->first]->id),-(it1->second)));
		
			}
			if ((it->second->fixed == 0) && (name_to_node_map[it1->first]->fixed==1))
                        {
                                struct edge* k = new edge;
                                k->source = it->first ;
                                k->target = it1->first;
                                //fixlist.push_back(pair<edge*, float> (k,it1->second));
                		m->pinlist.push_back(pair<int,float> ((name_to_node_map[it1->first]->id),it1->second));
				yco = (yco + ((it1->second) * (name_to_node_map[it1->first]->y))); 
		        }
			a = a + it1->second;
			m->data = yco;
		}
		m->movelist.push_back(pair<int,float> ((name_to_node_map[it->first]->id),a));
		
		//matrix_vecy.push_back(pair<int, matrix_y*> (id1, m));
		matrix_mapy[id1] = m;
		}
	}

}

void parse :: initsolution_x ()
{
int num1, num2 = 0;
int id2 = 0;
float ax, ap = 0;
float r, p;
float numerator, denominator, numerator1, alpha, beta = 0;
float xnew, rnew, pnew;
float xdiff, xdiff1;
ofstream myfile ("initsolution.txt");
//ofstream myfile1 ("ap.txt");
	for (int i = 0; i < countmove; i++)
        {
        	sol_map[i] = 0;
	}		

	
	for (map<int,matrix_x*> ::iterator it = matrix_mapx.begin(); it!= matrix_mapx.end(); it++)
	{
                ax = 0;
		for (list<pair<int, float> > :: iterator it1 = it->second->movelist.begin(); it1!= it->second->movelist.end(); ++it1)
		{
			ax = float (ax + ((it1->second)*(sol_map[it1->first])));
		}
		r = ((it->second->data) - ax);
		p = r;	
		struct param* g = new param;
                g->r = r;
		g->p = p;
		par_map[it->first] = g;
		//par_map[it->first] = g;
	}

float rmax = 0;
int iteration = 0;
int check = 4;
int iter = 0;
float rtemp = 0;
float res = 0;
do
{
	++iteration;
	rmax = 1;
	check = 0;
	for (map<int,matrix_x*>::iterator it = matrix_mapx.begin(); it!= matrix_mapx.end(); it++)
	{
		check = check + 1;
		ap = 0;
		for (list<pair<int, float> > :: iterator it1 = it->second->movelist.begin(); it1!= it->second->movelist.end(); ++it1)
                {
                        ap = float (ap + ((it1->second)*(par_map[it1->first]->p)));
		}
		temp_map[it->first] = ap;
		//temp_map[it->first] = ap;
			//if (it->first == it1->first)
			//{
			//cout << "par_map = " << par_map[it1->first]->p << endl;
                	//}
	}

	numerator = 0;
	denominator = 0;
	for (map<int,param*> ::iterator it = par_map.begin(); it!= par_map.end(); it++)
	{
		numerator = numerator + ((it->second->r)*(it->second->r));
		denominator = denominator + ((it->second->p)*(temp_map[it->first]));
	}
	alpha = float (numerator/denominator);
	
       
	for (map<int,float>::iterator it=sol_map.begin(); it!= sol_map.end(); ++it)
        {
		xnew = 0;
                xnew = float ((it->second) + (alpha*(par_map[it->first]->p)));
        	sol_map[it->first] = xnew;
	}
	
	for (map<int,param*> ::iterator it=par_map.begin(); it!= par_map.end(); ++it)
        {
		rnew = 0;
		rnew = float ((it->second->r) - (alpha*(temp_map[it->first])));
		par_map[it->first]->r = rnew;
		res = float (rnew/((matrix_mapx[it->first]->data)));
		res_vec.push_back(pair <int, float> (it->first, res));
                if (fabs(res) < rmax)
                {
                     rmax = fabs(res);
                }

	}

	numerator1 = 0;
	for (map<int,param*> ::iterator it = par_map.begin(); it!= par_map.end(); it++)
	{
		numerator1 = numerator1 + ((it->second->r)*(it->second->r));
	}

	beta = float (numerator1/numerator);


	for (map<int,param*> ::iterator it=par_map.begin(); it!= par_map.end(); ++it)
        {
		pnew = 0;
                pnew = float ((it->second->r) + (beta*(it->second->p)));
        	par_map[it->first]->p = pnew;
	}
	cout << " rmax is : " << rmax << endl;

	temp_map.clear();
	res_vec.clear();

}while ((iteration < 50)|| (rmax > 1e-5));

	cout << "iteration = " << iteration << endl;


        for (map<int,float> ::iterator it=sol_map.begin(); it!= sol_map.end(); ++it)
        {
		myfile << "Node is: " << id_to_node_map[it->first]->name << "\t";
		myfile << "x = " << it->second << endl;
	}


	for (map<string,node*>::iterator it=name_to_node_map.begin(); it!= name_to_node_map.end(); ++it)
	{
		if (it->second->fixed == 0)
		{
		it->second->x = sol_map[it->second->id];
		}
		it->second->adjmap_x.clear();
	}

	for (map<int,node*>::iterator it=id_to_node_map.begin(); it!= id_to_node_map.end(); ++it)
	{
		it->second->adjmap_x.clear();
		if (it->second->fixed == 0)
		{
		it->second->x = sol_map[it->second->id];
		}
	}

	alist.clear();
	fixlist.clear();
	matrix_mapx.clear();
	sol_map.clear();
	par_map.clear();
}


bool sort_nodes(double a, double b) { return a > b; }
void parse :: HPWL_x()
{
	float HPWL_x = 0;
	int degree;
	float min, max;
	int num;
	int k, present = 0;

        for (map<int,net*>::iterator it=hpl.begin(); it!= hpl.end(); ++it)
        {
                num = 0;
		degree = it->second->lent_x.size();
		vector<double> nodexpos;
		for (map<string,float> :: iterator it1=it->second->lent_x.begin(); it1!= it->second->lent_x.end(); ++it1)
		{
			it1->second = name_to_node_map[it1->first]->x;	
			nodexpos.push_back(it1->second);
		}

		sort(nodexpos.begin(), nodexpos.end(), sort_nodes);
		present = 0;
		
		for (vector<double> :: iterator it2 = nodexpos.begin(); it2!= nodexpos.end(); ++it2)
		{
			num = num + 1;
			if (num == 1)
			{
			HPWL_x = HPWL_x;
			}
			else
			{
			HPWL_x = HPWL_x + fabs(present - *it2);
			}

			present = *it2;
		}
		nodexpos.clear();

	}

	cout << "HPWL_x = " << HPWL_x << endl;
}


void parse :: initsolution_y ()
{
int num1, num2 = 0;
int id2 = 0;
float ay, ap = 0;
float r, p;
float numerator, denominator, numerator1, alpha, beta = 0;
float ynew, rnew, pnew;
//float xdiff, xdiff1;
ofstream myfile ("initsolution_y.txt");
//ofstream myfile1 ("ap.txt");
	for (int i = 0; i < countmove; i++)
        {
        	sol_map[i] = 0;
	}		

	
	for (map<int,matrix_y*> ::iterator it = matrix_mapy.begin(); it!= matrix_mapy.end(); it++)
	{
                ay = 0;
		for (list<pair<int, float> > :: iterator it1 = it->second->movelist.begin(); it1!= it->second->movelist.end(); ++it1)
		{
			ay = float (ay + ((it1->second)*(sol_map[it1->first])));
		}
		r = ((it->second->data) - ay);
		p = r;	
		struct param* g = new param;
                g->r = r;
		g->p = p;
		par_map[it->first] = g;
		//par_map[it->first] = g;
	}

float rmax = 0;
int iteration = 0;
int check = 4;
int iter = 0;
float rtemp = 0;
float res = 0;
do
{
	++iteration;
	rmax = 1;
	check = 0;
	for (map<int,matrix_y*> ::iterator it = matrix_mapy.begin(); it!= matrix_mapy.end(); it++)
	{
		check = check + 1;
		ap = 0;
		for (list<pair<int, float> > :: iterator it1 = it->second->movelist.begin(); it1!= it->second->movelist.end(); ++it1)
                {
                        ap = float (ap + ((it1->second)*(par_map[it1->first]->p)));
		}
		temp_map[it->first] = ap;
	}

	numerator = 0;
	denominator = 0;
	for (map<int,param*> ::iterator it = par_map.begin(); it!= par_map.end(); it++)
	{
		numerator = numerator + ((it->second->r)*(it->second->r));
		denominator = denominator + ((it->second->p)*(temp_map[it->first]));
	}
	alpha = float (numerator/denominator);
	//cout << "numerator is : " << numerator << endl;
	//cout << "denominator is : " << denominator << endl;
	//cout << "alpha is : " << alpha << endl;
	
       
	for (map<int,float> ::iterator it=sol_map.begin(); it!= sol_map.end(); ++it)
        {
		ynew = 0;
                ynew = float ((it->second) + (alpha*(par_map[it->first]->p)));
        	sol_map[it->first] = ynew;
	}
	
	for (map<int,param*> ::iterator it=par_map.begin(); it!= par_map.end(); ++it)
        {
		rnew = 0;
		rnew = float ((it->second->r) - (alpha*(temp_map[it->first])));
		par_map[it->first]->r = rnew;
		res = float (rnew/((matrix_mapy[it->first]->data)));
		res_vec.push_back(pair <int, float> (it->first, res));
                if (fabs(res) > rmax)
                {
                     rmax = fabs(res);
                }

	}

	numerator1 = 0;
	for (map<int,param*> ::iterator it = par_map.begin(); it!= par_map.end(); it++)
	{
		numerator1 = numerator1 + ((it->second->r)*(it->second->r));
	}

	beta = float (numerator1/numerator);


	for (map<int,param*> ::iterator it=par_map.begin(); it!= par_map.end(); ++it)
        {
		pnew = 0;
                pnew = float ((it->second->r) + (beta*(it->second->p)));
        	par_map[it->first]->p = pnew;
	}
	cout << " rmax is : " << rmax << endl;

	temp_map.clear();
	res_vec.clear();

}while ((iteration < 50) || (rmax > 1e-5) );

	cout << "iteration = " << iteration << endl;


        for (map<int,float> ::iterator it=sol_map.begin(); it!= sol_map.end(); ++it)
        {
		myfile << "Node is: " << id_to_node_map[it->first]->name << "\t";
		myfile << "y = " << it->second << endl;
	}


	for (map<string,node*>::iterator it=name_to_node_map.begin(); it!= name_to_node_map.end(); ++it)
	{
		if (it->second->fixed == 0)
		{
		it->second->y = sol_map[it->second->id];
		}
		it->second->adjmap_y.clear();
	}

	for (map<int,node*>::iterator it=id_to_node_map.begin(); it!= id_to_node_map.end(); ++it)
	{
		it->second->adjmap_y.clear();
		if (it->second->fixed == 0)
		{
		it->second->y = sol_map[it->second->id];
		}
	}

	alist.clear();
	fixlist.clear();
	matrix_mapy.clear();
	sol_map.clear();
	par_map.clear();
}


//bool sort_nodes(double a, double b) { return a > b; }
void parse :: HPWL_y()
{
	float HPWL_y = 0;
	int degree;
	float min, max;
	int num;
	int k, present = 0;

        for (map<int,net*>::iterator it=hpl.begin(); it!= hpl.end(); ++it)
        {
                num = 0;
		degree = it->second->lent_y.size();
		vector<double> nodexpos;
		for (map<string,float> :: iterator it1=it->second->lent_y.begin(); it1!= it->second->lent_y.end(); ++it1)
		{
			it1->second = name_to_node_map[it1->first]->y;	
			nodexpos.push_back(it1->second);
		}

		sort(nodexpos.begin(), nodexpos.end(), sort_nodes);
		present = 0;
		
		for (vector<double> :: iterator it2 = nodexpos.begin(); it2!= nodexpos.end(); ++it2)
		{
			num = num + 1;
			if (num == 1)
			{
			HPWL_y = HPWL_y;
			}
			else
			{
			HPWL_y = HPWL_y + fabs(present - *it2);
			}

			present = *it2;
		}

	}

	cout << "HPWL_y = " << HPWL_y << endl;
}

int main (int argc, char* argv[])
{
	parse file;
//	file.parseinputnodes(argv[1]);
	double total_time = 0;
	file.parseinputpl(argv[1]);
	for (int counter = 0; counter < 1; counter = counter + 1)
	{
	file.parsematrix(argv[2]);	
	double begin = cpuTime();
	file.initsolution_x();
	file.HPWL_x();
	 double end_time = cpuTime() - begin;
//	file.initsolution_y();
//	file.HPWL_y();
	//double end_time = cpuTime() - begin;
	total_time = total_time + end_time;
	cout << "time taken = " << total_time << endl;
	}
}

