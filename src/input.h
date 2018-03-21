#ifndef INPUT
#define INPUT
#include <fstream>
#include <map>
#include "trajectory.h"
#include <string>
#include <cstdio>
#include "RTree.h"
using namespace std;
//template <class T>

typedef string ValueType;

struct Rect
{
	Rect()  {}

	Rect(int a_minX, int a_minY, int a_maxX, int a_maxY)
	{
		min[0] = a_minX;
		min[1] = a_minY;

		max[0] = a_maxX;
		max[1] = a_maxY;
	}
	int min[2];
	int max[2];
};



class Input{
	private:
		map <string, Trajectory<double> *> dataset;
		RTree<string, double, 4, double> se_tree;// the R-tree of the starting and ending points of each input trajectory
		string filename;
		bool state;
		double minx;
		double maxx;
		double miny;
		double maxy;
		double min(double a, double b){
			if (a < b)
				return a;
			else
				return b;
		}

		double max(double a, double b){
			if (a > b)
				return a;
			else
				return b;
		}

	public:
		Input(string fn="dataset.txt"){
			double minx=0;
			double maxx=0;
			double miny=0;
			double maxy=0;
			double firsty=0;
			double firstx=0;
			state = true;
			filename = fn;
			ifstream f;
			f.open(filename.c_str());
			if(f.is_open()){
				time_t before_input;
				time(&before_input);
				int d=0;
				while(!f.eof()){
					string dataname;
					getline(f,dataname);
					if(dataname.size()<1){
						break;
					}
					ifstream t;
					t.open(dataname.c_str());

					if(t.is_open()){
						int i = 0;
						Trajectory<double>* tra =new Trajectory<double>;
						while(!t.eof()){
							string line="";
							if(i == 0){
								getline(t,line);
								i=1;
								continue;
							}
							getline(t,line);
							if(line.size()<1){
								break;
							}
							string x_string="";
							string y_string="";
							int j=0;
							while(line[j]!=' '){
								x_string=x_string+line[j];
								j++;
							}
							j++;
							while(line[j]!=' '){
								y_string=y_string+line[j];
								j++; 
							}
							double x=std::stod(x_string);
							double y=std::stod(y_string);
							minx=min(minx,x);
							maxx=max(maxx,x);
							miny=min(miny,y);
							maxy=max(maxy,y);
							if(d==0&&i==1){
								firstx=x;
								firsty=y;
								i=2;
							}
							tra->push_point(x-firstx,y-firsty);
						}
						t.close();
						tra->calculate_max_distance();
						dataset.insert(std::pair<string,Trajectory<double> * >(dataname,tra));
						double spoint[4];
						spoint[0] = tra->get_point(0).x_coo;
						spoint[1] = tra->get_point(0).y_coo;
						int tra_size=tra->get_size();
						spoint[2] = tra->get_point(tra_size-1).x_coo;
						spoint[3] = tra->get_point(tra_size-1).y_coo;
                        
                        //Insert the starting point and the ending point of each input trajectory into se_tree.
						se_tree.Insert(spoint,spoint,dataname);
					}
					else{
						cout<<"Error open "<<dataname<<"\n";
					}
					d++;
				}
				f.close();
				time_t after_input;
				time(&after_input); 
				double seconds_i = difftime(after_input,before_input);
				printf ("%.f seconds for load data.\n", seconds_i);
			}
			else{
				state = false;
				cout<<"Error opening input file.\n";
			}
		}

		~Input(){
			dataset.clear();
		}

		map <string, Trajectory<double> *>* get_dataset(){
			return &dataset;
		}

		bool get_state(){
			return state;
		}

		RTree<string, double, 4, double>* get_se_tree(){
			return &(se_tree);
		}
};

#endif
