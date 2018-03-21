#ifndef QUERY
#define QUERY
#include <fstream>
#include <map>
#include "trajectory.h"
#include <string>
#include <iostream>
#include <cmath>
#include "freespace.h"

using namespace std;
//template <class T>
vector<string> result_data;
bool MySearchCallback(string id, void* arg)
{
	result_data.push_back(id);
	return true; // keep going
}
vector<string> temp_result;
bool MySearchCallback2(string id, void* arg)
{
	temp_result.push_back(id);
	return false; // keep going
}
class Query{
	private: 
		map<string, Trajectory<double> *> * tra_list; // the list of the input trajectories
		string dataname;
		double bound;
		RTree<string, double, 4, double>* se_tree; // the R-tree of the starting and ending points of each input trajectory
		RTree<string, double, 2, double> my_tree; // the R-tree of the queried trajectory
		bool state;
		int output;
		double max(double a, double b){
			if (a > b)
				return a;
			else
				return b;
		}
	public:
		Query(map<string, Trajectory<double> *> * t_l, string dn, double b,
				RTree<string, double, 4, double>* se){
			se_tree=se;
			tra_list = t_l;
			state = true;
			bound=b;
			dataname=dn;
			output=0;
			if(tra_list -> find(dataname)==tra_list->end()){
				cout<<"Error query. No such trajectory "<<dataname<<"\n";
				state = false;
			}
			else{
				Trajectory<double>* tra_queried = (tra_list -> find(dataname))->second;
				int index=0;
                
                // Constuct the R-tree of the queried trajectory.
				while(index<tra_queried->get_size()){
					double a_point[2];
					a_point[0] = tra_queried ->get_point(index).x_coo;
					a_point[1] = tra_queried ->get_point(index).y_coo;
					my_tree.Insert(a_point,a_point,dataname);
					index++;
				}
                
				int tra_size = tra_queried->get_size();
				
                // Define the searching range of se_tree.
                double mi_s[4];
				mi_s[0]=tra_queried->get_point(0).x_coo-bound;
				mi_s[1]=tra_queried->get_point(0).y_coo-bound;
				mi_s[2]=tra_queried->get_point(tra_size-1).x_coo-bound;
				mi_s[3]=tra_queried->get_point(tra_size-1).y_coo-bound;
				double ma_s[4];
				ma_s[0]=tra_queried->get_point(0).x_coo+bound;
				ma_s[1]=tra_queried->get_point(0).y_coo+bound;
				ma_s[2]=tra_queried->get_point(tra_size-1).x_coo+bound;
				ma_s[3]=tra_queried->get_point(tra_size-1).y_coo+bound;
                
                // Use se_tree to select trajectories whose starting and ending points satisfying the condition.
				se_tree->Search(mi_s,ma_s,MySearchCallback,NULL);

				string result_number_str = "";
				int i = 0;
				while(true){
					//Error
					if (dataname[i] == '-'||dataname[i]=='\0')
						break;
					i++;
				}
				i++;
				while(true){
					//Error
					if (dataname[i] == '.'||dataname[i]=='\0')
						break;
					result_number_str = result_number_str + dataname[i];
					i++;
				}
				string result="result-"+result_number_str+".txt";
				ofstream ofs(result.c_str());
				int k=0;
				while(k < result_data.size()){
					Trajectory<double>* tra_input=tra_list->find(result_data[k])->second;
					if(tra_input==tra_queried){
						k++;
						continue;
					}
					int input_size = tra_input->get_size();
					double max_distance = max(tra_input->get_max_distance(),tra_queried->get_max_distance());
					double new_bound = bound + max_distance;
					int l=0;
					while(l<5){
						l++;
						int index2=input_size*l/6;
						double temp_x=tra_input->get_point(index2).x_coo;
						double temp_y=tra_input->get_point(index2).y_coo;
                        
                        // Define the searching range of my_tree.
						double mi_s2[2];
						mi_s2[0]=temp_x-new_bound;
						mi_s2[1]=temp_y-new_bound;
						double ma_s2[2];
						ma_s2[0]=temp_x+new_bound;
						ma_s2[1]=temp_y+new_bound;
                        
                        // Use my_tree to and five points of each input trajectory to decide whether the input trajectory is one of the expected results.
						my_tree.Search(mi_s2,ma_s2,MySearchCallback2,NULL);
						if(temp_result.size()==0){
							state = false;
							break;
						}
						else{
							temp_result.clear();
						}
					}
					if (state==true){
						Freespace <double> test(tra_queried, tra_input, bound, 
								tra_queried->get_size(),
								tra_input->get_size());
						bool result = test.algorithm();
						if(result){
							output++;
							ofs << tra_list->find(result_data[k])->first << '\n';
							//cout<<dataname<<'\n';
						} 
					}
					else{
						state = true;
					}
					k++;
				}
				if(output==0){
					ofs << "No trajectory fulfilled the given query range" << '\n';
					//cout<<dataname<<'\n';
				}
				ofs.close();
				result_data.clear();
			} 
		} 
};

#endif
