#include "freespace.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include "input.h"
#include "query.h"
#include <vector>
# include <omp.h>
# include <thread>

using namespace std;
struct query_node{
string dataname;
double bound;
};

int main(int argv, char **args)
{
	string query_file;
	string input_file;
	if(argv==3){
		query_file = args[2];
		input_file = args[1];
	}
	else if(argv > 3 || argv == 2){
		cout<<"Error Argument. Please view Readme.\n";
		return 0;
	}
	else{
		input_file = "dataset.txt";
		query_file = "queries.txt";
	} 
	Input input(input_file); 
	map <string, Trajectory<double> *> * dataset=input.get_dataset(); // load the list of input trajectories
	RTree<string, double, 4, double>* se_tree=input.get_se_tree(); // load the se_tree
	bool query_state = true;
	ifstream q;
	time_t before_query;
	time(&before_query);
	q.open(query_file.c_str());
	vector<query_node> query_list;
	if(q.is_open()){
		while(!q.eof()){
			string line;
			string dataname = "";
			double bound;
			getline(q,line);
			if(line.size()<1){
				break;
			}
			int i = 0;
			while(line[i] != ' '){
				dataname = dataname + line[i];
				i++;
			}
			while(line[i] == ' '){
				i++;
			}
			
			string bound_str = "";
			while(line[i] != '\n'&&line[i]!='\0'&&' '){
				bound_str = bound_str + line[i];
				i++;
			}
			bound = std::stod(bound_str);
			struct query_node node = {.dataname=dataname, .bound=bound};
			query_list.push_back(node);
			//cout<<bound<<'\n';
			//parallel
		}
		q.close();
	}
	else{
		//query_state = false;
		cout<<"Error opening query file.\n";
		return 0;
	}
	int k;
	unsigned cpucore = std::thread::hardware_concurrency();
	if (cpucore==0){
		cpucore=16;
		printf ("Core number not find, use 16 as Threads number.\n");
	}
	else{
		printf ("Find %d Core. Use %d as Threads number.\n",cpucore,cpucore);
	}
	#pragma omp parallel shared(dataset,query_list,se_tree,chunk) private(k)
	{
		#pragma omp for schedule(guided) nowait num_threads(cpucore)
		for(k=0;k<query_list.size();k++) {
			Query query(dataset, query_list[k].dataname, query_list[k].bound, se_tree);
			cout<<query_list[k].dataname<<" is done.\n";
		}
	}
	time_t after_query;
	time(&after_query); 
	double seconds_q = difftime(after_query,before_query);
	printf ("%.f seconds for query.\n", seconds_q);
	return 0;
}



