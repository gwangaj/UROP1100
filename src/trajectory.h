#ifndef TRAJECTORY
#define TRAJECTORY

#include <vector>
#include "point.h"
#include <cmath>
using namespace std;

template <class T>
class Trajectory{
	private:
		vector < Point<T> > points;
		T max_distance_sq; // The square of the longest line segment of the Trajectory object
		T max_distance; // The longest line segment of the Trajectory object
		T max(T a, T b){
			if (a > b)
				return a;
			else
				return b;
		}
	public:
		Trajectory(){
			max_distance_sq=0;
			max_distance=0;
		}
		~Trajectory(){
			points.clear();
		}
		void push_point(T x, T y){
			if(points.size()==0){
				Point<T> p;
				p.x_coo=x;
				p.y_coo=y;
				//p.number=k;
				//p.t_id=t;
				points.push_back(p);}
			else{
				T back_x = points.back().x_coo;
				T back_y = points.back().y_coo;
				Point<T> p;
				p.x_coo=x;
				p.y_coo=y;
				T distance_sq=(back_x-x)*(back_x-x)+(back_y-y)*(back_y-y);
				max_distance_sq=max(max_distance_sq, distance_sq);
				//p.number=k;
				//p.t_id=t;
				points.push_back(p);
			}

		}
    
        // Get the Point object given the index.
		Point<T> get_point(int index){
			return points[index];
		}
    
        // Get the total number of Point objects.
		int get_size(){
			return points.size();
		}
    
		T get_max_distance(){
			return max_distance;
		}
		void calculate_max_distance(){
			max_distance=sqrt(max_distance_sq);
		}
};

#endif
