#ifndef CELL
#define CELL


#include "point.h"
#include <cmath>
template <class T>
class Cell{
	private:
		int s_index; // the horizontal index of the Cell object
		int t_index; // the vertical index of the Cell object
		T c_bound; // the bound needed to be compared with the Fréchet distance
		Point<T> q_start; // the starting point of the queried line segment
		Point<T> q_end; // the ending point of the queried line segment
		Point<T> i_start; // the starting point of the input line segment
		Point<T> i_end; // the ending point of the input line segment

		double UF1; // the left side of the freespace in the upper edge of the cell
		double UF2; // the right side of the freespace in the upper edge of the cell
		double UR1; // the left side of the reachable space in the upper edge of the cell
		double UR2; // the right side of the reachable space in the upper edge of the cell
		double RF1; // the lower side of the freespace in the right edge of the cell
		double RF2; // the upper side of the freespace in the right edge of the cell
		double RR1; // the lower side of the reachable space in the right edge of the cell
		double RR2; // the upper side of the reachable space in the right edge of the cell

		double LF1;
		double LF2;
		double LR1;
		double LR2;
		double BF1;
		double BF2;
		double BR1;
		double BR2;

		bool test_up; // record whether the upper side of the cell has been tested
		bool test_right; // record whether the right side of the cell has been tested

        // These variables are used to hold values for easier manipulation.
		double a;
		double b;
		double c;
		double d;
		double e;
		double f;

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

        // The function is used to get the determinant of the quadratic equation with one unknown
		double calculate_sqrt_delta(double A, double B, double C){
			if (B*B - 4*A*C >= 0)
				return sqrt(B*B - 4*A*C);
			else
				return -1;
		}
    
        // The function is used to get the first solution of the quadratic equation with one unknown
		double solve_equation_1(double A, double B, double C, double D){
			if(A == 0){
				return -1;
			}

			if(D >= 0)
				return (-B + D)/(2*A);
			else
				return -1;
		}

        // The function is used to get the second solution of the quadratic equation with one unknown
		double solve_equation_2(double A, double B, double C, double D){
			if(A == 0){
				return -1;
			}

			if(D >= 0)
				return (-B - D)/(2*A);
			else
				return -1;
		}

	public:
		Cell(Point<T> q_s, Point<T> q_e, Point<T> i_s, Point<T> i_e, 
				int s_in, int t_in, T c_b){
			s_index = s_in;
			t_index = t_in;
			c_bound=c_b;
			q_start = q_s;
			q_end = q_e;
			i_start = i_s;
			i_end = i_e;
			UF1 = -1; // -1 means there is no free space
			UF2 = -1;
			UR1 = -1; // -1 means there is no reachable space
			UR2 = -1;
			RF1 = -1;
			RF2 = -1;
			RR1 = -1;
			RR2 = -1;
			LF1 = -1;
			LF2 = -1;
			LR1 = -1;
			LR2 = -1;

			BF1 = -1;
			BF2 = -1; 
			BR1 = -1;
			BR2 = -1;

			test_right = false;
			test_up = false;

            
			a = (double)((q_s.x_coo - q_e.x_coo) * (q_s.x_coo - q_e.x_coo)+
					(q_s.y_coo - q_e.y_coo) * (q_s.y_coo - q_e.y_coo));
			b = (double)((i_s.x_coo - i_e.x_coo) * (i_s.x_coo - i_e.x_coo)+
					(i_s.y_coo - i_e.y_coo) * (i_s.y_coo - i_e.y_coo));
			c = (double)(2 * (q_s.x_coo - q_e.x_coo) * (i_e.x_coo - i_s.x_coo)+
					2 * (q_s.y_coo - q_e.y_coo) * (i_e.y_coo - i_s.y_coo));
			d = (double)(2 * (q_e.x_coo - q_s.x_coo) * (q_s.x_coo - i_s.x_coo)+
					2 * (q_e.y_coo - q_s.y_coo) * (q_s.y_coo - i_s.y_coo));
			e = (double)(2 * (i_s.x_coo - i_e.x_coo) * (q_s.x_coo - i_s.x_coo)+
					2 * (i_s.y_coo - i_e.y_coo) * (q_s.y_coo - i_s.y_coo));
			f = (double)((q_s.x_coo - i_s.x_coo) * (q_s.x_coo - i_s.x_coo)+
					(q_s.y_coo - i_s.y_coo) * (q_s.y_coo - i_s.y_coo)-
					c_b * c_b);
		}

		void calculate_LF(){
			double g = calculate_sqrt_delta(b,e,f);
			double t1 = solve_equation_1(b,e,f,g);
			double t2 = solve_equation_2(b,e,f,g);
			double maxv=max(t1,t2);
			double minv=min(t1,t2);
			if(maxv>=0&&minv<=1){
				LF1=max(0,minv);
				LF2=min(1,maxv);
			}
		}

		void calculate_BF(){
			double g = calculate_sqrt_delta(a,d,f);
			double s1 = solve_equation_1(a,d,f,g);
			double s2 = solve_equation_2(a,d,f,g);
			double maxv=max(s1,s2);
			double minv=min(s1,s2);
			if(maxv>=0&&minv<=1){
				BF1=max(0,minv);
				BF2=min(1,maxv);
			}
		}
		void calculate_UF(){
			double g = calculate_sqrt_delta(a,c+d,b+e+f);
			double s1 = solve_equation_1(a,c+d,b+e+f,g);
			double s2 = solve_equation_2(a,c+d,b+e+f,g);
			double maxv=max(s1,s2);
			double minv=min(s1,s2);
			if(maxv>=0&&minv<=1){
				UF1=max(0,minv);
				UF2=min(1,maxv);
			}
		}
		void calculate_RF(){
			double g = calculate_sqrt_delta(b,c+e,a+d+f);
			double t1 = solve_equation_1(b,c+e,a+d+f,g);
			double t2 = solve_equation_2(a,c+d,a+d+f,g);
			double maxv=max(t1,t2);
			double minv=min(t1,t2);
			if(maxv>=0&&minv<=1){
				RF1=max(0,minv);
				RF2=min(1,maxv);
			}
		}
    
        // Initialize UR.
		void calculate_UR_start(){
			UR1 = UF1;
			UR2 = UF2;
		}

        // Initialize RR.
		void calculate_RR_start(){
			RR1 = RF1;
			RR2 = RF2;
		}

        // Uptate UR when going upward
		void calculate_UR_up(){
			calculate_UF();
			if(UF2 >= BR1&&BR1!=-1){
				UR1 = max(UF1, BR1);
				UR2 = UF2;
			}
		}

        // Uptate RR when going to the right
		void calculate_RR_up(){
			calculate_RF();
			RR1 = RF1;
			RR2 = RF2;
		}
    
        // Uptate UR when going to the right
		void calculate_UR_right(){
			calculate_UF();
			UR1 = UF1;
			UR2 = UF2;
		}

        // Uptate RR when going upward
		void calculate_RR_right(){
			calculate_RF();
			if(RF2 >= LR1&&LR1!=-1){
				RR1 = max(RF1, LR1);
				RR2 = RF2;
			}
		}
    
    
        // Set fuctions.
    
		void set_test_up_true(){
			test_up = true;
		}

		void set_test_right_true(){
			test_right = true;
		}

		void set_BR(double br1, double br2){
			BR1 = br1;
			BR2 = br2;
		}

		void set_LR(double lr1, double lr2){
			LR1 = lr1;
			LR2 = lr2;
		}

		void set_UR(double ur1, double ur2){
			UR1 = ur1;
			UR2 = ur2;
		}

		void set_RR(double rr1, double rr2){
			RR1 = rr1;
			RR2 = rr2;
		}

    
        // Get fuctions.
    
		bool up_reach(){
			if (UR1 != -1 && UR2 != -1)
				return true;
			else
				return false;
		}

		bool right_reach(){
			if (RR1 != -1 && RR2 != -1)
				return true;
			else
				return false;
		}

		double get_LF1(){
			return LF1;
		}

		double get_LF2(){
			return LF2;
		}

		double get_BF1(){
			return BF1;
		}

		double get_BF2(){
			return BF2;
		}

		double get_UF1(){
			return UF1;
		}

		double get_UF2(){
			return UF2;
		}

		double get_RF1(){
			return RF1;
		}

		double get_RF2(){
			return RF2;
		}

		double get_LR1(){
			return LR1;
		}

		double get_LR2(){
			return LR2;
		}

		double get_BR1(){
			return BR1;
		}

		double get_BR2(){
			return BR2;
		}

		double get_UR1(){
			return UR1;
		}

		double get_UR2(){
			return UR2;
		}

		double get_RR1(){
			return RR1;
		}

		double get_RR2(){
			return RR2;
		}
		int get_s_index(){
			return s_index;
		}

		int get_t_index(){
			return t_index;
		}

		bool get_test_up(){
			return test_up;
		}

		bool get_test_right(){
			return test_right;
		}

};

#endif
