///////////////////
//This code is used to implement the idea of real-time replanning. 
//namely, replanning ahead. 
//assume the time for replanning is also dt
///////////////////////////////////

#define _CRT_RAND_S

#define DIM 2
#define X_DIM 4  //X = (x,y,theta, v)  theta: orientation, v speed.
#define U_DIM 2  //a, phi, a: acceleration, phi: steering wheel angel
#define Z_DIM 3  //it can observe position x and position y.  
#define INFTY 9e9


//define constratins:


const static double goal_radius = 0.5;
const static double plan_goal_radius = 0.6;
const static double dt = 0.5;  // asume the time for replan is also dt. 
const static double car_l = 1.0;
const static double MaxPlanTime = 1.0;


#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>
#include <windows.h>
#include <algorithm>
#include <ppl.h>
#include <omp.h>
#include <set>
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "RRT.h"
#include "LQGMP.h"

using namespace Concurrency;
//using namespace std;

static Matrix<2> goal;
static Matrix<3> start;
static Matrix<3, 3> P0; //covariance

/***********************Environment settings***************************/
int cal_environment;
int cal_line;
int cal_start;
int cal_obstacles;
int cal_goal;
int cal_rrt;
int cal_paths;
int cal_ellipse;
int cal_ellipse_trunc;
int cal_point;
int cal_cienvironment, cal_ciobstacles, cal_cipoint;
int cal_plan;
int cal_execute;
int cal_robot;
int cal_link1;
int cal_link2;
int cal_link3;
int cal_visual1;
int cal_visual2;
int cal_visual3;
static std::vector<Matrix<2,1>> obstaclesSet; //ax + by = b
/******************End of Environment settings************************/


double Random()
{
	double a;
	a = rand()%1000;
	return a / 1000;
}

/**********************init environment***************************/
void initEnvironment() 
{
	goal[0] = 4.0;
	goal[1] = -1.0;

	start[0] = 1/4.0 * M_PI; start[1] = 0.0; start[2] = 0;
	P0 = identity<3>() * 0.001;
	// callisto params
	CAL_SetViewParams(0, 0, 0, 10, 0, 0, 0, 0, 1, 0); 

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 0, 0, 1, 1);

	
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 1, 0, 1);
	
	CAL_CreateGroup(&cal_line, 0, false, "Line");
	CAL_SetGroupColor(cal_line, 0, 0, 0);

	float line[6] = {-5.0, 5.0, 0.0, 5.0, 5.0, 0.0};
	int np[1] = {2};
	CAL_CreatePolyline(cal_line, 1, np, line);
	float line1[6] = {5.0, 5.0, 0.0, 5.0, -5.0, 0.0};
	CAL_CreatePolyline(cal_line, 1, np, line1);
	float line2[6] = {5.0, -5.0, 0.0, -5.0, -5.0, 0.0};
	CAL_CreatePolyline(cal_line, 1, np, line2);
	float line3[6] = {-5.0, -5.0, 0.0, -5.0, 5.0, 0.0};
	CAL_CreatePolyline(cal_line, 1, np, line3);

	int obj1 = -1; int obj2 = -1; int obj3 = -1;
	CAL_CreateBox(cal_obstacles, 1.0, 1.0, 0.02, -2.5, 2.5, 0.0, &obj1);
	CAL_CreateBox(cal_obstacles, 1.0, 1.0, 0.02, -2.5, -2.5, 0.0, &obj2);
	CAL_CreateBox(cal_obstacles, 1.0, 0.5, 0.02, 2.5, 0.0, 0.0, &obj3);
	std::cout<<obj1<<" "<<obj2<<" "<<obj3<<std::endl;


	CAL_CreateGroup(&cal_paths, 0, true, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 1);

	CAL_CreateGroup(&cal_link1, 0, true, "link1");
	CAL_SetGroupColor(cal_link1, 1, 0, 0);

	CAL_CreateGroup(&cal_link2, 0, true, "link2");
	CAL_SetGroupColor(cal_link2, 1, 0, 0);

	CAL_CreateGroup(&cal_link3, 0, true, "link3");
	CAL_SetGroupColor(cal_link3, 1, 0, 0);
	CAL_SetGroupVisibility(cal_link1, 0, false);
	CAL_SetGroupVisibility(cal_link2, 0, false);
	CAL_SetGroupVisibility(cal_link3, 0, false);

	CAL_CreateGroup(&cal_visual1, 0, false, "visual1");
	CAL_SetGroupColor(cal_visual1, 1, 0, 0);

	CAL_CreateGroup(&cal_visual2, 0, false, "visual2");
	CAL_SetGroupColor(cal_visual2, 0, 1, 0);

	CAL_CreateGroup(&cal_visual3, 0, false, "visual3");
	CAL_SetGroupColor(cal_visual3, 0, 0, 1);

	//VisualManipulator(start, cal_visual1, cal_visual2, cal_visual3);
	CreateManipulator(start, cal_link1, cal_link2, cal_link3, true); //true for construction.

	Matrix<2> p = forward_k(start, 3, 1.5);
	std::cout<<p[0]<<" "<<p[1]<<std::endl;
	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0, 0, 0, 0);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 0, 0.0, 0.4);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "Ellipse_trunc");
	CAL_SetGroupColor(cal_ellipse_trunc, 1, 0, 0);

	CAL_CreateGroup(&cal_execute, 0, false, "execute");
	CAL_SetGroupColor(cal_execute, 1, 0, 0);
	
	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	CAL_SetGroupColor(cal_goal, 0, 1, 1, 0.5);
	CAL_CreateCylinder(cal_goal, (float) goal_radius, 0.05f, 0, 0, 0);
	CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], -0.025f);
	CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);

	//CAL_CreateSphere(cal_point, 0.2, 7.0/4.0*3.1416, 0.0, 0.0);
	
}


void showpath(const std::vector<RRT::PathNode>& path, int group_id1, int group_id2, int group_id3){
	
	
	for(int i = 0; i < (int)path.size(); i++){
		Matrix<3> state = path[i].T;
		Matrix<2> p1 = forward_k(state, 1, 1.5);
	//	std::cout<<p1[0]<<" "<<p1[1]<<std::endl;
		float line1[6] = {0.0, 0.0, 0.0, p1[0], p1[1], 0.0};
		int np[1] = {2};
		//CAL_SetGroupColor(cal_rrt, 1.0, 0.0, 0.0);
		CAL_CreatePolyline(group_id1, 1, np, line1);

		Matrix<2> p2 = forward_k(state, 2, 1.5);
		float line2[6] = {p1[0], p1[1], 0.0, p2[0], p2[1], 0.0};
		//CAL_SetGroupColor(cal_rrt, 0.0, 1.0, 0.0);
		CAL_CreatePolyline(group_id2, 1, np, line2);

		Matrix<2> p3 = forward_k(state, 3, 1.5);
		float line3[6] = {p2[0], p2[1], 0.0, p3[0], p3[1], 0.0};
		//CAL_SetGroupColor(cal_rrt, 0.0, 0.0, 1.0);
		CAL_CreatePolyline(group_id3, 1, np, line3);

	}

}

void showconfigurationpath(const std::vector<RRT::PathNode>& path)
{
	Dynamics dy(dt);
	for(int i = 0; i < (int)path.size()- 1; i++){
		Matrix<3> state = path[i].T;
		Matrix<3> u = path[i].u;
		Matrix<3> nextState = dy.dynamics_zero_noise(state, u);
		float line[6] = {state[0], state[1], state[2], nextState[0], nextState[1], nextState[2]};
		int np[1] = {2};
		CAL_CreatePolyline(cal_visual1, 1, np, line);
	}
}


/*********************End of environment init**************************/


int main()
{

	srand(1000);
	CAL_Initialisation (true, true, true);
	initEnvironment();
	
	double dt = 0.5;
	double plantime = 10;
	P0 = identity<3>() * 0.05*0.05;

	RRT rrt(start, goal, dt, plantime, goal_radius, cal_obstacles, cal_link1, cal_link2, cal_link3, cal_rrt);
	rrt.setPlannerDefaults();
	rrt.Plan_K_Seconds();
	std::cout<<rrt.pathSet.size()<<std::endl;
	showpath(rrt.pathSet[0], cal_visual1, cal_visual2, cal_visual3);
	showconfigurationpath(rrt.pathSet[0]);

	LQGMP lqgmp(rrt.pathSet[0],dt, P0); 
	double prob = lqgmp.computeProbability(P0, cal_obstacles, cal_environment, cal_point);
	lqgmp.draw_prior_distribution(cal_ellipse);
	std::cout<<"Probablity of success (LQGMP): "<<prob<<std::endl;

	double prob2 = lqgmp.EstimateAndTruncate(P0, cal_obstacles, cal_environment, cal_point);
	std::cout<<"Probability of success (Truncated LQGMP): "<<prob2<<std::endl;
	lqgmp.draw_truncate_distribution(cal_ellipse_trunc);
	int num;
	std::cin>>num;

	CAL_End();
	return 0;
}