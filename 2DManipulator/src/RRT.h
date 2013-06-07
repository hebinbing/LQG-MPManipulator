#ifndef _RRT_
#define _RRT_
#define _CRT_RAND_S
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
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"

class RRT{
	
public:
	struct TreeNode{
		Matrix<3> T; //x, y, theta, v
		Matrix<3> u; //v, phi
		int bp;
		bool marked;
		int attempts;

		double depth;
		double clearance;
		std::vector<int> children;
		TreeNode()
		{
			T.reset();
			u.reset();
			bp = -1;
			marked = false;
			depth = 0.0;
			clearance = 9999999;
			attempts = 0;
		}
	};

	struct PathNode{
		Matrix<3> T;
		Matrix<3> u;
	};
	Matrix<3> start;
	Matrix<2> goal;
	double dt;
	double rthreshold;

	double factor;
	double planbias;
	double maxstep;

	int maxchildren;
	int maxattempts;
	int maxiter;
	int maxNNiter;
	int maxtries;
	double plantime;
	double K;
	std::vector<PathNode> rrtpath;
	std::vector<TreeNode> rrttree;
	std::vector<int> rrtpaths;
	std::vector<std::vector<PathNode>> pathSet;

	double w_max;
	double w_min;
	int cal_obstacles;
	int cal_link1;
	int cal_link2;
	int cal_link3;
	int cal_rrt;
	double len; 
	double w;

	RRT(){};
	RRT(const Matrix<3>& Start, const Matrix<2>& gl, const double& timestep, const double& ptime, const double& gr, const int& cal_obs, const int& cal_l1, const int& cal_l2, const int& cal_l3, const int& cal_r){
		start = Start;
		goal = gl;
		dt = timestep;
		plantime = ptime;
		rthreshold = gr;
		pathSet.clear();
		cal_rrt = cal_r;

		Dynamics dy(dt);
		w_max = dy.w_max;
		w_min = dy.w_min;
		len = dy.len;
		w = dy.w;
		cal_obstacles = cal_obs;
		cal_link1 = cal_l1;
		cal_link2 = cal_l2;
		cal_link3 = cal_l3;

	}

	void setPlannerDefaults();
	void initTree(std::vector<TreeNode>& tree, const Matrix<3>& pose);
	double dist(const Matrix<3>& p1, const Matrix<3>& p2);
	int nearestNeighbor(const std::vector<TreeNode>& tree, const Matrix<3>& point);
	Matrix<3> generateControl(const Matrix<3>& p1, const Matrix<3>& p2);
	bool rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths, bool& found,  const clock_t& sc);
	bool buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, const clock_t& sc);
	bool executeRRT(const clock_t& sc);
	bool Plan_K_Seconds();
	void showPath(const int& cal_plan);
};


#endif