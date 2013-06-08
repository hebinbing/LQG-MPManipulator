#ifndef _LQGMP_
#define _LQGMP_

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
#include <set>
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Dynamics.h"
#include "RRT.h"

class LQGMP{

public:
	struct PathNodeLQG{
		Matrix<3> T;
		Matrix<3> u; 

		Matrix<3,3> A; 
		Matrix<3,3> B;
		Matrix<3,3> V;
		Matrix<3,3> L;
		Matrix<3,2> K;
		Matrix<2,3> H;

		Matrix<6,6> F;
		Matrix<6,5> G;
		Matrix<6,6> Sigma;

		Matrix<6> y;
		Matrix<6,6> R;

		PathNodeLQG(){
			H.reset();
			T.reset();
			u.reset();
			A.reset();
			B.reset();
			V.reset();
			L.reset();
			K.reset();
			F.reset();
			G.reset();
			Sigma.reset();
			y.reset();
			R.reset();
		}
	};

	struct informationNode{
		int k; //link
		double length; //length measured from the origin of the kth link
		Matrix<2> objpoint;
		int objid;
		double distance;
		Matrix<2> x;
		Matrix<2,3> A;
		Matrix<2> c;

		informationNode(){
			k = 0;
			length = 0;
			objpoint = zeros<2,1>();
			objid = -1;
			distance = -1;
			x = zeros<2,1>();
			A = zeros<2,3>();
			c = zeros<2,1>();
		}
	};

	double dt;
	Matrix<3,3> P0; //initial covariance.
	Matrix<2,2> N; //sense noise
	Matrix<3,3> M; //process noise;
	Matrix<3,3> C;
	Matrix<3,3> D;
	Matrix<5,5> Q;
	std::vector<PathNodeLQG> m_rrtpathLQG;
	std::vector<PathNodeLQG> pathlqg;
	std::vector<RRT::PathNode> m_rrtpath;
	std::vector<PathNodeLQG> m_solutionLQG;
	std::set<int> left_obj;
	double len;

	LQGMP(){};
	LQGMP(const std::vector<RRT::PathNode>& rawpath, const double& timestep, const Matrix<3,3>& initialCov){
		//m_rrtpath.clear();
		m_rrtpath = rawpath;
		dt = timestep;
		P0 = initialCov;

		Dynamics dy(dt);
		len = dy.len;
		
		int size = (int)m_rrtpath.size();
		pathlqg.clear();

		M = identity<3>() * 0.05*0.05;
		N = identity<2>() * 0.05*0.05;
		
		C = identity<3>();
		D = identity<3>();
		int obj1 = 4; int obj2 = 5; int obj3 =6;
		left_obj.insert(obj1);
		left_obj.insert(obj2);
		left_obj.insert(obj3);

		Q.reset();
		Q.insert(0,0, M);
		Q.insert(3,3, N);

	}

	void createABVLK();
	void draw_prior_distribution(const int& cal_ellipse);
	double computeWSConfidence(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	double computeWSConfidence1(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_point, bool& flag);
	double computeConfidence(const Matrix<3>& q, const Matrix<3,3>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	double computeProbability(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point);

	double query(const Matrix<2>& xpos, const Matrix<2,2>& xR, const int& cal_cloneobstacles, const int& cal_point, std::vector<std::pair<Matrix<2>, double>>& cvx);
	double LQGMPTruncation(const Matrix<3>& qpos, const Matrix<3,3>& qR, const int& cal_cloneobstacles, const int& cal_point, std::vector<std::pair<Matrix<3>, double>>& ccvx);
	double computeLQGMPTruncation(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point);

	void queryalldistance(const Matrix<3>& qpos, const Matrix<3,3>& qR, const int& cal_obstacles, const int& cal_environment, const int& cal_point, std::vector<informationNode>& Distance);
	int pick_min_dis(std::vector<informationNode>& Distance);
	bool check_contact_point(const std::vector<std::pair<Matrix<2>, double>>& cvx, const Matrix<2>& objpoint);
	void DistanceToConstraints(const Matrix<3>& qpos, const Matrix<3,3>& qR, std::vector<std::pair<Matrix<3>, double>>& qcvx, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	double EstimateAndTruncate(const Matrix<3,3> initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	void draw_truncate_distribution(const int& cal_trunc_ellipse);
};


#endif