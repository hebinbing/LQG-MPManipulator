#ifndef _DYNAMICS_
#define _DYNAMICS_


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

class Dynamics{

public:
	double tau;
	double w_max;
	double w_min;
	double len;
	double w;

	Dynamics(const double& dt)
	{
		len = 1.5;
		w = 0.3;
		tau = dt;
		w_max = 0.4;
		w_min = -0.4;
	}


	Matrix<3> dynamics_zero_noise(const Matrix<3>& X, const Matrix<3>& u){
		Matrix<3> input = u;
		if(input[0] > w_max)
			input[0] = w_max;
		if(input[1] > w_max)
			input[1] = w_max;
		if(input[2] > w_max)
			input[2] = w_max;

		if(input[0] < w_min)
			input[0] = w_min;
		if(input[1] < w_min)
			input[1] = w_min;
		if(input[2] < w_min)
			input[2] = w_min;

		Matrix<3> nextX = zeros<3,1>();
		nextX[0] = X[0] + tau * input[0];
		nextX[1] = X[1] + tau * input[1];
		nextX[2] = X[2] + tau * input[2];

		for(int i = 0; i < 3; i++){
			if(nextX[i] > 2*3.1416)
				nextX[i] = nextX[i] - 2*3.1416;
			if(nextX[i] < 0)
				nextX[i] = nextX[i] + 2*3.1416;
		}


		return nextX;
	} 

	Matrix<3> dynamics_noise(const Matrix<3>& X, const Matrix<3>& u, const Matrix<3,3> Noise){
		Matrix<3> input = u;
		if(input[0] > w_max)
			input[0] = w_max;
		if(input[1] > w_max)
			input[1] = w_max;
		if(input[2] > w_max)
			input[2] = w_max;

		if(input[0] < w_min)
			input[0] = w_min;
		if(input[1] < w_min)
			input[1] = w_min;
		if(input[2] < w_min)
			input[2] = w_min;

		Matrix<3,1> noise = sampleGaussian(zeros<3,1>(), Noise);
		input += noise;

		Matrix<3> nextX = zeros<3,1>();
		nextX[0] = X[0] + tau * input[0];
		nextX[1] = X[1] + tau * input[1];
		nextX[2] = X[2] + tau * input[2];

		for(int i = 0; i < 3; i++){
			if(nextX[i] > 2*3.1416)
				nextX[i] = nextX[i] - 2*3.1416;
			if(nextX[i] < 0)
				nextX[i] = nextX[i] + 2*3.1416;
		}

		return nextX;
	}

	Matrix<2> observe_zero_noise(const Matrix<3>& X)
	{
		Matrix<2> ob = zeros<2,1>();
		ob[0] = X[0];  //base;

		Matrix<2> p = forward_k(X, 3, len);
		ob[1] = p[1]; //y position of the endeffector.
		return ob;
	}
	Matrix<2> observe_noise(const Matrix<3>& X, const Matrix<2,2>& Noise)
	{
		Matrix<2> ob = zeros<2,1>();
		ob[0] = X[0];  //base;

		Matrix<2> p = forward_k(X, 3, len);
		ob[1] = p[1]; //y position of the endeffector.

		Matrix<2> noise = sampleGaussian(zeros<2,1>(), Noise);
		return ob + noise;
	}

};



#endif