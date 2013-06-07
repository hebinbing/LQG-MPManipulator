#include "LQGMP.h"


void LQGMP::createABVLK()
{
	int ns = (int) m_rrtpath.size();	
	pathlqg.clear();
	for(int i = 0; i < ns; i++){
		PathNodeLQG tmpnode;
		tmpnode.T = m_rrtpath[i].T;
		tmpnode.u = m_rrtpath[i].u;
		pathlqg.push_back(tmpnode);
	}
	for(int i = 0; i < ns; i++){
		RRT::PathNode tmp = m_rrtpath[i];
		Matrix<3> u = tmp.u;
		Matrix<3> X = tmp.T;

		Matrix<3,3> A;
		A.reset();
		Matrix<3,3> B;
		B.reset();
		Matrix<3,3> V;
		V.reset();
		A = identity<3>();
		B = dt * identity<3>();

		Matrix<2,3> H; H.reset();
		H(0,0) = X[0];			H(0,1) = 0;				H(0,2) = 0; 
		H(1,0) = len*cos(X[0]); H(1,1) = len*cos(X[1]); H(1,2) = len*cos(X[2]);

		V = B;
		pathlqg[i].A = A;
		pathlqg[i].B = B;
		pathlqg[i].V = V;
		pathlqg[i].H = H;
	}

	Matrix<3,3> S = C;
	int length = (int)pathlqg.size() - 1;
	for(int k = length - 1; k != -1; k--){
		pathlqg[k].L = -!(~pathlqg[k].B*S*pathlqg[k].B + D)*~pathlqg[k].B*S*pathlqg[k].A;
		S = C + ~pathlqg[k].A*S*pathlqg[k].A + ~pathlqg[k].A*S*pathlqg[k].B*pathlqg[k].L;
	}

	Matrix<3,3> P = P0;
	for(int k = 1; k <= length; k++){
		P = pathlqg[k-1].A*P*~pathlqg[k-1].A + pathlqg[k-1].V * M *~pathlqg[k-1].V;
		
		pathlqg[k].K = P * ~pathlqg[k].H* !(pathlqg[k].H*P*~pathlqg[k].H + N);
		P = (identity<3>() - pathlqg[k].K * pathlqg[k].H) * P;
	}
	pathlqg[0].Sigma = zeros<6,6>();
	pathlqg[0].Sigma.insert(0,0, P0);

	for(int k = 1; k <= length; k++){
		Matrix<6,6> F;
		F.insert(0,0, pathlqg[k-1].A);
		F.insert(0,3, pathlqg[k-1].B * pathlqg[k-1].L);
		F.insert(3,0, pathlqg[k].K * pathlqg[k].H * pathlqg[k-1].A);
		F.insert(3,3, pathlqg[k-1].A + pathlqg[k-1].B * pathlqg[k-1].L - pathlqg[k].K * pathlqg[k-1].H * pathlqg[k-1].A);
		pathlqg[k-1].F = F;

		Matrix<6, 5> G;
		G.insert(0,0, pathlqg[k-1].V);
		G.insert(0,3, zeros<3,2>());
		G.insert(3,0, pathlqg[k].K * pathlqg[k-1].H * pathlqg[k-1].V);
		G.insert(3,3, pathlqg[k].K);
		pathlqg[k-1].G = G;
		
		Matrix<5,5> Q;
		Q.reset();
		Q.insert(0,0, M);
		Q.insert(3,3, N);
		pathlqg[k].Sigma = F * pathlqg[k-1].Sigma * ~F + G * Q * ~G;
	}
}



void LQGMP::draw_prior_distribution(const int& cal_ellipse){

	createABVLK();
	for(int i = 0; i < (int)pathlqg.size(); i++){
		drawEllipse3d(pathlqg[i].T, pathlqg[i].Sigma.subMatrix<3,3>(0,0), cal_ellipse, true);
	}

}

//return workspace d
double LQGMP::computeWSConfidence(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	Matrix<2,2> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~EVec);
	Matrix<4,1> q = quatFromRot(Temp);
	for(int i = 0; i < 2; i++){
		EVal(i,i) += 0.0001;
	}

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

	Matrix<2,2> invScale = zeros<2,2>();
	invScale(0,0) = 1/(float)sqrt(EVal(0,0));
	invScale(1,1) = 1/(float)sqrt(EVal(1,1));
	Matrix<2> transPos =  invScale * ~EVec * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);

	int num_pairs = -1;
	CAL_GetClosestPairs(cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	if(num_pairs == 0){
		CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
		CAL_SetGroupScaling(cal_environment, 1,1,1);
		return -9999;
	}

	double distance = results[0].distance;
	delete [] results;

	CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
	CAL_SetGroupScaling(cal_environment, 1,1,1);

	return distance;
}

double LQGMP::computeWSConfidence1(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_point, bool& flag)
{
	Matrix<2,2> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~EVec);
	Matrix<4,1> q = quatFromRot(Temp);
	for(int i = 0; i < 2; i++){
		EVal(i,i) += 0.0001;
	}

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_obstacles, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

	Matrix<2,2> invScale = zeros<2,2>();
	invScale(0,0) = 1/(float)sqrt(EVal(0,0));
	invScale(1,1) = 1/(float)sqrt(EVal(1,1));
	Matrix<2> transPos =  invScale * ~EVec * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);

	int num_pairs = -1;
	CAL_GetClosestPairs(cal_point, cal_obstacles, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	if(num_pairs == 0){
		CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
		CAL_SetGroupScaling(cal_obstacles, 1,1,1);
		flag = false; //no obstacles.
		return 999999999.0;
	}

	double distance = results[0].distance;
	delete [] results;

	CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
	CAL_SetGroupScaling(cal_obstacles, 1,1,1);

	return distance;
}

double LQGMP::computeConfidence(const Matrix<3>& q, const Matrix<3,3>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	double mindis = 10000;
	int mlink = -1; 
	int mpoint = -1;
	double seg = 5;
	double seglen = len / seg;
	//find the maximum likelihodd collision point. the one that has the minimum times of standard deviation (1 after transformation.)
	for(int link = 1; link <=3; link++){
		for(int point = 0; point < seg; point++){
			double length = (point+1)*seglen; //length on the particular link.
			Matrix<2,3> Joc = zeros<2,3>();
			forward_jocabi(q, link, length, Joc);
			Matrix<2> c = forward_k(q, link, length) - Joc * q;
			Matrix<2,3> A = Joc; //A*q + c
			Matrix<2> xmean = forward_k(q, link, length);
			Matrix<2,2> xcov = A * R * ~A;
			double tmpdis = computeWSConfidence(xmean, xcov, cal_obstacles, cal_environment, cal_point); 
			if(tmpdis < mindis){
				mindis = tmpdis;
				mlink = link;
				mpoint = point;
			}
		} 
	}
	return mindis;
}

double LQGMP::computeProbability(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	P0 = initialCov;
	createABVLK();
	double log_prob = 0;
	for(int i = 0; i < (int)pathlqg.size(); i++)
	{
		Matrix<3,1> q = pathlqg[i].T;
		Matrix<3,3> qR = pathlqg[i].Sigma.subMatrix<3,3>(0,0);
		double conf = computeConfidence(q, qR, cal_obstacles, cal_environment, cal_point);
		double p = log(incompletegamma(1, 0.5*conf*conf));
		log_prob += p;
	}

	return exp(log_prob);
}

//this query only query for one xpos (the mostlikely-in-collision sample on the links), only query once. clone obstacles are the cloned and left obstacles.
double LQGMP::query(const Matrix<2>& xpos, const Matrix<2,2>& xR, const int& cal_cloneobstacles, const int& cal_point, std::vector<std::pair<Matrix<2>, double>>& cvx)
{
	int ps = 1; //return the probablity of collisioin free for this constraint.
	Matrix<2,2> EVec, EVal;
	jacobi(xR, EVec, EVal);
	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~EVec);
	Matrix<4,1> q = quatFromRot(Temp);
	for(int i = 0; i < 2; i++){
		EVal(i,i) += 0.0001; //smooth it.
	}
	//transform the cloned obstacles.
	CAL_SetGroupQuaternion(cal_cloneobstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_cloneobstacles, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);
	Matrix<2,2> invScale = zeros<2,2>();
	invScale(0,0) = 1/(float)sqrt(EVal(0,0));
	invScale(1,1) = 1/(float)sqrt(EVal(1,1));
	Matrix<2,2> T = invScale * ~EVec;
	Matrix<2,2> invT = !T;
	Matrix<2> transPos =  invScale * ~EVec * xpos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0.0);
	int num_pairs;
	CAL_GetClosestPairs(cal_point, cal_cloneobstacles, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	Matrix<2> a = zeros<2,1>();
	double b = 0;
	if(num_pairs != 0){
		//just need the first result.
		Matrix<2> op = zeros<2,1>();
		op[0] = results[0].vector0[0];
		op[1] = results[0].vector0[1];
		Matrix<2> objp = zeros<2,1>();
		objp[0] = results[0].vector1[0];
		objp[1] = results[0].vector1[1];
		op = invT*op; objp = invT*objp;
		Matrix<2> a = objp - op;
		a = a / sqrt(tr(~a*a));
		b = tr(~a * objp);
		
		//statistics for ax
		double mean = tr(~a*xpos);
		double variance = tr(~a*xR*a) + 0.00001;
		ps = cdf((b - mean) / sqrt(variance));

		cvx.push_back(std::make_pair(a, b)); //this is the workspace constraint.

		//j = 0: truncate the object itself first. 
		for(int j = 0; j < num_pairs; j++)
		{
			Matrix<2> objpos = zeros<2,1>();
			objpos[0] = results[j].vector1[0];
			objpos[1] = results[j].vector1[1];
			objpos = invT*objpos;
			if(tr(~a*objpos) > b - 0.00001)
				CAL_DestroyObject(results[j].objID1);
		}
	}
	delete [] results;
	CAL_SetGroupQuaternion(cal_cloneobstacles, 0, 0, 0, 1);
	CAL_SetGroupScaling(cal_cloneobstacles, 1, 1, 1);
	return ps;
}

//this is for a configuration spacce truncation for one time step.it needs to keep truncating until there is no obstalces.
double LQGMP::LQGMPTruncation(const Matrix<3>& qpos, const Matrix<3,3>& qR, const int& cal_cloneobstacles, const int& cal_point, std::vector<std::pair<Matrix<3>, double>>& ccvx)
{
	double ps = 0;
	
	while(true){
		
		bool flag = true;
		
		//find the mostlikly-in-collision sample on the links.
		double mindis = 10000; //measure the number of standard deviations.
		int mlink = -1; 
		int mpoint = -1;
		double seg = 5;
		double seglen = len / seg;
		Matrix<2> minmean = zeros<2,1>();
		Matrix<2,2> mincov = zeros<2,2>();
		Matrix<2> minc = zeros<2,1>();
		Matrix<2,3> minA = zeros<2,3>();

		std::vector<std::pair<Matrix<2>, double>> cvx; //this is for storing a constraint in workspace.
		cvx.clear();
		//find the maximum likelihodd collision point. the one that has the minimum times of standard deviation (1 after transformation.)
		for(int link = 1; link <=3; link++){
			for(int point = 0; point < seg; point++){
				double length = (point+1)*seglen; //length on the particular link.
				Matrix<2,3> Joc = zeros<2,3>();
				forward_jocabi(qpos, link, length, Joc);
				Matrix<2> c = forward_k(qpos, link, length) - Joc * qpos;
				Matrix<2,3> A = Joc; //A*q + c
				Matrix<2> xmean = forward_k(qpos, link, length);
				Matrix<2,2> xcov = A * qR * ~A; 
				//compute this sample's "distance" to the nearest obstacle.
				bool fflag = true;
				double tmpdis = computeWSConfidence1(xmean, xcov, cal_cloneobstacles, cal_point, fflag); 
				if(fflag == false){ //this means there is no cloned obstacles. all cloned obstacles are pruned away. here we can stop compute constraints.
					flag = false;
					break;
				}
				if(tmpdis < mindis){
					mindis = tmpdis;
					mlink = link;
					mpoint = point;
					minmean = xmean;
					mincov = xcov;
					minc = c;
					minA = A;
				}
			}
			if(flag == false)
				break;
		}
		if(flag == false)
			break;
		//this computes the probablity of this mostlikely-in-collision sample being collision free. and return the workspace constraint. 
		double tmpps = query(minmean, mincov, cal_cloneobstacles, cal_point, cvx);
		ps += (1 - tmpps);
		if((int)cvx.size() == 1){
			Matrix<2> a = cvx[0].first; double b = cvx[0].second; //convert the workspace constraint into a configuration space constraint.
			Matrix<3> ca = ~minA*a; double cb = b - tr(~a*minc);
			ccvx.push_back(std::make_pair(ca, cb));
		}
	}

	return 1 - ps;
}

//compute the probablity of success and also do the truncations.
double LQGMP::computeLQGMPTruncation(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	double newuniMean = 0, newuniVar = 0;
	P0 = initialCov;
	createABVLK();
	int l = (int)pathlqg.size(); 
	std::cout<<l<<std::endl;
	pathlqg[0].y.reset();
	pathlqg[0].R.reset();
	pathlqg[0].R.insert(0,0, P0);

	std::vector<int> cal_cloneobs;
	cal_cloneobs.resize(l-1);

	double ps = 1.0;
	for(int i = 1; i < l; i++){
		//std::cout<<i<<std::endl;
		//clone the obstacles for later truncations. For each step we need to clone the obstacles. only do truncation in the cloned obstacles.
		CAL_CreateGroup(&cal_cloneobs[i-1], 0, true);
		CAL_CloneGroup(&cal_cloneobs[i-1], cal_obstacles, 0, true);
		CAL_SetGroupVisibility(cal_cloneobs[i-1], 0, false);

		Matrix<6,1> augqpos = zeros<6,1>(); augqpos.insert<3,1>(0,0,pathlqg[i-1].T); augqpos.insert<3,1>(3,0, pathlqg[i-1].T);
		augqpos += pathlqg[i-1].y;
		Matrix<6,6> augqcov = pathlqg[i-1].R;
		std::vector<std::pair<Matrix<3>, double>> ccvx;
		ccvx.clear();   //return the constraints in the configuration space. ca * q < cb
		double pps = LQGMPTruncation(augqpos.subMatrix<3,1>(0,0), augqcov.subMatrix<3,3>(0,0), cal_cloneobs[i-1], cal_point, ccvx);
		std::cout<<pps<<std::endl;
		ps *= pps;
		//do truncation against constraints.
		for(int c = 0; c < (int)ccvx.size(); c++){
			Matrix<3> ca = ccvx[c].first;
			double cb = ccvx[c].second;
			Matrix<6> aa = zeros<6,1>();  aa.insert<3,1>(0,0, ca);
			double bb = cb;
			double uniMean = tr(~aa * augqpos);
			double uniVar = tr(~aa * augqcov * aa) + 0.000001; //smooth it.
			truncate(bb, uniMean, uniVar, newuniMean, newuniVar);
			Matrix<6,1> xyCov = augqcov * aa;
			Matrix<6> L = xyCov / uniVar;

			pathlqg[i-1].y -= L*(uniMean - newuniMean);
			pathlqg[i-1].R -= L*(uniVar - newuniVar)*(~L);
		}

		pathlqg[i].y = pathlqg[i-1].F * pathlqg[i-1].y;
		pathlqg[i].R = pathlqg[i-1].F * pathlqg[i-1].R * ~pathlqg[i-1].F + pathlqg[i-1].G * Q * ~pathlqg[i-1].G;

		CAL_DestroyGroup(cal_cloneobs[i-1]);
	}
	return ps;
}



void LQGMP::queryalldistance(const Matrix<3>& qpos, const Matrix<3,3>& qR, const int& cal_obstacles, const int& cal_environment, const int& cal_point, std::vector<informationNode>& Distance)
{
	double seg = 5;
	double seglen = len / seg;
	Distance.clear();
	for(int link = 1; link <=3; link++){
		for(int point = 0; point < seg; point++){
			
			double k = link;
			double length = (point+1)*seglen; //length on the particular link.
			Matrix<2,3> Joc = zeros<2,3>();
			forward_jocabi(qpos, link, length, Joc);
			Matrix<2> c = forward_k(qpos, link, length) - Joc * qpos;
			Matrix<2,3> A = Joc; //A*q + c
			Matrix<2> xmean = forward_k(qpos, link, length);
			Matrix<2,2> xcov = A * qR * ~A; 
			Matrix<2,2> EVec, EVal;
			jacobi(xcov, EVec, EVal);
			Matrix<3,3> Temp = identity<3>();
			Temp.insert(0,0, ~EVec);
			Matrix<4,1> q = quatFromRot(Temp);
			for(int i = 0; i < 2; i++){
				EVal(i,i) += 0.0001;
			}
			//transform environment.
			CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
			CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

			Matrix<2,2> invScale = zeros<2,2>();
			invScale(0,0) = 1/(float)sqrt(EVal(0,0));
			invScale(1,1) = 1/(float)sqrt(EVal(1,1));
			Matrix<2,2> T = invScale * ~EVec;
			Matrix<2,2> invT = !T;
			Matrix<2> transPos =  invScale * ~EVec * xmean;
			CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);
			int num_pairs = 0;
			CAL_GetClosestPairs(cal_point, cal_environment, &num_pairs);
			SCALResult* results = new SCALResult[num_pairs];
			CAL_GetResults(results);
			for(int i = 0; i < num_pairs; i++){
				Matrix<2> objp = zeros<2,1>(); 
				objp[0] = results[i].vector1[0]; objp[1] = results[i].vector1[1];
				int objid = results[i].objID1;
				objp = invT * objp;
				informationNode tmpnode;
				tmpnode.k = k;
				tmpnode.length = length;
				tmpnode.objid = objid;
				tmpnode.objpoint = objp;
				tmpnode.distance = results[i].distance;
				tmpnode.A = A;
				tmpnode.x = xmean;
				tmpnode.c = c;
				Distance.push_back(tmpnode);
			}
			delete []results;

		}
	}
}

int LQGMP::pick_min_dis(std::vector<informationNode>& Distance){
	int minnum = -1;
	double mindis = 999999999999.0;
	
	for(int i = 0; i < (int)Distance.size(); i++){
		double tmpdis = Distance[i].distance;
		if(tmpdis < mindis){
			mindis = tmpdis;
			minnum = i;
		}
	}

	return minnum;
}

bool LQGMP::check_contact_point(const std::vector<std::pair<Matrix<2>, double>>& cvx, const Matrix<2>& objpoint)
{
	for(int i = 0; i < (int)cvx.size(); i++){
		Matrix<2> a = cvx[i].first;
		double b = cvx[i].second;
		if(tr(~a*objpoint) > b) //means it is already outside of the current convex region.
			return false;
	}
	return true; //
}


void LQGMP::DistanceToConstraints(const Matrix<3>& qpos, const Matrix<3,3>& qR, std::vector<std::pair<Matrix<3>, double>>& qcvx, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	std::vector<informationNode> Distances;
	Distances.clear();
	queryalldistance(qpos, qR, cal_obstacles, cal_environment, cal_point, Distances);

	int numdis = (int)Distances.size();
	std::cout<<numdis<<std::endl;

	std::vector<std::pair<Matrix<2>, double>> cvx;
	cvx.clear();
	int obj_left = 3;
	while(obj_left > 0){
		int num = pick_min_dis(Distances);
		if(num == -1)
			return;

		informationNode minnode = Distances[num];
		bool in = check_contact_point(cvx, minnode.objpoint);
		if(in == true){ 
			//if it is stll in the current convex region, we need construnct a new half plane for the current convex region.
			Matrix<2> op = minnode.x;
			Matrix<2> objp = minnode.objpoint;
			Matrix<2> a = objp - op;
			a = a / sqrt(tr(~a*a));
			double b = tr(~a * objp);
			cvx.push_back(std::make_pair(a, b));
			Matrix<3> qa = ~(minnode.A)*a;
			double qb = b - tr(~a*(minnode.c));
			qcvx.push_back(std::make_pair(qa, qb));
		}
		Distances.erase(Distances.begin() + num);
	}
	cvx.clear();
	Distances.clear();
}

double LQGMP::EstimateAndTruncate(const Matrix<3,3> initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	double ps = 1.0;
	P0 = initialCov;
	double newuniMean = 0, newuniVar = 0;
	P0 = initialCov;
	createABVLK();
	int l = (int)pathlqg.size(); 
	std::cout<<l<<std::endl;
	pathlqg[0].y.reset();
	pathlqg[0].R.reset();
	pathlqg[0].R.insert(0,0, P0);
	
	for(int i = 0; i < l; i++){
		std::vector<std::pair<Matrix<3>, double>> qcvx;
		qcvx.clear();
		Matrix<6> augpos = zeros<6,1>(); 
		augpos.insert<3,1>(0,0, pathlqg[i].T); augpos.insert<3,1>(3,0, pathlqg[i].T);
		augpos += pathlqg[i].y;
		Matrix<6,6> augcov = zeros<6,6>();
		augcov = pathlqg[i].R;
		Matrix<3,1> qpos = augpos.subMatrix<3,1>(0,0);
		Matrix<3,3> qcov = augpos.subMatrix<3,3>(0,0);
		DistanceToConstraints(augpos.subMatrix<3,1>(0,0), augcov.subMatrix<3,3>(0,0), qcvx, cal_obstacles, cal_environment, cal_point);

		//compute the probablity of collision for this step.
		double ps_c = 0;
		int lencvx = (int)qcvx.size();
		for(int c=0; c < lencvx; c++){
			Matrix<3,1> a = qcvx[c].first;
			double b = qcvx[c].second;
			double alpha = (b - tr(~a*qpos)) / sqrt(tr(~a*qcov*a));
			ps_c += (1 - cdf(alpha));
		}
		ps *= (1 - ps_c);

		//truncations;
		for(int c = 0; c < lencvx; c++){
			Matrix<6,1> aa = zeros<6,1>(); aa.insert<3,1>(0,0, qcvx[c].first);
			double bb = qcvx[c].second;
			double oldMean = tr(~aa * augpos);
			double oldVar = tr(~aa*augcov*aa);
			truncate(bb, oldMean, oldVar, newuniMean, newuniVar);
			Matrix<6,1> xyCov = augcov * aa;
			Matrix<6> L = xyCov / oldVar;
			pathlqg[i].y -= L*(oldMean - newuniMean);
			pathlqg[i].R -= L*(oldVar - newuniVar) * (~L);
		}

		if( i < l - 1){
			pathlqg[i+1].y = pathlqg[i].F * pathlqg[i].y;
			pathlqg[i+1].R = pathlqg[i].F * pathlqg[i].R * (~pathlqg[i].F) + pathlqg[i].G * Q * (~pathlqg[i].G);
		}
	}
	return ps;
}




