#include "RRT.h"

void RRT::setPlannerDefaults()
{
	factor = 2.0;
	K = 0.35 * 3.1415926;
	//w_max = k / ( v + 1).

	planbias = 0.1;
	maxiter = 1000;
	maxNNiter = 1000;

	maxtries = 3000;

	maxstep = 1;
	maxchildren = 150;
	maxattempts = 100;
}


void RRT::initTree(std::vector<TreeNode>& tree, const Matrix<3>& pose)
{
	TreeNode n;
	n.T = pose;

	n.bp = -1;

	n.depth = 0.0;
	n.marked = false;
	n.attempts = 0;

	tree.push_back(n);
}

double RRT::dist(const Matrix<3>& p1, const Matrix<3>& p2)
{
	double dis = 0;
	dis = sqrt(tr(~(p1 - p2) * (p1 - p2)));

	return dis;
}

int RRT::nearestNeighbor(const std::vector<TreeNode>& tree, const Matrix<3>& point)
{
	int closest = -1;
	double mindist = 10000;

	for(int i = 0; i < tree.size(); i++){
		Matrix<3> tmpnode = tree[i].T;
		double tmpdis = dist(tmpnode, point);
		if(tmpdis < mindist && tree[i].attempts < maxattempts && tree[i].marked == false){
			closest = i;
			mindist = tmpdis;
		}
	}

	return closest;
}


Matrix<3> RRT::generateControl(const Matrix<3>& p1, const Matrix<3>& p2)
{
	int maxi = 0;
	double val = -1;
	for(int i = 0; i < 3; i++){
		if(abs(p1[i] - p2[i]) > val){
			maxi = i;
			val = abs(p1[i] - p2[i]);
		}
	}

	Matrix<3> w = zeros<3,1>();
	w[maxi] = (p1[maxi] - p2[maxi]) / (abs(p1[maxi] - p2[maxi])) * w_max;
	for(int i = 0; i < 3; i++){
		if( i != maxi){
			w[i] = w[maxi] * (p1[i] - p2[i]) / (p1[maxi] - p2[maxi]);
		}
	}

	double alpha = random();
	w = w * alpha;
	return w;
}


bool RRT::rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths, bool& found, const clock_t& sc){
	
	found = false;
	if((double)((clock() - sc)/CLOCKS_PER_SEC) > plantime)
		return false;

	int node = -1;
	int tries = 0;

	Matrix<3> point;
	point.reset();
	double rgoal = rthreshold;

	tries = 0;
	do{

		if(random() < 0.05){
			point[0] = 7.0 / 4.0 * 3.1416 - 0.04 + random() * 0.08;
			point[1] = - 0.04 + random() * 0.08;
			point[2] = - 0.04 + random() * 0.08;
		} 
		else{
			point[0] = 0 + random() * 2*3.1416;
			point[1] = 0 + random() * 2*3.1416;
			point[2] = 0 + random() * 2*3.1416;
		}

		bool success = pointCollision(point, cal_link1, cal_link2, cal_link3, cal_obstacles);

		if(success == true){
			node = nearestNeighbor(tree, point);
			if(node != -1){
				tree[node].attempts ++;
			}
		}

	} while ((node == -1) && (++tries < maxNNiter));
	
	if(tries == maxtries){
			return false;  //return false when cannot sample
	}

	if(node == -1)
		return false;

	Matrix<3> x_new, x_old;
	x_new.reset(); x_old.reset();
	Matrix<3> rancontrol;
	rancontrol.reset();
	x_old = tree[node].T;

	rancontrol = generateControl(point, x_old);
	//rancontrol[0] = w_min + random() * (w_max - w_min);
	//rancontrol[1] = w_min + random() * (w_max - w_min);
	//rancontrol[2] = w_min + random() * (w_max - w_min);

	bool valid = false;
	bool success = true;

	Dynamics dyn(dt);
	x_new = dyn.dynamics_zero_noise(x_old, rancontrol);
		
	success = lineCollision(x_old, x_new, cal_link1, cal_link2, cal_link3, cal_obstacles);
	//x_old = x_new;
	
	if(success == true){

		//float line[6] = {x_old[0], x_old[1], x_old[2], x_new[0], x_new[1], x_new[2]};
		//int np[1] = {2};
		//CAL_CreatePolyline(cal_rrt, 1, np, line);

		TreeNode newnode;
		newnode.T = x_new;
		newnode.bp = node;
		newnode.u = rancontrol;
		newnode.marked = false;
		newnode.attempts = 0;
		int newid = (int) tree.size();
		tree[node].children.push_back(newid);
		tree.push_back(newnode);
	
		Matrix<2,1> endp = forward_k(x_new, 3,len);
		//std::cout<<endp[0]<<" "<<endp[1]<<std::endl;
		double dg = tr(~(endp - goal) * (endp - goal));
		if(dg < rthreshold * rthreshold)
		{
			TreeNode& tnode = tree[newid];
			tnode.marked = true;
			paths.push_back(newid);
			found = true;
			return true;
		}
	}
	return true;
}

bool RRT::buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, const clock_t& sc)
{
	bool stepRet = true;
	// Initialize tree
	initTree(tree, start);

	bool found = false;
	for (int i = 0; stepRet ; ++i) {
		stepRet = rrtStep(tree, paths, found, sc);
		if(found == true)
			break;
	}
	if (stepRet && !paths.empty()) 
	{
		//std::cout << "\nRRT build: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
		int numpaths = (int)paths.size();
		//std::cout << "Num paths found: " << numpaths << std::endl;

	} 
	else {
		std::cout << "Unable to find solution, reusing previous solution" << std::endl;
	}

	return (stepRet && !paths.empty());
}

bool RRT::executeRRT(const clock_t& sc)
{
	rrttree.clear();
	rrtpath.clear();
	rrtpaths.clear();

	return buildRRT(rrttree, rrtpaths, sc);
}

bool RRT::Plan_K_Seconds()
{
	clock_t startCLK = clock();
	
	while(((double)(clock() - startCLK) / CLOCKS_PER_SEC) < plantime){
	
		bool s = executeRRT(startCLK);
		if(s == true){
			
			std::vector<PathNode> onePath;
			onePath.clear();
			int tmp = rrtpaths[0];
			PathNode tmpnode;
			tmpnode.T = rrttree[tmp].T;
			tmpnode.u = rrttree[tmp].u;
			onePath.push_back(tmpnode);

			int parent = -1;
			while(parent != 0){
				parent = rrttree[tmp].bp;
				tmpnode.T = rrttree[parent].T;
				tmpnode.u = rrttree[parent].u;
				onePath.push_back(tmpnode);
				tmp = parent;
			}
			std::reverse(onePath.begin(), onePath.end());

			for(int i = 0; i < (int) onePath.size() - 1; i++){
				onePath[i].u = onePath[i+1].u;
			}
			onePath[(int)onePath.size() - 1].u[0] = 0;
			onePath[(int)onePath.size() - 1].u[1] = 0;

			pathSet.push_back(onePath);
			onePath.clear();
		}
	}

	if(pathSet.size() == 0){
		std::cout<<"No Path find in "<<plantime<<" seconds"<<std::endl;
		return false;
	}
	//std::cout<<pathSet.size()<<std::endl;
	return true;
}