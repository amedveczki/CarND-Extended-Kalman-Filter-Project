#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::Vector4d;
using Eigen::MatrixXd;
using Eigen::Matrix4d;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	assert(!estimations.empty() && estimations.size() == ground_truth.size());

	for (int i=0; i < estimations.size(); ++i) {
		VectorXd x = (estimations[i] - ground_truth[i]);
		x = x.array() * x.array();
		rmse += x/estimations.size();

	}

	return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);

	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float p2 = px*px + py*py;
	float sp2 = sqrt(p2);
	float sp32 = sp2*sp2*sp2;

	if (p2 == 0)
	{
        std::cout << "Division by zero" << std::endl;
		return Hj;
	}
	Hj << px/sp2, py/sp2, 0, 0,
	   -py/p2, px/p2, 0, 0,
	   py*(vx*py-vy*px)/sp32, px*(vy*px-vx*py)/sp32, px/sp2, py/sp2;

	return Hj;
}

