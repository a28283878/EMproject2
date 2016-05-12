#pragma once
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "vector"
class NewtonMethod {
public:
	NewtonMethod() {

	};
	struct hessian {
		std::vector<std::vector<double>> contains;
	};
	double h = 0.0000001;
	double function(double x, double y);
	double gradientx(double x, double y);
	double gradienty(double x, double y);
	hessian gethessian(double x, double y);
	hessian gethessianinverse(hessian x);
	std::vector<double> d(double x, double y, hessian inverse);
};