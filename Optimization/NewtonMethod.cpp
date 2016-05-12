#include "NewtonMethod.h"
#include "stdlib.h"
#include "stdio.h"
#include "iostream"
using namespace std;
double NewtonMethod::function(double x,double y) {
	return (pow(x, 2) + 3.25*pow(y, 2) - 3 * x*y - 4 * y + 7);
}

double NewtonMethod::gradientx(double x, double y) {
	double k = function(x + h, y);
	double j = function(x - h, y);
	//cout << k << "   " << j << "   " << k - j << "   " << 2 * h << endl;
	return (k - j) / (2 * h);
}

double NewtonMethod::gradienty(double x, double y) {
	double k = function(x, y + h);
	double j = function(x, y - h);
	return (k - j) / (2 * h);
}
NewtonMethod::hessian NewtonMethod::gethessian(double x, double y) {
	NewtonMethod::hessian output;
	for (int i = 0; i < 2; i++) {
		vector<double> temp;
		for (int j = 0; j < 2; j++) {
			if (i == 0 && j == 0) {
				double k = gradientx(x + h, y);
				double j = gradientx(x - h, y);
				temp.push_back((k - j) / (2 * h));
			}
			else if (i == 0 && j == 1) {
				double k = gradientx(x, y + h);
				double j = gradientx(x, y - h);
				temp.push_back((k - j) / (2 * h));
			}
			else if (i == 1 && j == 0) {
				double k = gradienty(x + h, y);
				double j = gradienty(x - h, y);
				temp.push_back((k - j) / (2 * h));
			}
			else {
				double k = gradienty(x, y + h);
				double j = gradienty(x, y - h);
				temp.push_back((k - j) / (2 * h));
			}
		}
		output.contains.push_back(temp);
	}
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << output.contains[i][j] << "   ";
		}
		cout << endl;
	}
	return output;
}

NewtonMethod::hessian NewtonMethod::gethessianinverse(NewtonMethod::hessian x) {
	double det = x.contains[0][0] * x.contains[1][1] - x.contains[0][1] * x.contains[1][0];
	NewtonMethod::hessian output;
	vector<double> temp;
	temp.push_back(x.contains[1][1] / det);
	temp.push_back(-1 * x.contains[0][1] / det);
	output.contains.push_back(temp);
	temp[0] = -1 * x.contains[1][0] / det;
	temp[1] = x.contains[0][0] / det;
	output.contains.push_back(temp);
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			cout << output.contains[i][j] << "   ";
		}
		cout << endl;
	}
	return output;
}

vector<double> NewtonMethod::d(double x, double y, NewtonMethod::hessian inverse) {
	double d[2];//微分一次的值
	d[0] = gradientx(x, y);
	d[1] = gradienty(x, y);
	vector<double> output;
	for (int i = 0; i < 2; i++) {
		double sum = 0;
		for (int j = 0; j < 2; j++) {
			sum += inverse.contains[i][j] * d[j];
		}
		output.push_back(sum);
	}
	return output;
}