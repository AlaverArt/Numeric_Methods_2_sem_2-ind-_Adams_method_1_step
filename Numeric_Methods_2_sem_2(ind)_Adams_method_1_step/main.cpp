#include<iostream>
#include<functional>
#include<vector>
#include<fstream>
//ex1
//constexpr double u_0 = 0.062233835;
//constexpr double t_0 = 0.5;
//constexpr double t_end = 2.0;

//ex2
constexpr double u_0 = 0.062233835;
constexpr double t_0 = 0.5;
constexpr double t_end = 2.0;

constexpr double tau = 0.05;
constexpr double b0 = 0.5;
constexpr double b1 = 0.5;

constexpr double epsSimpleIter = 0.0000001;
constexpr double epsNewton = 0.0000001;

double f(double tn, double yn) {
	//ex1
	//return yn * (3 + 2 * tn) / (tn * tn);//y'=(3y+2ty)/(t^2)

	//ex2
	return -(yn + 2 * tn * yn * yn * log(tn)) / tn;//y'=-(y+2ty*lnt)/t
}

double f_newton(double yn, double tn, double yn_sub_1, double fn_sub_1) {
	//yn - variable
	//tn, yn_sub_1, fn_sub_1 - constant
	return yn - tau * b0 * f(tn, yn) - yn_sub_1 - tau * b1 * fn_sub_1;
}

double diff_f_newton(std::function<double(double, double, double, double)> func,
	double y, double yn_sub_1, double fn_sub_1, double tn) {

	constexpr double dy = 0.0000001;
	return (func(y + dy, tn, yn_sub_1, fn_sub_1) - func(y, tn, yn_sub_1, fn_sub_1)) / dy;
}

double newtonesMethod(std::function<double(double, double, double, double)> func,
						double yn0, double tn, double yn_sub_1, double fn_sub_1){

	//x0 - initial approximation(начальное приближение)

	double yk = yn0;
	double yk_prev;
	do {
		yk_prev = yk;

		yk = yk_prev - func(yk_prev, tn, yn_sub_1, fn_sub_1) / diff_f_newton(func, yk_prev, yn_sub_1, fn_sub_1, tn);

	} while (abs(yk - yk_prev) > epsNewton);
	return yk;
}

std::vector<double> adamsMethod_newton(std::function<double(double, double)> f, double y0, double t0, double end_t) {
	int numb_yn = (int)(abs(end_t - t0) / tau) + 1;
	std::vector<double> y(numb_yn);
	double fn_sub_1;
	y[0] = y0;
	fn_sub_1 = f(t0, y[0]);

	double tn = t0;
	for (int n = 1; n < numb_yn; n++)
	{
		tn += tau;
		y[n] = newtonesMethod(f_newton, y[n - 1], tn, y[n - 1], fn_sub_1);
		fn_sub_1 = f(tn, y[n]);
	}

	return y;
}

std::vector<double> adamsMethod_simpleIters(std::function<double(double, double)> f, double y0, double t0, double end_t) {
	int numb_yn = (int)(abs(end_t - t0) / tau) + 1;
	std::vector<double> y(numb_yn);
	y[0] = y0;

	double y_pred;
	double fn_sub_1;
	double fn = f(t0, y0);
	double tn = t0;
	for (int n = 1; n < numb_yn; n++) {
		tn += tau;

		//computing yn
		fn_sub_1 = f(tn - tau, y[n - 1]);
		y[n] = tau * b1 * fn_sub_1 + y[n - 1];

		do {
			y_pred = y[n];
			fn_sub_1 = fn;
			fn = f(tn, y_pred);
			y[n] = tau * (b0 * fn + b1 * fn_sub_1) + y[n - 1];

		} while (abs(y[n] - y_pred) > epsSimpleIter);
	}

	return y;
}

int main(int argc, char* argv)
{
	std::ofstream fout_newton("output_newton.txt");
	std::ofstream fout_simpleIters("output_simpleIters.txt");

	std::vector<double> u_newton = adamsMethod_newton(f, u_0, t_0, t_end);
	for (double item : u_newton) {
		//std::cout << item << std::endl;
		fout_newton << item << "\n";
	}

	std::vector<double> u_simpleIters = adamsMethod_simpleIters(f, u_0, t_0, t_end);
	for (double item : u_simpleIters) {
		//std::cout << item << std::endl;
		fout_simpleIters << item << "\n";
	}

	return 0;
}