#include<iostream>
#include<functional>
#include<vector>
#include<fstream>
constexpr double u_0 = 0.062233835;
constexpr double t_0 = 0.5;
constexpr double t_end = 0.9;

constexpr double tau = 0.05;
constexpr double b0 = 0.5;
constexpr double b1 = 0.5;

double f(double tn, double yn) {
	return (3 * yn + 2 * tn*yn) / (tn * tn);//y'=(3y+2ty)/(t^2)
}

double f_newton(double yn, double tn, double yn_sub_1, double fn_sub_1) {
	//yn - variable
	//tn, yn_sub_1, fn_sub_1 - constant
	return yn - tau * b0 * f(tn, yn) - yn_sub_1 - b1 * fn_sub_1;
}

double diff_f_newton(std::function<double(double, double, double, double)> func,
	double x, double yn_sub_1, double fn_sub_1, double tn) {

	constexpr double dx = 0.0000001;
	return (func(x + dx, tn, yn_sub_1, fn_sub_1) - func(x, tn, yn_sub_1, fn_sub_1)) / dx;
}

double newtonesMethod(std::function<double(double, double, double, double)> func,
						double yn0, double tn, double yn_sub_1, double fn_sub_1){

	//x0 - initial approximation(начальное приближение)
	constexpr double newtonsEps = 0.00001;//accuracy(погрешность)

	double yk = yn0;
	double yk_prev;
	do {
		yk_prev = yk;

		yk = yk_prev - func(yk_prev, tn, yn_sub_1, fn_sub_1) / diff_f_newton(func, yk_prev, yn_sub_1, fn_sub_1, tn);
	} while (abs(yk - yk_prev) > newtonsEps);
	return yk;
}

std::vector<double> adamsMethod(std::function<double(double, double)> f, double y0 ,double t0, double end_t) {
	int numb_yn = (int)((end_t - t0) / tau + 1);
	std::vector<double> y(numb_yn);
	y[0] = y0;

	double fn_sub_1 = f(t0, y0);
	double tn = t0;

	//yn=tau(b0*fn+b1*fn_1) + yn_1
	for (int n = 1; n < numb_yn; n++) {
		//computing yn

		tn += tau;
		y[n] = newtonesMethod(f_newton, y[n - 1], tn, y[n - 1], fn_sub_1);

		fn_sub_1 = f(tn, y[n]);
	}
	return y;
}

int main(int argc, char* argv) {
	std::ofstream fout("output.txt");

	std::vector<double> u = adamsMethod(f, u_0, t_0, t_end);
	for (double item : u) {
		std::cout << item << " ";
		fout << item << "\n";
	}

	return 0;
}