// Copyright
#include "RungeKutta.hpp"

RungeKutta::RungeKutta() {
  fxpoints = fopen("fxpoints.txt", "w");
  fxpoints2 = fopen("fxpoints2.txt", "w");
}

RungeKutta::~RungeKutta() {}

float RungeKutta::f(float x, float y) {
  float fx;
  // fx = g - (c/m)*y;
  fx = -2*pow(x, 3) + 12*pow(x, 2) - 20*x + 8.5;
  return fx;
}

float RungeKutta::real(float x) {
  float y;
  // y = (g*m)/c * (1 - exp((-c/m)*x));
  y = -0.5*pow(x, 4) + 4*pow(x, 3) - 10*pow(x,2) + 8.5*x +1;
  return y;
}

float RungeKutta::RKfourth(float x, float y) {
  float yk, k1, k2, k3, k4;
  k1 = f(x, y);
  k2 = f(x + 0.5*h, y + 0.5*k1*h);
  k3 = f(x + 0.5*h, y + 0.5*k2*h);
  k4 = f(x + h, y + k3*h);
  yk = y + ((k1 + 2*k2 + 2*k3 + k4)/6)*h;
  return yk;
}

void RungeKutta::getValues(float step, int n) {
  float x, y, yk, yreal, error;
  h = step;
  x = xi;
  y = yi;
  error = 0;
  std::cout << "Fourth-order Runge-Kutta's method with step size: " << h << std::endl;
  printHeader();
  printElement(x, y, y, error);
  savePoints(x, y, n);
  while (x < xf) {
    yk = RKfourth(x, y);
    x += h;
    y = yk;
    yreal = real(x);
    error = 100*fabs(yreal-y)/yreal;
    printElement(x, yreal, y, error);
    savePoints(x, y, n);
  }
  std::cout << std::endl;
  if (n == 1)
    plot();
  else
    replot();
}

void RungeKutta::RKMethod(float step1, float step2) {
  getValues(step1, 1);
  getValues(step2, 2);
}

void RungeKutta::printElement(float x, float yr, float y, float e) {
  const int precision = 5;
  const int& width = 10;
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << x << std::fixed << std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << yr << std::fixed << std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << y << std::fixed <<  std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << e << std::fixed <<  std::setw(width);
  std::cout << '|' << std::endl;
}

void RungeKutta::printHeader() {
  const int& width = 10;
  std::cout << '|' << std::setw(width) << "x" << std::setw(width);
  std::cout << '|' << std::setw(width) << "yReal" << std::setw(width);
  std::cout << '|' << std::setw(width) << "yRK" << std::setw(width);
  std::cout << '|' << std::setw(width) << "error" << std::setw(width);
  std::cout << '|' << std::endl;
}

void RungeKutta::savePoints(float x, float y, int n) {
  if (n == 1) {
    fprintf(fxpoints, "%lf %lf \n", x, y);
  }
  else {
    fprintf(fxpoints2, "%lf %lf \n", x, y);
  }
}

void RungeKutta::plot() {
  gnuplot = popen("gnuplot -persist", "w");
  fprintf(gnuplot, "set title \"Fourth-order Runge-Kutta for y' = -2x^3 + 12x^2 -20x + 8.5\" font \"Times-Roman,14\"\n");
  fprintf(gnuplot, "f(x) = %s\n", "-0.5*x**4 + 4*x**3 - 10 *x**2 + 8.5*x +1");
  // fprintf(gnuplot, "set title \"Fourth-order Runge-Kutta for y' = g - c/m *y\" font \"Times-Roman,14\"\n");
  // fprintf(gnuplot, "f(x) = (%f*%f)/%f * (1 - exp((-%f/%f)*x))\n", g, m, c, c, m);
  fprintf(gnuplot, "set xrange[%d:%d]\n", xi, xf);
  fprintf(gnuplot, "plot f(x) title \"real solution\", \"fxpoints.txt\" u 1:2 w l lc rgb '#0072bd' title \"h = %.1f\", '' u 1:2 w p ls 7 ps 2 lc rgb '#0072bd' notitle\n", h);
  fclose(fxpoints);
}

void RungeKutta::replot() {
  fprintf(gnuplot, "replot \"fxpoints2.txt\" u 1:2 w l lc rgb '#DC143C' title \"h = %.1f\", '' u 1:2 w p ls 7 ps 1.5 lc rgb '#DC143C' notitle\n", h);
  fclose(fxpoints2);
}

int main() {
  RungeKutta rungeKutta;
  rungeKutta.RKMethod(0.2, 1);
  return 0;
}
