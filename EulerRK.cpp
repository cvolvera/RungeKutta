// Copyright
#include "EulerRK.hpp"

EulerRK::EulerRK() {
  fxpoints = fopen("fxpoints.txt", "w");
  fxpoints2 = fopen("fxpoints2.txt", "w");
}

EulerRK::~EulerRK() {}

float EulerRK::f(float x, float y) {
  float fx;
  // fx = g - (c/m)*y;
  fx = -2*pow(x, 3) + 12*pow(x, 2) - 20*x + 8.5;
  return fx;
}

float EulerRK::real(float x) {
  float y;
  // y = (g*m)/c * (1 - exp((-c/m)*x));
  y = -0.5*pow(x, 4) + 4*pow(x, 3) - 10*pow(x,2) + 8.5*x +1;
  return y;
}

float EulerRK::Eulerformula(float x, float y) {
  float yk;
  yk = y + f(x, y)*h;
  return yk;
}

float EulerRK::RKfourth(float x, float y) {
  float yk, k1, k2, k3, k4;
  k1 = f(x, y);
  k2 = f(x + 0.5*h, y + 0.5*k1*h);
  k3 = f(x + 0.5*h, y + 0.5*k2*h);
  k4 = f(x + h, y + k3*h);
  yk = y + ((k1 + 2*k2 + 2*k3 + k4)/6)*h;
  return yk;
}

void EulerRK::getValues(float step, int n) {
  float x, y, yk, yRK, yreal, error, errorRK;
  h = step;
  x = xi;
  y = yi;
  yRK = yi;
  error = 0;
  std::cout << "First and Fourth order Runge-Kutta methods with step size: " << h << std::endl;
  printHeader();
  printElement(x, y, y, y, error, error);
  savePoints(x, y, yRK, n);
  while (x < xf) {
    yk = Eulerformula(x, y);
    y = yk;
    yk = RKfourth(x, yRK);
    yRK = yk;
    x += h;
    yreal = real(x);
    error = 100*fabs(yreal-y)/yreal;
    errorRK = 100*fabs(yreal-yRK)/yreal;
    printElement(x, yreal, y, yRK, error, errorRK);
    savePoints(x, y, yRK, n);
  }
  std::cout << std::endl;
  if (n == 1)
    plot();
  else
    replot();

}

void EulerRK::EulerRKmethods(float step1, float step2) {
  getValues(step1, 1);
  getValues(step2, 2);
}

void EulerRK::printElement(float x, float yr, float y, float yRK, float e, float eRK) {
  const int precision = 5;
  const int& width = 10;
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << x << std::fixed << std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << yr << std::fixed << std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << y << std::fixed <<  std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << yRK << std::fixed <<  std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << e << std::fixed <<  std::setw(width);
  std::cout << '|' << std::setw(width) << std::setprecision(precision) << eRK << std::fixed <<  std::setw(width);
  std::cout << '|' << std::endl;
}

void EulerRK::printHeader() {
  const int& width = 10;
  std::cout << '|' << std::setw(width) << "x" << std::setw(width);
  std::cout << '|' << std::setw(width) << "yReal" << std::setw(width);
  std::cout << '|' << std::setw(width) << "yEuler" << std::setw(width);
  std::cout << '|' << std::setw(width) << "yRK" << std::setw(width);
  std::cout << '|' << std::setw(width) << "e Euler" << std::setw(width);
  std::cout << '|' << std::setw(width) << "e RK" << std::setw(width);
  std::cout << '|' << std::endl;
}

void EulerRK::savePoints(float x, float y, float yRK, int n) {
  if (n == 1) {
    fprintf(fxpoints, "%lf %lf %lf \n", x, y, yRK);
  }
  else {
    fprintf(fxpoints2, "%lf %lf %lf \n", x, y, yRK);
  }
}

void EulerRK::plot() {
  gnuplot = popen("gnuplot -persist", "w");
  fprintf(gnuplot, "set title \"First and Fourth order RK Methods for y' = -2x^3 + 12x^2 -20x + 8.5\" font \"Times-Roman,14\"\n");
  fprintf(gnuplot, "f(x) = %s\n", "-0.5*x**4 + 4*x**3 - 10 *x**2 + 8.5*x +1");
  // fprintf(gnuplot, "set title \"First and Fourth order RK Methods for y' = g - c/m *y\" font \"Times-Roman,14\"\n");
  // fprintf(gnuplot, "f(x) = (%f*%f)/%f * (1 - exp((-%f/%f)*x))\n", g, m, c, c, m);
  fprintf(gnuplot, "set xrange[%d:%d]\n", xi, xf);
  fprintf(gnuplot, "plot f(x) title \"real solution\", \"fxpoints.txt\" u 1:2 w l lc rgb '#DC143C' title \"Euler h = %.1f\", '' u 1:2 w p ls 7 ps 2 lc rgb '#DC143C' notitle, '' u 1:3 w l lc rgb '#0072bd' title \"RK h = %.1f\", '' u 1:3 w p ls 7 ps 2 lc rgb '#0072bd' notitle\n", h, h);
  fclose(fxpoints);
}

void EulerRK::replot() {
  fprintf(gnuplot, "replot \"fxpoints2.txt\" u 1:2 w l lc rgb '#edb120' title \"Euler h = %.1f\", '' u 1:2 w p ls 7 ps 2 lc rgb '#edb120' notitle, '' u 1:3 w l lc rgb '#77ac30' title \"RK h = %.1f\", '' u 1:3 w p ls 7 ps 2 lc rgb '#77ac30' notitle\n", h, h);
  fclose(fxpoints2);
}

int main() {
  EulerRK eulerRK;
  eulerRK.EulerRKmethods(0.2,1);
  return 0;
}
