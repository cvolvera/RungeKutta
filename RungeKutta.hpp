#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#define xi 0
#define xf 4
#define yi 1
#define g 9.8
#define c 12.5
#define m 68.1

class RungeKutta {
 private:
  float h;
  FILE *fxpoints;
  FILE *fxpoints2;
  FILE *gnuplot;

 public:
  RungeKutta();
  ~RungeKutta();

  float f(float x, float y);
  float real(float x);
  float RKfourth(float x, float y);
  void getValues(float step, int n);
  void RKMethod(float step1, float step2);

  void printElement(float x, float yr, float y, float e);
  void printHeader();

  void savePoints(float x, float y, int n);
  void plot();
  void replot();
};
