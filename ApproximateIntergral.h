#ifndef INTEGRAL_H
#define INTEGRAL_H


#include <string>
#include <vector>
<<<<<<< HEAD
#include <sstream>
#include <cmath>

using namespace std;

class ApprIntegral {
=======

using namespace std;

class ApprIntegral
{
>>>>>>> master
  static const string trig[];
  static vector<string> TermSeperater(string str);
  static double CalcPolynomial(string term, double num);
  static double CalcTrig(string term, double num, int trigFuncIndex);
  static double termCalc(string term, double num);
<<<<<<< HEAD
  ApprIntegral() {};
public:
   static double calcEquation(string equation, double num);
   static double approximateIntegral(string equation, double a, double b, int numOfRects);
=======
  ApprIntegral(){};
public:
  static double calcEquation(string equation, double num);
  static double approximateIntegral(string equation, double a, double b, int numOfRects);
>>>>>>> master
};


#endif