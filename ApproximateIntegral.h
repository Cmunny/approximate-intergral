#ifndef INTEGRAL_H
#define INTEGRAL_H


#include <string>
#include <vector>

using namespace std;

class ApprIntegral {
  static const string trig[];
  static vector<string> SeperateTermByAdd(string);
  static vector<string> SeperateTermByMulti(string);
  static double CalcPolynomialAndExp(string , double);
  static double CalcTrig(string, double);
  static double CalcTerm(string, double);
  static string SubArgument(string);
  static double CalcLog(string, double);
  ApprIntegral(){}
public:
   static double CalcEquation(string equation, double num);
   static double ApproximateIntegral(string equation, double a, double b, int numOfRects);
};


#endif
