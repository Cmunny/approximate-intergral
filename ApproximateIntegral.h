#ifndef INTEGRAL_H
#define INTEGRAL_H


#include <string>
#include <vector>

using namespace std;

class ApprIntegral {
  static const string trig[];
  static vector<string> TermSeperater(string );
  static double CalcPolynomial(string , double );
  static double CalcTrig(string, double , const vector<int>& );
  static double termCalc(string, double );
  static string SubArgument(string);
  ApprIntegral() {};
public:
   static double calcEquation(string equation, double num);
   static double approximateIntegral(string equation, double a, double b, int numOfRects);
};


#endif
