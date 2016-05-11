#ifndef INTEGRAL_H
#define INTEGRAL_H


#include <string>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

class ApprIntegral {
    static const string trig[];
  

  ApprIntegral() {};
public:

  static vector<string> TermSeperater(string str);

   static double CalcPolynomial(string term, double num);

   static double CalcTrig(string term, double num, int trigFuncIndex);

   static double termCalc(string term, double num);



   static double calcEquation(string equation, double num);



   static double approximateIntegral(string equation, double a, double b, int numOfRects);
};


#endif
