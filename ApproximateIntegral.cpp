#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include "ApproximateIntegral.h"
using namespace std;

const string ApprIntegral::trig[] = {"sin","cos","tan","csc","sec","cot", "arccos" "arcsin", "arctan", "cosh", "sinh", "tanh", "arccosh", "arcsinh", "arctanh",};
const double ApprIntegral::PI = 3.14159265359;
vector<string> ApprIntegral::SeperateTermByAdd(string str)
{
  vector<string> terms;
  size_t found = str.find_first_of("+-");
  size_t startPos = 0;
  size_t firstLeftParenth = str.find('(');
  size_t lastRightParenth = str.rfind(')');

  //ignores the the sign if it is the very beginning of term. i.e. -cosx + 1
  if (found == 0)
    found = str.find_first_of("+-", found + 1);

  //seperates the equation into seperate terms and stores them the terms vector
  //the leading operator is stored with term. i.e -3x or +x;
  while (found != string::npos)
  {
    //skips the operators in between parentheses.
    //i.e x + cos(1+x) gets seperated into x and +cos(1+x).
    if (lastRightParenth != string::npos && firstLeftParenth != string::npos && found > firstLeftParenth && found < lastRightParenth)
      found = str.find_first_of("+-", found + 1);

    if (found != string::npos)
    {
      string s = str.substr(startPos, found - startPos);
      terms.push_back(s);
      startPos = found;
      found = str.find_first_of("+-", found + 1);
    }
  }
  terms.push_back(str.substr(startPos, string::npos));
  return terms;
}

vector<string> ApprIntegral::SeperateTermByMulti(string str)
{
  vector<string> terms;
  size_t found = str.find_first_of("*/");
  size_t startPos = 0;

  //seperates the equation into seperate terms and stores them the terms vector
  //the leading operator is stored with term. i.e -3x or +x;
  while (found != string::npos)
  {
    //skips the operators in between parentheses.
    //i.e x + cos(1+x) gets seperated into x and +cos(1+x).
    size_t currentTermLeftParenth = str.find("(", startPos);
    size_t currentTermRightParenth = str.rfind(")", found);
    vector<size_t> leftParenths;
    vector<size_t> rightParenths;
    while (currentTermLeftParenth != string::npos || currentTermRightParenth != string::npos && currentTermLeftParenth < found && currentTermRightParenth < found)
    {
      if (currentTermLeftParenth != string::npos && currentTermLeftParenth < found)
        leftParenths.push_back(currentTermLeftParenth);
      if (currentTermRightParenth != string::npos && currentTermRightParenth < found)
        rightParenths.push_back(currentTermRightParenth);

      currentTermLeftParenth = str.find("(", currentTermLeftParenth + 1);
      currentTermRightParenth = str.rfind(")", currentTermRightParenth - 1);
    }
    if (leftParenths.size() != rightParenths.size())
      found = str.find_first_of("*/", found + 1);


    if (found != string::npos)
    {
      string s = str.substr(startPos, found - startPos);
      terms.push_back(s);
      startPos = found;
      found = str.find_first_of("*/", found + 1);
      currentTermLeftParenth = str.find(found + 1, '(');
      currentTermRightParenth = str.rfind(found, ')');
    }
  }
  terms.push_back(str.substr(startPos, string::npos));
  return terms;
}

double ApprIntegral::CalcPolynomialAndExp(string term, double num)
{
  double termConstant = 1;
  //removes the leading operator
  //if the operator is negative the termConstanticient is made negative;
  if (term.at(0) == '+' || term.at(0) == '*' || term.at(0) == '/')
    term.erase(0, 1);

  else if (term.at(0) == '-')
  {
    term.erase(0, 1);
    termConstant = -1;
  }

  double exponent = 1;
  double x = 1; //The result of the exponent calulation. Similar to the x in the term but simply left as 1 if the term is a constant 
  //determines if raised to a power. If yes then the exponent is extracted.
  size_t foundExponent = term.find('^');
  if (foundExponent != string::npos)
  {
    string exponentStr = term.substr(foundExponent + 1, string::npos);
    string subArg = SubArgument(exponentStr);

    if (exponentStr == subArg && term.find('x', foundExponent) == string::npos)
    {
      string str = term.substr(foundExponent + 1, string::npos);
      stringstream(str) >> exponent;
    }
    else
      exponent = CalcEquation(subArg, num);
    term.erase(foundExponent, string::npos);
  }

  size_t found = term.find('x');
  size_t foundE = term.find('e');
  //If 'x' is not found then the term is a constant so the termConstant is assigned the constant. 
  if (found == string::npos && foundE == string::npos)
    termConstant = CalcEquation(term, num);
  else if (foundE != string::npos)
    x = exp(exponent);
  else
  {
    //calculates the substituted value
    string subArg = SubArgument(term);
    if (subArg != term && term != "x")
      x = pow(CalcEquation(subArg, num), exponent);
    else
      x = pow(num, exponent);

    //The value preceding 'x'
    stringstream(term.substr(0, found)) >> termConstant;
  }
  if (found == string::npos && foundExponent != string::npos && foundE == string::npos)
    return pow(termConstant, exponent);

  x *= termConstant;
  return x;
}

double ApprIntegral::CalcTrig(string term, double num)
{
  int trigIndex = 0;
  double termConstant = 1;
  if (term.at(0) == '+' || term.at(0) == '*' || term.at(0) == '/')
    term.erase(0, 1);

  else if (term.at(0) == '-')
  {
    term.erase(0, 1);
    termConstant = -1;
  }

  //Gets the constant of the term
  //This is multiplication in order to preserve the potenital negative quality of termConstant
  termConstant *= GetConstantMultiplier(term);

  size_t currentTrigTerm;
  size_t argBegin = term.find('(');
  for (int i = 0; i < trig->size(); i++)
  {
    currentTrigTerm = term.rfind(trig[i], argBegin);
    trigIndex = i;
    if (currentTrigTerm != string::npos)
      break;
  }
  //calculates the value of the argument of the trig function.
  string subTerm = term.substr(currentTrigTerm, string::npos);
  string trigArg = SubArgument(subTerm);
  if (trigArg != "x" && trigArg != term)
    num = CalcEquation(trigArg, num);

  double trig = 1;
  switch (trigIndex)
  {
    case 0: trig = sin(num);
      break;
    case 1: trig = cos(num);
      break;
    case 2: trig = tan(num);
      break;
    case 3: trig = 1 / sin(num);
      break;
    case 4: trig = 1 / cos(num);
      break;
    case 5: trig = 1 / tan(num);
      break;
    case 6: trig = acos(num);
      break;
    case 7: trig = asin(num);
      break;
    case 8: trig = atan(num);
      break;
    case 9: trig = cosh(num);
      break;
    case 10: trig = sinh(num);
      break;
    case 11: trig = tanh(num);
      break;
    case 12: trig = acosh(num);
      break;
    case 13: trig = asinh(num);
      break;
    case 14: trig = atanh(num);
      break;
  }
  return termConstant * trig;
}

double ApprIntegral::CalcTerm(string term, double num)
{
  //Determines what kind of term the string term is and calls the respective function to calculate the term.
  //Current supported types are:
  //Polynomials
  //Trig Functions
  //Logarithms and Exponential functions
  double result = 1;
  double x;
  auto multiTerms = SeperateTermByMulti(term);
  for (int i = 0; i < multiTerms.size(); i++)
  {
    bool isPlainTrig = false;
    size_t funcWithExp = term.find('^');
    for (int j = 0; j < trig->size() && !isPlainTrig; j++)
    {
      if (term.find(trig[j]) != string::npos)
        isPlainTrig = true;
    }
    
    if (multiTerms[i] == "x")
      x = num;
    //if the term has an exponent it is routed through CalcPolynomialAndExp
    else if (funcWithExp == string::npos)
    {
      if (isPlainTrig)
        x = CalcTrig(multiTerms[i], num);
      else if (multiTerms[i].find("log") != string::npos || multiTerms[i].find("ln") != string::npos)
        x = CalcLog(multiTerms[i], num); // Calculates the natural log or the common log
      else if (multiTerms[i].find('e') == string::npos && multiTerms[i].find('x') == string::npos)
        x = CalcConstant(multiTerms[i]);// Calculates a constant expression witout an exponent such as 2pi or e
                                        //Constant multipliers of functions are handled in their respective functions
      else
        x = CalcPolynomialAndExp(multiTerms[i], num); // Calculates a term with an x in it without an exponent
    }
    else
      x = CalcPolynomialAndExp(multiTerms[i], num); 

    if (term.find('/') != string::npos)
      x = 1 / x; // the value of x is made reciprocal of its current value.

    result *= x;
  }
  return result;
}

//returns a substring of the argument of a function.
//if there is no argument then the orignial term is returned.
//in the case of function that cannot be simplified the first argument is returned. i.i sin(x + 1)cos(x) returns x + 1.
string ApprIntegral::SubArgument(string term)
{
  size_t leftParenth = term.find('(');
  if (leftParenth != string::npos)
  {
    size_t rightParenth = term.rfind(')');
    return term.substr(leftParenth + 1, rightParenth - leftParenth - 1);
  }
  return term;
}

double ApprIntegral::CalcLog(string term, double num)
{
  double termConstant = GetConstantMultiplier(term);
  string subArg = SubArgument(term);
  double argumentNum = CalcEquation(subArg, num);
  if (term.find("ln") != string::npos)
    return termConstant * log(argumentNum);
  return termConstant * log10(argumentNum);
}

double ApprIntegral::CalcConstant(string term)
{
  double coefficient = 1;
  if (term.at(0) == '+' || term.at(0) == '*' || term.at(0) == '/')
    term.erase(0, 1);
  else if (term.at(0) == '-')
  {
    term.erase(0, 1);
    coefficient = -1;
  }

  string constantArg;
  if (term.find('(') != string::npos)
    constantArg = SubArgument(term);
  else
    constantArg = term;

  bool endtermCoefficient = false;
  int i = 0;
  for (; !endtermCoefficient && i < constantArg.length(); i++)
    if (!isdigit(constantArg.at(i)))
      endtermCoefficient = true;

  stringstream(constantArg.substr(0, i + 1)) >> coefficient;

  if (term.find("pi") != string::npos)
    return coefficient * PI;
  
  return coefficient;

}

double ApprIntegral::GetConstantMultiplier(string term)
{
  double termConstant  = 1;
  //Gets rid of parentheses if they surrounbd the whole term
  if (term[0] == '(')
    term = SubArgument(term);

  if (isdigit(term.at(0)))
  {
    bool hasTermConstant = false;
    int i = 0;
    for (; !hasTermConstant && i < term.length(); i++)
      if (!isdigit(term.at(i)))
        hasTermConstant = true;
    stringstream(term.substr(0, i + 1)) >> termConstant;
  }
  return termConstant;
}

double ApprIntegral::CalcEquation(string equation, double num)
{
  //Seperates the equation into terms and sums the result of each seperate term.
  vector<string> terms = SeperateTermByAdd(equation);
  double result = 0;
  for (int i = 0; i < terms.size(); i++)
    result += CalcTerm(terms.at(i), num);
  return result;
}


double ApprIntegral::ApproximateIntegral(string equation, double a, double b, int numOfRects)
{
  double deltaX = (b - a) / numOfRects;
  double result = 0;
  int n = 0;
  //Approximates the integral using Simpsons Rule
  for (double i = a; i <= b || n <= numOfRects; i += deltaX , n++)
  {
    if (n == 0 || n == numOfRects)
      result += CalcEquation(equation, i);
    else if (n % 2 == 0)
      result += 2 * CalcEquation(equation, i);
    else
      result += 4 * CalcEquation(equation, i);
  }

  return (deltaX / 3) * result;
}
