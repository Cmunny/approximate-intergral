
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include "ApproximateIntegral.h"
using namespace std;

const string ApprIntegral::trig[] = { "sin","cos","tan","csc","sec","cot" };

vector<string> ApprIntegral::TermSeperater(string str)
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

double ApprIntegral::CalcPolynomial(string term, double num)
{
  double coeff = 1;
  //removes the leading operator
  //if the operator is negative the coefficient is made negative;
  if (term.at(0) == '+')
    term.erase(0, 1);

  else if (term.at(0) == '-')
  {
    term.erase(0, 1);
    coeff = -1;
  }

  int exp = 1;
  double x = 1; //The result of the exponent calulation. Similar to the x in the term but simply left as 1 if the term is a constant 
  //determines if raised to a power. If yes then the exponent is extracted.
  size_t found = term.find('^');
  if (found != string::npos)
  {
    string expStr = term.substr(found + 1, string::npos);
    string subArg = SubArgument(expStr);

    if (expStr == subArg)
    {
      string str = term.substr(found + 1, string::npos);
      stringstream(str) >> exp;
    }
    else
    {
      exp = calcEquation(subArg, num);
    }
  }

  found = term.find('x');
  //If 'x' is not found then the term is a constant so the coefficient is assigned the constant. The final result will the the coefficient.
  if (found == string::npos)
    coeff = stof(term);
  else
  {
    //calculates the substituted value
    //if there are no parentheses then the argument is assumed to be simply 'x'
    string subArg = SubArgument(term);
    if (subArg != term)
      x = pow(calcEquation(subArg, num), exp);
    else
      x = pow(num, exp);
    
    //if '/' is found then the term is qoutient like 1/x. 
    size_t divPos = term.find('/');
    if (divPos != string::npos)
    {
      found = divPos;
      // the value of x is made reciprol of its current value.
      x = 1 / x;
    }
    //The value preceding either 'x' or '/' in the term. i.e. 10x or 10/x will both result coeff == 10. 
    stringstream(term.substr(0, found)) >> coeff;
  }

  x *= coeff;
  return x;
}

double ApprIntegral::CalcTrig(string term, double num, int trigFuncIndex)
{
  double coeff = 1; //coefficient
  if (term.at(0) == '-')
    coeff *= -1;

  //Gets the coefficient of the term
  bool endCoeff = false;
  int i = 0;
  for (; !endCoeff || i < term.length(); i++)
  {
    if (!isdigit(term.at(i)))
    {
      endCoeff = true;
    }
  }
  stringstream(term.substr(0, i + 1)) >> coeff;

  //finds the term in the argument of the trig function
  
  //if there are no parentheses then the argument is assumed to be simply 'x'
  string trigArg = SubArgument(term);
    if (trigArg != "x" && trigArg != term)
      num = calcEquation(trigArg, num);
  

  double trig = 1;
  switch (trigFuncIndex)
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
  }

  return coeff * trig;
}

double ApprIntegral::termCalc(string term, double num)
{
  double x = 0;

  //Determines what kind of term the string term is and calls the respective function to calculate the term.
  //Current supported types are:
  //Polynomials
  //Trig Functions

  bool isTrig = false;
  for (int i = 0; i < 6 && !isTrig; i++)
  {
    if (term.find(trig[i]) != string::npos)
    {
      isTrig = true;
      x = CalcTrig(term, num, i);
    }
  }
  if (!isTrig)
  {
    x = CalcPolynomial(term, num);
  }
  return x;
}

//returns a substring of the argument of a function.
//if there is no argument then the orignial term is returned.
string ApprIntegral::SubArgument(string term)
{
  size_t leftParenth = term.find('(');
  if (leftParenth != string::npos)
  {
    size_t rightParenth = term.find(')');
    string subArg = term.substr(leftParenth + 1, rightParenth - leftParenth - 1);
    return subArg;
  }
  return term;
}

double ApprIntegral::calcEquation(string equation, double num)
{
  //Seperates the equation into terms and sums the result of each seperate term.
  vector<string> terms = TermSeperater(equation);
  double result = 0;
  for (int i = 0; i < terms.size(); i++)
    result += termCalc(terms.at(i), num);
  return result;
}


double ApprIntegral::approximateIntegral(string equation, double a, double b, int numOfRects)
{

  double deltaX = (b - a) / numOfRects;
  double result = 0;
  int n = 0;
  //Approximates the integral using Simpsons Rule
  for (double i = a; i <= b || n <= numOfRects; i += deltaX, n++)
  {
    if (n == 0 || n == numOfRects)
      result += calcEquation(equation, i);
    else if (n % 2 == 0)
      result += 2 * calcEquation(equation, i);
    else
      result += 4 * calcEquation(equation, i);
  }

  return (deltaX / 3) * result;
}
