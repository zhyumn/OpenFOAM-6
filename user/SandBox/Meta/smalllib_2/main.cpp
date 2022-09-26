#include <iostream>
using namespace std;
#define USEOPT
#include "meta_opt.H"

struct A
{
    int x;

    declare_opt_i(int, sqrx) 
    inline int sqrx();

    declare_opt_i(int, cbx) 
    inline int cbx();

    declare_opt_i(int, sqrtx) 
    inline int sqrtx();
};

define_opt_i(int, A, sqrx) 
inline int A::sqrx()
{
    cout << "sqrx is called!" << endl;
    return x * x;
}

define_opt_i(int, A, cbx) 
inline int A::cbx()
{
    cout << "cbx is called!" << endl;
    return ICOPT(sqrx)() * ICOPT(sqrx)() / x;
}

define_opt_i(int, A, sqrtx) 
inline int A::sqrtx()
{
    cout << "cbx is called!" << endl;
    return sqrt(x);
}

int calculate(A &a)
{
    return a.OPT(sqrx)() + a.OPT(cbx)();
}

int main()
{
    A a;
    a.x = 10;
    cout << a.sqrx() << endl;
    cout << a.cbx() << endl;
    cout << calculate(a) << endl;

    return 0;
}
