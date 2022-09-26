#include <iostream>
using namespace std;
#define USEOPT
#include "meta_opt.H"

struct A
{
    int x;
    declare_opt(inline, int, sqrx, ()) inline int sqrx();

    declare_opt(inline, int, cbx, ()) inline int cbx();
};
define_opt(inline, int, A, sqrx, ()) inline int A::sqrx()
{
    cout << "sqrx is called!" << endl;
    return x * x;
}
define_opt(inline, int, A, cbx, ()) inline int A::cbx()
{
    cout << "cbx is called!" << endl;
    return ICOPT(sqrx)() * ICOPT(sqrx)() / x;
}

int calculate(A &a)
{
    return a.OPT(sqrx)() + a.OPT(cbx)();
}

int main()
{
    A a;
    a.x = 10;
    //cout << a.sqrx() << endl;
    cout << calculate(a) << endl;

    return 0;
}
