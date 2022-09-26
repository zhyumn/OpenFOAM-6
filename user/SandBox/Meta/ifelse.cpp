#include <iostream>
using namespace std;
template <bool c, typename Then, typename Else>
class IF_
{
};

template <typename Then, typename Else>
class IF_<true, Then, Else>
{
public:
    typedef Then reType;
};

template <typename Then, typename Else>
class IF_<false, Then, Else>
{
public:
    typedef Else reType;
};

struct A_base
{
};

struct A_true : public A_base
{
    void inline print() { cout << "true" << endl; }
};

struct A_false : public A_true
{
    void inline print() { cout << "false" << endl; }
};

struct A : public A_false
{
    template <bool c>
    void inline print() { IF_<c, A_true, A_false>::reType::print(); }
};

int main()
{
    A a;
    a.print<true>();
    a.print<false>();
    return 0;
}