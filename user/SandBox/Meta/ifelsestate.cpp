#include <iostream>
using namespace std;

#define OPT(func) func<func##_ns::f()>
namespace fun_ns
{
    constexpr int flag(int);
    template <class Tag>
    struct writer
    {
        friend constexpr int flag(Tag) { return 0; }
    };

    template <bool B, class Tag = int>
    struct dependent_writer : writer<Tag>
    {
    };

    template <
        bool B = noexcept(flag(0)),
        int = sizeof(dependent_writer<B>)>
    constexpr int f() { return B; }
}

struct A
{
    void inline fun();
    template <bool a>
    void inline fun() {}
};

template <>
void inline A::fun<true>() { cout << "true" << endl; }
template <>
void inline A::fun<false>() { cout << "false" << endl; }

int main()
{
    A a;
    a.OPT(fun)();
    a.OPT(fun)();
    return 0;
}