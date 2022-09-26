#include <iostream>
using namespace std;

#define OPT(func) func<func##_ns::f()>
/*namespace fun_ns
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
}*/

#define declare_opt(Inline, Type, Fun, Arg) \
    Inline Type Fun Arg;                    \
    Type __##Fun##__opt;                    \
    template <bool c>                       \
    Inline Type Fun Arg

struct A
{
    declare_opt(inline, int, fun, (int));
    /*
    inline int fun(int);

    int __fun_opt;

    template <bool a>
    inline int fun(int);
    */
};

#define define_opt(Inline, Type, Class, Fun, Arg)                         \
    template <>                                                           \
    Inline Type Class::Fun<true> Arg { return __##Fun##__opt; }           \
    template <>                                                           \
    Inline Type Class::Fun<false> Arg { return __##Fun##__opt = Fun(x); } \
    namespace Fun##_ns                                                    \
    {                                                                     \
        constexpr int flag(int);                                          \
        template <class Tag>                                              \
        struct writer                                                     \
        {                                                                 \
            friend constexpr int flag(Tag) { return 0; }                  \
        };                                                                \
                                                                          \
        template <bool B, class Tag = int>                                \
        struct dependent_writer : writer<Tag>                             \
        {                                                                 \
        };                                                                \
                                                                          \
        template <                                                        \
            bool B = noexcept(flag(0)),                                   \
            int = sizeof(dependent_writer<B>)>                            \
        constexpr int f() { return B; }                                   \
    }

inline int A::fun(int x)
{
    return x * x;
}

define_opt(inline, int, A, fun, (int x));
//template <>
//inline int A::fun<true>(int x) { return __fun__opt; }
//template <>
//inline int A::fun<false>(int x) { return __fun__opt = fun(x); }

int main()
{
    A a;
    cout << a.OPT(fun)(5) << endl;
    cout << a.OPT(fun)(15) << endl;
    return 0;
}