#ifndef META_OPT_H
#define META_OPT_H

class Normal_Call
{
};
class Dummy_Class
{
};
#ifdef USEOPT
#define ICOPT(func) func<func##_ns::f<T>()>
#define OPT(func) func<func##_ns::f<Dummy_Class>()>

#define declare_opt(Inline, Type, Fun, Arg) \
    Type __##Fun##__opt;                    \
    template <bool c>                       \
    Inline Type Fun Arg;                    \
    template <class T>

#define define_opt(Inline, Type, Class, Fun, Arg)                                       \
    template <>                                                                         \
    Inline Type Class::Fun<true> Arg { return __##Fun##__opt; }                         \
    template <>                                                                         \
    Inline Type Class::Fun<false> Arg { return __##Fun##__opt = Fun<Dummy_Class> Arg; } \
    namespace Fun##_ns                                                                  \
    {                                                                                   \
        constexpr int flag(int);                                                        \
        template <class Tag>                                                            \
        struct writer                                                                   \
        {                                                                               \
            friend constexpr int flag(Tag) { return 0; }                                \
        };                                                                              \
                                                                                        \
        template <bool B, class Tag = int>                                              \
        struct dependent_writer : writer<Tag>                                           \
        {                                                                               \
        };                                                                              \
                                                                                        \
        template <class T,                                                              \
                  bool B = noexcept(flag(0)),                                           \
                  int = sizeof(dependent_writer<B>)>                                    \
        constexpr int f() { return B; }                                                 \
    }                                                                                   \
    template <class T>
#else
#define OPT(func) func
#define ICOPT(func) func

#define declare_opt(Inline, Type, Fun, Arg)
#define define_opt(Inline, Type, Class, Fun, Arg)
#endif

#endif
