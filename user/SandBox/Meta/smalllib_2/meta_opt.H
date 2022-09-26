#ifndef META_OPT_H
#define META_OPT_H

struct Normal_Call
{
    constexpr static int f() { return false; }
};
#ifdef USEOPT
#define ICOPT(func) func<T::f()>
#define OPT(func) func<func##_ns::OPT_Call::f()>

#define declare_opt_i(Type, Fun) \
    Type __##Fun##__opt;         \
    template <bool c>            \
    inline Type Fun();           \
    template <class T = Normal_Call>

#define define_opt_i(Type, Class, Fun)                                                  \
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
        struct OPT_Call                                                                 \
        {                                                                               \
            template <bool B = noexcept(flag(0)), int N = sizeof(dependent_writer<B>)>  \
            constexpr static int f() { return B; }                                      \
        };                                                                              \
    }                                                                                   \
    template <>                                                                         \
    inline Type Class::Fun<true> () { return __##Fun##__opt; }                         \
    template <>                                                                         \
    inline Type Class::Fun<false> () { return __##Fun##__opt = Fun<Fun##_ns::OPT_Call> (); } \
    template <class T>
#else
#define OPT(func) func
#define ICOPT(func) func

#define declare_opt_i(Type, Fun)
#define define_opt_i(Type, Class, Fun)
#endif

#endif
