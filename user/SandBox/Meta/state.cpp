#include <iostream>
using namespace std;
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

int main()
{
    constexpr bool a = f();
    constexpr bool b = f();
    cout << a << "\n"
         << b << endl;
}