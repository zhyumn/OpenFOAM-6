#include<iostream>
using namespace std;

struct MyInt
{
    int a_;
    inline operator int& ()
    {
        cout<<"MyInt::int"<<endl;
        return a_;
    }
    inline int  operator =(int a)
    {
        return a_ = a;
    }
    inline void* operator new(std::size_t size)
    {
        cout << "haha:" << size << endl;
        return ::operator new(size);

    }

    inline void* operator new[](std::size_t size)
    {
        cout << "haha:" << size << endl;
        return ::operator new(size);
    }
        inline void operator delete(void* p)
    {
        cout << "hee:" << endl;
        return ::operator delete(p);
    }

    inline void operator delete[](void* p, size_t a)
    {
        cout << "hee:[]" << a << endl;
        return ::operator delete(p);
    }
};

struct MyInt_2 :public MyInt
{
    int b;
    int fun() { return 0; }

};

class intlist
{
public:
    int* s_;
    intlist(int a) { s_ = new int[a]; }
};

class Intlist
{
public:
    MyInt_2* s_;
    Intlist(int a) { s_ = new MyInt_2[a]; }
    ~Intlist() { delete[] s_; }
};


int main()
{
    Intlist a(100);
    a.s_[0].a_ = 100.1;
    a.s_[0].b = 12;
    cout << a.s_[0] << "," << a.s_[0] << endl;
    return 0;

}