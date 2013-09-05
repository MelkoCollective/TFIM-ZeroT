#ifndef HEAD_PROJ_H
#define HEAD_PROJ_H

#include <iostream>
using namespace std;

#include <iomanip>
#include <vector>
#include <stack>
#include "MersenneTwister.h"
#include "matrix.h"

//Global constants here

#define alpha 2 //This is the number of replicas

class index2{

    public:

        int A;
        int B;

        index2() {A=0; B=0;};
        index2(const int a, const int b) {A=a; B=b;};
        index2 operator=(const index2 &z) {A=z.A; B=z.B; return *this;};
        void set(const int a, const int b) {A=a; B=b;};
        void print(){ cout<<"("<<A<<","<<B<<")"<<endl; };
};

class index3{

    public:

        int A;
        int B;
        int C;

        index3() {A=0; B=0; C=0;};
        index3(const int a, const int b, const int c)
              {A=a; B=b; C=c;};
        index3 operator=(const index3 &z) 
              {A=z.A; B=z.B; C=z.C; return *this;};
        void set (const int a, const int b, const int c)
              {A=a; B=b; C=c;};
        void print(){ cout<<"("<<A<<","<<B<<","<<C<<")"<<endl; };
};

class index4{

    public:

        int A;
        int B;
        int C;
        int D;

        index4() {A=0; B=0; C=0; D=0;};
        index4(const int a, const int b, const int c, const int d)
              {A=a; B=b; C=c; D=d;};
        index4 operator=(const index4 &z) 
              {A=z.A; B=z.B; C=z.C; D=z.D; return *this;};
        void set (const int a, const int b, const int c, const int d)
              {A=a; B=b; C=c; D=d;};
        void print(){ cout<<"("<<A<<","<<B<<","<<C<<","<<D<<")"<<endl; };
};




#endif
