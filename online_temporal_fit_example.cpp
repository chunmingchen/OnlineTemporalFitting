#include <algorithm>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <cstdio>
#include <string>
#include <string.h>
#include <iostream>
#include "VectorMatrix.h"
#include "omp.h"

#include "online_quad_bezier_fit.h"
using namespace std;


inline VECTOR3 &operator +=(VECTOR3 & v0, const VECTOR3 & v1)
{
    v0[0]+=v1[0]; v0[1]+=v1[1]; v0[2]+=v1[2];
    return v0;
}

inline VECTOR3 operator*(const VECTOR3 & v0, const VECTOR3 & v1)
{
    return VECTOR3(v0[0]*v1[0], v0[1]*v1[1], v0[2]*v1[2]);
}

inline VECTOR3 sqrt(const VECTOR3 &x)
{
    return  VECTOR3(sqrt(x[0]), sqrt(x[1]), sqrt(x[2]));
}

inline VECTOR3 max(const VECTOR3 &x, const VECTOR3 &y)
{
    return VECTOR3(max(x[0], y[0]), max(x[1], y[1]), max(x[2], y[2]));
}

inline std::ostream& operator<<(std::ostream& os, const VECTOR3 &vec)
{
    os <<  "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
    return os;
}


////////////////////////////////////////////////////////


class ErrorModeling{

public:
    void test_online_scalar()
    {
        float seq[3][10] = {{1, 2, 3, 100, 5, 6, 7, 8, 9, 10},
                          {1, 8, 27, 100, 125, 196, 343, 512, 729, 1111},
                          {1, 4, 9, 16, 25, 36, 49, 64, 81, 100} };
        int n = 10;
        int sampling = n-1;
        int i, s;

        vector<OnlineQuadBezierFit<float> > onlineAry(3);
        vector<float> rmsErrAry;  // mean square error;  xyz, dim
        vector<float> quadBezierAry; // control point


        //////////////////////////////////////
        /// test online version
        //////////////////////////////////////
        for (s = 0; s<3; s++)
            for (i=0; i<n; i++)
            {
                onlineAry[s].addData(seq[s][i], (double)i/(n-1));
            }

        fitOnlineQuadBezier(quadBezierAry, rmsErrAry, onlineAry, n);
        for (s = 0; s < 3; s++)
        {
            cout << "online fitting result for scalar sequence:  " << s << endl;
            cout << "P0 = " << seq[s][0] << ",  P1 = " << seq[s][n-1] << endl;
            cout << "Bezier control value=" <<  quadBezierAry[s] << endl;
            cout << "Stderr=" << rmsErrAry[s] << endl;
            cout << "============================" << endl;
        }

    }


    void test_online_vector()
    {
        float seq[3][10] = {{1, 2, 3, 100, 5, 6, 7, 8, 9, 10},
                          {1, 8, 27, 100, 125, 196, 343, 512, 729, 1111},
                          {1, 4, 9, 16, 25, 36, 49, 64, 81, 100} };
        int n = 10;
        int sampling = n-1;
        int i, s;

        vector<OnlineQuadBezierFit<VECTOR3> > onlineAry(1);
        vector<VECTOR3> rmsErrAry;  // mean square error;  xyz, dim
        vector<VECTOR3> quadBezierAry; // control point


        //////////////////////////////////////
        /// test online version
        //////////////////////////////////////
        for (i=0; i<n; i++)
        {
            onlineAry[0].addData(VECTOR3(seq[0][i],seq[1][i],seq[2][i]), (double)i/(n-1));
        }

        fitOnlineQuadBezier(quadBezierAry, rmsErrAry, onlineAry, n);
        cout << "online fitting result for the first sequence:  " << endl;
        cout << "P0 = " << VECTOR3(seq[0][0], seq[1][0], seq[2][0]) << ",  P1 = " << VECTOR3(seq[0][n-1], seq[1][n-1], seq[2][n-1]) << endl;
        cout << "Bezier control value=" <<  quadBezierAry[0] << endl;
        cout << "Stderr=" << rmsErrAry[0] << endl;

    }
};

int main(int argc, const char **argv) {

    // ErrorModeling<float, 1>em; // scalar fields
    ErrorModeling em; // vector with 3 elements
    em.test_online_scalar();
    em.test_online_vector();
    return 0;

}
