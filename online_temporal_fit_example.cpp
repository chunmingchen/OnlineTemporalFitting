#include <algorithm>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <cstdio>
#include <string>
#include <string.h>
#include <iostream>
#include <math.h>

#include "online_quad_bezier_fit.h"
using namespace std;


////////////////////////////////////////////////////////


class ErrorModeling{

public:
    void test_online_scalar()
    {
        // generate 3 scalar sequences of length 10
        float seq[3][10] = {{1, 2, 3, 100, 5, 6, 7, 8, 9, 10},
                          {1, 8, 27, 100, 125, 196, 343, 512, 729, 1111},
                          {1, 4, 9, 16, 25, 36, 49, 64, 81, 100} };
        int n = 10;
        int sampling = n-1;
        int t, s;

        vector<OnlineQuadBezierFit<float> > onlineAry(3);
        vector<float> rmsErrAry;  // mean square error;  xyz, dim
        vector<float> quadBezierAry; // control point


        //////////////////////////////////////
        /// incremental fitting
        //////////////////////////////////////
        // for each time step 
        for (t=0; t<n; t++)
            // assign for each data sequence
            for (s = 0; s<3; s++)
            {
                onlineAry[s].addData(seq[s][t], (double)t/(n-1));
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

};

int main(int argc, const char **argv) {

    // ErrorModeling<float, 1>em; // scalar fields
    ErrorModeling em; // vector with 3 elements
    em.test_online_scalar();
    return 0;

}
