#ifndef MDCT_HPP
#define MDCT_HPP

#include <math.h>
#include <vector>

void easy_IMDCT(vector<double> &X, vector<double> &y, vector<double> &window, int n) {
    for (int p = 0; p < n; ++p) {
        y[p] = 0.0;
        for (int m = 0; m < n / 2; ++m) {
            y[p] += X[m] * cos(M_PI / 2 / n * (2 * p + 1 + n / 2) * (2 * m + 1));
        }
        y[p] *= window[p];
    }
}

#endif
