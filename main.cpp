#include <iostream>
#include "treecode.h"
#include <chrono>
#include "omp.h"

using namespace std;

int main() {

    // using int32_t as ordering, maximum (2^7)^4 elements in matrix.
    // otherwise causes overflow.
    auto qt = new treecode(0, 0, 1., 7);
    for (auto it : qt->root->points) {it->attribute = 1.0;}
    qt->root->populate();

    scalar_t theta = 0.5;

    size_t n = qt->size * qt->size;
    auto matrix = new scalar_t[n * n];

    auto i = 0;
    auto t0 = chrono::system_clock::now();

    omp_set_num_threads(4);
#pragma omp parallel for private(i) shared(qt, theta, n, matrix)
    for (i = 0; i < n; i++) {
        traversal(qt, theta, qt->root->points[i], qt->root, n, matrix);
    }
    auto t1 = chrono::system_clock::now();

    std::cout << "matrix setting uses " <<
            chrono::duration_cast<chrono::milliseconds>(t1 - t0).count() << " milliseconds"
    << std::endl;

    delete[] matrix;
    delete qt;
    return 0;

}