#include <iostream>
#include "treecode.h"

using namespace std;

int main() {

    auto qt = new treecode(0, 0, 1., 6);
    for (auto it : qt->root->points) {it->attribute = 1.0;}
    qt->root->populate();

    scalar_t theta = 0.2;
    int n = qt->size * qt->size;

    auto matrix = new scalar_t[n * n];

    for (auto i = 0; i < n; i++) {

        traversal(qt, theta, qt->root->points[i], qt->root, n, matrix);

    }

    delete[] matrix;
    delete qt;
    return 0;

}