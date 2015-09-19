//
// Created by lurker on 9/19/15.
//

#include "treecode.h"

attribute_t treecode::getAttribute(scalar_t x0, scalar_t y0) noexcept {
    auto col = int((x0 - x) * size);
    auto row = int((y0 - y) * size);
    return root->points[row * size + col]->attribute;
}

void treecode::setAttribute(attribute_t attr, scalar_t x0, scalar_t y0) noexcept {
    auto col = int((x0 - x) * size);
    auto row = int((y0 - y) * size);
    root->points[row * size + col]->attribute = attr;
}


