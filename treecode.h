//
// Created by lurker on 9/19/15.
//

#ifndef QUADTREE_TREECODE_H
#define QUADTREE_TREECODE_H

#include "quadtree.h"

typedef int level_t;

/*
 * treecode structure
 * @members
 *
 * @x : x coordinate of left bottom corner
 * @y : y coordinate of left bottom corner
 * @length : length of side
 * @size : number of points per slice, it is 2^n
 * @root : root pointer of quadtree
 *
 */
typedef struct treecode {
    scalar_t  x;
    scalar_t  y;
    scalar_t  length;
    level_t   size;
    quadtree*  root;

    /*
     * treecode constructor
     * @params
     *
     */
    explicit treecode(scalar_t _x_start, scalar_t _y_start,
              scalar_t _length, level_t _level)  noexcept :
            x(_x_start), y(_y_start), length(_length), size(1 << _level){
        ord_t order = 0;
        root = new quadtree(length, x, y);
        root->setStatus(Status::ROOT);
        for (auto i = 1; i <= size; i++) {
            for (auto j = 1; j <= size; j++) {
                // assign points with
                auto ptr = make_shared<point>(
                        x + (2.0 * i - 1.) * length/(2.0 * size),
                        y + (2.0 * j - 1.) * length/(2.0 * size)
                );
                ptr->id = order++;
                root->addPoint(ptr);
            }
        }
    }
    /*
     * treecode destructor
     */
    ~treecode() noexcept {delete root;}

    /*
     * return attribute at grid containing (x0, y0)
     */
    attribute_t getAttribute(scalar_t x0, scalar_t y0) noexcept ;

    /*
     * set attribute at the grid containing (x0, y0)
     */
    void setAttribute(attribute_t attr, scalar_t x0, scalar_t y0) noexcept ;

} treecode;

/*
 *  return distance between two locations
 *  @params
 *
 *  @x0 : x coordinate
 *  @y0 : y coordinate
 *  @x1 : x coordinate
 *  @y1 : y coordinate
 */
inline scalar_t distance(scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

/*
 * return line integral between two locations
 *
 * todo: provide accurate line integral
 *
 * @params
 *
 * @tree : treecode structure
 * @x0 : x coordinate
 * @y0 : y coordinate
 * @x1 : x coordinate
 * @y1 : y coordinate
 */
inline scalar_t integral(treecode *tree, scalar_t x0, scalar_t y0, scalar_t x1, scalar_t y1) noexcept {
    return distance(x0, y0, x1, y1) * tree->root->attribute;
}

/*
 * return evaluation of kernel over a grid
 *
 * todo: provide accurate evaluation function
 *
 * @params
 */
inline scalar_t eval(treecode *tree, scalar_t x0, scalar_t y0, quadtree* branch_ptr) noexcept {
    return exp(-integral(tree, x0, y0, branch_ptr->x, branch_ptr->y ))
           * branch_ptr->length * branch_ptr->length
           /distance(x0, y0, branch_ptr->x, branch_ptr->y);
}
/*
 * traversal quadtree
 *
 * todo: pre-order traversal
 *
 * @params
 *
 * @tree
 * @theta
 * @point_ptr
 * @branch_ptr
 * @n
 * @matrix_ptr
 */
inline void traversal(treecode *tree, scalar_t& theta,
                      shared_ptr<point> point_ptr, quadtree* branch_ptr,
                      int& n, scalar_t* matrix_ptr) noexcept {

    auto d = distance(point_ptr->x, point_ptr->y, branch_ptr->x, branch_ptr->y);

    if (branch_ptr->status == Status::LEAF || branch_ptr->length/ d < theta) {
        // if reaches a leaf or in far field
        auto ret = eval(tree, point_ptr->x, point_ptr->y, branch_ptr);
        for (auto point : branch_ptr->points) {
            if (point->id != point_ptr->id) {
                matrix_ptr[n * point->id + point_ptr->id] = ret / branch_ptr->points.size();
            }
            else {
                matrix_ptr[n * point->id + point_ptr->id] =
                        2 * M_PI * (1 - exp(-branch_ptr->attribute * branch_ptr->length))/branch_ptr->attribute;
            }
        }
    }
    else {
        // else recursively traverse.
        for (auto child : branch_ptr->children) {
            if (child->status != Status::EMPTY) {
                traversal(tree, theta, point_ptr, child.get(), n, matrix_ptr);
            }
        }
    }
}

#endif //QUADTREE_TREECODE_H
