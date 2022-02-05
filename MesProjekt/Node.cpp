#include "Node.h"
Node::Node() :x(0.0), y(0.0)
{
}
Node::Node(double x, double y)
{
   this-> x = x;
   this-> y = y;
   
}

void Node::print()
{
    std::cout << std::setprecision(4);
    std::cout << "X: " << x << "\t Y: " << y << "\n";
}