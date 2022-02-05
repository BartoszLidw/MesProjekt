#pragma once
#include <iostream>
#include <iomanip>
class Node
{
public:
	double x, y;
	int brzeg;
	double temp = 100.0;
	Node();
	Node(double x, double y);
	void print();
};

