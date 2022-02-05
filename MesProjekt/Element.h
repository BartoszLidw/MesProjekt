#pragma once
#include <iostream>

#include "MatrixH.h"
class Element
{
public:
	int ID[4];
	double Hwynik[4][4] = { 0.0 };
	double HCwynik[4][4] = { 0.0 };
	MatrixH H[4];
	MatrixH Hbc[4];

	//Element();
	//Element(int id1, int id2, int id3, int id4);
	void print();
};

