#pragma once
#include <math.h>
#include "Element.h"
#include "Node.h"
#include <string.h>
#include <string>

class Grid
{
	public:
		//Dane 
		double INITIAL_TEMPERATURE;
		double SYMULATION_TIME;
		double STMULATION_STEP_TIME;
		double AMBIENT_TEMPERATURE;
		double ALFA;
		double SPECIFIC_HEAT;
		double CONDUCTIVITY;
		double DENSITY;

	double height, lenght;
	int nHeight, nLenght;
	int numberOfElements, numberOfNods;
	Node* marks;
	Element* elements;
	double pointKsiBok[4][4] = { 0.0 };
	double pointEtaBok[4][4] = { 0.0 };
	double pointKsztalt[4][4] = { 0.0 };
	double weight[2] = { 1, 1 };
	double* vectorP;
	double* Pbrzeg;
	double** Hg;
	double** HBCg;
	double** Cg;
	double** Hg_Cg;

	Grid(double H, double L, int nH, int nB);
	Grid(std::string nazwa);
	void calculateDerivative();
	void HgCgcalculate();
	void calculateLocalMatrix(int i);
	void calculateHbc_VectorP(int i);
	void agregate();
	void Pcalculate();
	void iterateSolution();
	void getTemperature( double** A, double* T);
	void print();
	void printTemperatures(int x);
};

