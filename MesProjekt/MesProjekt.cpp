
#include <iostream>
#include "Grid.h"

int main()
{
    Grid Proba("Test1_4_4.txt");
   // Grid Proba("Test2_4_4_MixGrid.txt");
   // Grid Proba(0.1, 0.1, 4, 4);
    //Grid Proba("Test3_31_31_kwadrat.txt");
    
    //
    Proba.calculateDerivative();

    
    for (int i = 0; i < Proba.numberOfElements; i++) {
        Proba.calculateLocalMatrix(i);
    }
    for (int i = 0; i < Proba.numberOfElements; i++) {
        Proba.calculateHbc_VectorP(i);
    }
  
    Proba.agregate();
    Proba.iterateSolution();
   // Proba.HgCgcalculate();
   // Proba.Pcalculate();

   
}


