#include "MatrixH.h"
void MatrixH::print()
{
    
    std::cout << std::endl;
    std::cout << " To jest elem : \n";
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++) {
            std::cout << "  " << H[i][j];
        }
        std::cout << std::endl;
    }
}