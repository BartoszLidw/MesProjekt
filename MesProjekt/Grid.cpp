#include "Grid.h"
#include <math.h>
#include <fstream>
Grid::Grid(double H, double L, int nH, int nL)
{
    this->height = H;
    this->lenght = L;
    this->nHeight = nH;
    this->nLenght = nL;
    this->numberOfNods = nH * nL;
    this->numberOfElements = (nH - 1) * (nL - 1);
    marks = new Node[numberOfNods];
    int x = 0;
    for (int i = 0; i < nLenght; i++)
    {
        for (int j = 0; j < nHeight; j++)
        {
            marks[x] =  Node(i * (lenght / (nLenght - 1)), j * (height / (nHeight - 1)));
            marks[x].temp = 100.0;
            if (i == 0 || j == 0 || i == nLenght - 1 || j == nHeight - 1)
            {
                marks[x].brzeg = 1;
            }
            else
            {
                marks[x].brzeg = 0;
            }
                 x++;
        }

    }
    elements = new Element[numberOfElements];
    int alfa = 0;
    for (int i = 0; i < numberOfElements; i++)
    {
        if (!(i % (nH - 1))) {
            alfa++;
        }
        elements[i].ID[0] = i + alfa;
        elements[i].ID[1] = elements[i].ID[0] + nH;
        elements[i].ID[2] = elements[i].ID[1] + 1;
        elements[i].ID[3] = elements[i].ID[0] + 1;
      //  std::cout << "(" << elements[i].ID[0] - 1 << " " << elements[i].ID[1] - 1 << " " << elements[i].ID[2] -1 << " " << elements[i].ID[3] -1 << ")\n";
    }
    Hg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        Hg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            Hg[i][j] = 0.0;
    }
    HBCg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        HBCg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            HBCg[i][j] = 0.0;
    }
    Cg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        Cg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            Cg[i][j] = 0.0;
    }
    Hg_Cg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        Hg_Cg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            Hg_Cg[i][j] = 0.0;
    }

    vectorP = new double[numberOfNods];
    for (int j = 0; j < numberOfNods; j++)
    {
        vectorP[j] = 0.0;
    }

    Pbrzeg = new double[numberOfNods];
    for (int j = 0; j < numberOfNods; j++)
    {
       Pbrzeg[j] = 0.0;
    }
}
Grid::Grid(std::string nazwa)
{
  
    std::fstream odczyt(nazwa);
   
    std::string pomocniczy;
    
        odczyt >> pomocniczy >> this->SYMULATION_TIME ;
        std::cout << pomocniczy << SYMULATION_TIME <<" \n";
        odczyt >> pomocniczy >> this->STMULATION_STEP_TIME;          
        std::cout << pomocniczy << STMULATION_STEP_TIME << " \n";
        odczyt >> pomocniczy >> this->CONDUCTIVITY;
        std::cout << pomocniczy << CONDUCTIVITY << " \n";
        odczyt >> pomocniczy >> this->ALFA;
        std::cout << pomocniczy << ALFA << " \n";
        odczyt >> pomocniczy >> this->AMBIENT_TEMPERATURE;
        std::cout << pomocniczy << AMBIENT_TEMPERATURE << " \n";
        odczyt >> pomocniczy >> this->INITIAL_TEMPERATURE;
        std::cout << pomocniczy << INITIAL_TEMPERATURE << " \n";
        odczyt >> pomocniczy >> this->DENSITY;
        std::cout << pomocniczy << DENSITY << " \n";
        odczyt >> pomocniczy >> this->SPECIFIC_HEAT;
        std::cout << pomocniczy << SPECIFIC_HEAT << " \n";
        odczyt >> pomocniczy >> this->numberOfNods;
        std::cout << pomocniczy << numberOfNods << " \n";
        odczyt >> pomocniczy >> this->numberOfElements;
        std::cout << pomocniczy << numberOfElements << " \n";
        marks = new Node[numberOfNods];
        int x = 0;

        odczyt >> pomocniczy;
        std::cout << pomocniczy << " \n";
        for (int i = 0; i < numberOfNods; i++)
        {
            std::string x, y;
            std::string linia;
       

            odczyt >> linia >> x >> y;
            x.pop_back();
            y.pop_back();
            double num_double_x = std::stod(x);
            double num_double_y = std::stod(y);
            std::cout<<std::setprecision(10) << linia << " " << x << " " << y << " \n";
                marks[i] = Node(num_double_x, num_double_y);
                marks[i].temp = INITIAL_TEMPERATURE;
                marks[i].brzeg = 0;
        }
        elements = new Element[numberOfElements];
  
        odczyt >> pomocniczy;
        std::cout << pomocniczy << " \n";
        for (int i = 0; i < numberOfElements; i++)
        {
            std::string as, bs, cs, ds;
            int a, b, c, d;
            odczyt >> pomocniczy>>  as >> bs >> cs >>ds;
            a = std::stoi(as);
            b = std::stoi(bs);
            c = std::stoi(cs);
            d = std::stoi(ds);
          //  std::cout << a << " " << b << " " << c << " " << d << "\n";
            elements[i].ID[0] = a;
            elements[i].ID[1] = b;
            elements[i].ID[2] = c;
            elements[i].ID[3] = d;
              std::cout << "(" << elements[i].ID[0] - 1 << " " << elements[i].ID[1] - 1 << " " << elements[i].ID[2] -1 << " " << elements[i].ID[3] -1 << ")\n";
        }
        odczyt >> pomocniczy;
        std::cout << pomocniczy << " \n";
     
        while (true) //pêtla nieskoñczona
        {
            std::string as;
            int a;
            odczyt >> as;
            if (!odczyt.fail()) {
                a = std::stoi(as);
                marks[a -1].brzeg = 1;
                std::cout << a << " \n";
            }
            else
                break;
           

        }

    Hg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        Hg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            Hg[i][j] = 0.0;
    }
    HBCg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        HBCg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            HBCg[i][j] = 0.0;
    }
    Cg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        Cg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            Cg[i][j] = 0.0;
    }
    Hg_Cg = new double* [numberOfNods];
    for (int i = 0; i < numberOfNods; i++)
    {
        Hg_Cg[i] = new double[numberOfNods];
        for (int j = 0; j < numberOfNods; j++)
            Hg_Cg[i][j] = 0.0;
    }

    vectorP = new double[numberOfNods];
    for (int j = 0; j < numberOfNods; j++)
    {
        vectorP[j] = 0.0;
    }

    Pbrzeg = new double[numberOfNods];
    for (int j = 0; j < numberOfNods; j++)
    {
        Pbrzeg[j] = 0.0;
    }
}
void Grid::calculateDerivative()
{   
    double pointksi[4] = { -0.5773502692, 0.5773502692 ,0.5773502692 ,-0.5773502692 };
    double pointeta[4] = { -0.5773502692, -0.5773502692 ,0.5773502692 ,0.5773502692 };
   
    for (int i = 0; i < 4; i++)
    {
        pointEtaBok[i][0] = -0.25 * (1.0 - pointeta[i]);
        pointEtaBok[i][1] = -0.25 * (1.0 + pointeta[i]);
        pointEtaBok[i][2] = 0.25 * (1.0 + pointeta[i]);
        pointEtaBok[i][3] = 0.25 * (1.0 - pointeta[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        pointKsiBok[i][0] = -0.25 * (1.0 - pointksi[i]);
        pointKsiBok[i][1] = 0.25 * (1.0 - pointksi[i]);
        pointKsiBok[i][2] = 0.25 * (1.0 + pointksi[i]);
        pointKsiBok[i][3] = -0.25 * (1.0 + pointksi[i]);
    }
    for (int i = 0; i < 4; i++)
    {
        pointKsztalt[i][0] = 0.25 * (1.0 - pointksi[i]) * (1.0 - pointeta[i]);
        pointKsztalt[i][1] = 0.25 * (1.0 + pointksi[i])* (1.0 - pointeta[i]);
        pointKsztalt[i][2] = 0.25 * (1.0 + pointksi[i])* (1.0 + pointeta[i]);
        pointKsztalt[i][3] = 0.25 * (1.0 - pointksi[i])* (1.0 + pointeta[i]);
    }
}
void Grid::calculateLocalMatrix(int index)
{

    double x[4] = { marks[elements[index].ID[0] - 1].x  , marks[elements[index].ID[1] - 1].x, marks[elements[index].ID[2] - 1].x,marks[elements[index].ID[3] - 1].x };
    double y[4] = { marks[elements[index].ID[0] - 1].y  , marks[elements[index].ID[1] - 1].y, marks[elements[index].ID[2] - 1].y,marks[elements[index].ID[3] - 1].y };
    double dxdKsi, dydKsi, dxdEta, dydEta;
    double detJ[4];
    double reservedetJ[4];
    double dx[4];
    double dy[4];
    double Hp[4][4] = { 0.0 };
    double Hc[4][4] = { 0.0 };
    for (int k = 0; k < 4; k++) 
    {
       
        for (int i = 0; i < 4; i++)
        {
            dxdKsi = pointKsiBok[k][0] * x[0] + pointKsiBok[k][1] * x[1] + pointKsiBok[k][2] * x[2] + pointKsiBok[k][3] * x[3];
            dydKsi = pointKsiBok[k][0] * y[0] + pointKsiBok[k][1] * y[1] + pointKsiBok[k][2] * y[2] + pointKsiBok[k][3] * y[3];
            dxdEta = pointEtaBok[k][0] * x[0] + pointEtaBok[k][1] * x[1] + pointEtaBok[k][2] * x[2] + pointEtaBok[k][3] * x[3];
            dydEta = pointEtaBok[k][0] * y[0] + pointEtaBok[k][1] * y[1] + pointEtaBok[k][2] * y[2] + pointEtaBok[k][3] * y[3]; 

            detJ[k] = (dxdKsi * dydEta) - (dydKsi * dxdEta); //JACOBIAN 
            reservedetJ[i] = 1 / ((dxdKsi * dydEta) - (dydKsi * dxdEta)); //JACOBIAN ODWROTNY
            //pochodne po dEta i dKsi w punktach ca³kowania
            dx[i] = reservedetJ[i] * ((dydEta * pointKsiBok[k][i]) - (dydKsi * pointEtaBok[k][i]));
            dy[i] = reservedetJ[i] * ((-dxdEta * pointKsiBok[k][i]) + (dxdKsi * pointEtaBok[k][i]));     
        }
        
            for (int zmienna = 0; zmienna < 4; zmienna++)
            {
                for (int j = 0; j < 4; j++)
                {
                    Hp[zmienna][j] += detJ[k] * (dx[zmienna] * dx[j] + dy[zmienna] * dy[j]) * CONDUCTIVITY;
                    Hc[zmienna][j] += (SPECIFIC_HEAT * DENSITY * detJ[k] * (pointKsztalt[k][j] * pointKsztalt[k][zmienna]));
                }
            }
        
      
    }
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++) {
            
            elements[index].HCwynik[i][j] = Hc[i][j];
            elements[index].Hwynik[i][j] = Hp[i][j];    
        }
     
    }
}
void Grid::calculateHbc_VectorP(int index)
{
    double x[4] = { marks[elements[index].ID[0] - 1].x  , marks[elements[index].ID[1] - 1].x, marks[elements[index].ID[2] - 1].x,marks[elements[index].ID[3] - 1].x };
    double y[4] = { marks[elements[index].ID[0] - 1].y  , marks[elements[index].ID[1] - 1].y, marks[elements[index].ID[2] - 1].y,marks[elements[index].ID[3] - 1].y };
    double pointKsiBok[4][2] = { {-0.5773502692, 0.5773502692}, {1.0, 1.0 }, {0.5773502692,-0.5773502692 },{-1.0, -1.0 } };
    double pointEtaBok[4][2] = { {-1, -1}, {0.5773502692, -0.5773502692}, {1.0,1.0 },{0.5773502692, -0.5773502692 } };
    double N[4][2][4];
    double detJ[4];

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            N[i][j][0] = 0.25 * (1.0 - pointKsiBok[i][j]) * (1 - pointEtaBok[i][j]);
            N[i][j][1] = 0.25 * (1.0 + pointKsiBok[i][j]) * (1 - pointEtaBok[i][j]);
            N[i][j][2] = 0.25 * (1.0 + pointKsiBok[i][j]) * (1 + pointEtaBok[i][j]);
            N[i][j][3] = 0.25 * (1.0 - pointKsiBok[i][j]) * (1 + pointEtaBok[i][j]);
        }
    }
  
    detJ[0] = sqrt(((x[1] - x[0]) * (x[1] - x[0])) + ((y[1] - y[0]) * (y[1] - y[0]))) / 2;
    detJ[1] = sqrt(((x[2] - x[1]) * (x[2] - x[1])) + ((y[2] - y[1]) * (y[2] - y[1]))) / 2;
    detJ[2] = sqrt(((x[3] - x[2]) * (x[3] - x[2])) + ((y[3] - y[2]) * (y[3] - y[2]))) / 2;
    detJ[3] = sqrt(((x[0] - x[3]) * (x[0] - x[3])) + ((y[0] - y[3]) * (y[0] - y[3]))) / 2;
    for (int side = 0; side < 4; side++)
    {
        bool sidecheck = false;
        if (side == 0) {
            if (marks[elements[index].ID[0] - 1].brzeg == 1 && marks[elements[index].ID[1] - 1].brzeg == 1)
            {
                sidecheck = true;
            }    
        }
        if (side == 1) {
            if (marks[elements[index].ID[1] - 1].brzeg == 1 && marks[elements[index].ID[2] - 1].brzeg == 1)
            {
                sidecheck = true;
            }
        }
        if (side == 2) {
            if (marks[elements[index].ID[2] - 1].brzeg == 1 && marks[elements[index].ID[3] - 1].brzeg == 1)
            {
                sidecheck = true;
            }
        }
        if (side == 3) {
            if (marks[elements[index].ID[3] - 1].brzeg == 1 && marks[elements[index].ID[0] - 1].brzeg == 1)
            {
                sidecheck = true;
            }
        }

        if (sidecheck == true)
        {
            double pc1[4][4] = { 0.0 };
            double pc2[4][4] = { 0.0 };
            double Plocal[4] = { 0.0, 0.0 ,0.0 ,0.0 };
            for (int j = 0; j < 4; j++)
            {
                Plocal[j] +=  N[side][0][j] * ALFA * AMBIENT_TEMPERATURE
                * detJ[side];
                Plocal[j] +=  N[side][1][j] * ALFA * AMBIENT_TEMPERATURE
                * detJ[side];
                for (int k = 0; k < 4; k++)
                {
                    pc1[j][k] = pc1[j][k] + (N[side][0][j] * N[side][0][k] * ALFA * detJ[side]);
                    pc2[j][k] = pc2[j][k] + (N[side][1][j] * N[side][1][k] * ALFA * detJ[side]);
                }
            }  
           /* for (int j = 0; j < 4; j++)
            {
                std::cout << Plocal[j] << " ";
            }
            std::cout << " lokal \n";*/
            for (int j = 0; j < 4; j++)
            {             
                Pbrzeg[(elements[index].ID[j] - 1)] +=  Plocal[j];          
                for (int k = 0; k < 4; k++)
                {
                    HBCg[elements[index].ID[j] - 1][elements[index].ID[k] - 1] = HBCg[elements[index].ID[j] - 1][elements[index].ID[k] - 1] + pc1[j][k] + pc2[j][k];
                }
              
            }

        }
    } 
  
}
void Grid::agregate()
{
    for (int i = 0; i < numberOfNods; i++)
    {
        for (int j = 0; j < numberOfNods; j++)
        {
            Hg[i][j] = Hg[i][j] + HBCg[i][j];

        }

    }
    for (int i = 0; i < numberOfElements; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                Hg[(elements[i].ID[j]) - 1][(elements[i].ID[k]) -1 ] += elements[i].Hwynik[j][k];
                Cg[(elements[i].ID[j]) - 1][(elements[i].ID[k]) - 1] += elements[i].HCwynik[j][k];
            }
        }
    }
  /*  
    for (int i = 0; i < numberOfNods; i++)
    {
        for (int j = 0; j < numberOfNods; j++)
        { 
                std::cout<< std::setprecision(4) << Hg[i][j] << " ";
             
        }
        std::cout << std::endl;
    }*/

 
    //std::cout << " \nHg + HBc\n";
    //for (int i = 0; i < numberOfNods; i++)
    //{
    //    for (int j = 0; j < numberOfNods; j++)
    //    {
    //        std::cout << std::setprecision(4) << Hg[i][j] << " ";

    //    }
    //    std::cout << std::endl;
    //}
  /*  for (int i = 0; i < numberOfNods; i++)
    {
        for (int j = 0; j < numberOfNods; j++)
        {
            std::cout << std::setprecision(4) << Cg[i][j] << " ";

        }
        std::cout << std::endl;
    }
  */
}
void Grid::HgCgcalculate()
{
    for (int i = 0; i < numberOfNods; i++)
    {
        for (int j = 0; j < numberOfNods; j++)
        {
            Hg_Cg[i][j] = Hg[i][j] + Cg[i][j] / STMULATION_STEP_TIME;
        }
    }
  /*  for (int i = 0; i < numberOfNods; i++)
    {
        for (int j = 0; j < numberOfNods; j++)
        {
            std::cout << std::setprecision(4) << Hg_Cg[i][j] << " ";

        }
        std::cout << std::endl;
    }*/
}
void Grid::Pcalculate()
{
    for (int i = 0; i < numberOfNods; i++)
    {
        double tmp = 0.0;
        for (int j = 0; j < numberOfNods; j++)
        {
            tmp +=  marks[j].temp * (Cg[i][j] / STMULATION_STEP_TIME);  
        }
        vectorP[i] =  tmp + Pbrzeg[i];  
    }
}
void Grid::getTemperature( double** MATRIXH_P, double* temp)//metoda eliminacji gaussa
{

    const double eps = 1e-12;

    int i, j, k;
    double m, s;

   

    for (i = 0; i < numberOfNods - 1; i++)
    {
        for (j = i + 1; j < numberOfNods; j++)
        {
            if (fabs(MATRIXH_P[i][i]) < eps) 
                break;
            m = -MATRIXH_P[j][i] / MATRIXH_P[i][i];
            for (k = i + 1; k <= numberOfNods; k++)
                MATRIXH_P[j][k] += m * MATRIXH_P[i][k];
        }
    }
    for (i = numberOfNods - 1; i >= 0; i--)
    {
        s = MATRIXH_P[i][numberOfNods];    
        for (j = numberOfNods - 1; j >= i + 1; j--) {           
            s -= MATRIXH_P[i][j] * temp[j];
        }
        if (fabs(MATRIXH_P[i][i]) < eps)
            break;
        temp[i] = s / MATRIXH_P[i][i];
    }
 
}
void Grid::iterateSolution()
{
    int currentTime = 0;
    printTemperatures(currentTime);
    for (int index = 0; index < SYMULATION_TIME / STMULATION_STEP_TIME; index++) {
        HgCgcalculate();
        Pcalculate();

        double* temp1 = new double[numberOfNods];
        for (int i = 0; i < numberOfNods; i++)
            temp1[i] = 0.0;

        double** H_P = new double* [numberOfNods];      
        for (int i = 0; i < numberOfNods; i++)
        {
            H_P[i] = new double[numberOfNods + 1];
            for (int j = 0; j < numberOfNods; j++)
                H_P[i][j] = 0.0;
        }

        for (int i = 0; i < numberOfNods; i++)
        {
            for (int j = 0; j < numberOfNods; j++)
            {
                H_P[i][j] = Hg_Cg[i][j];
            }
            H_P[i][numberOfNods] = vectorP[i];
        }

        currentTime += STMULATION_STEP_TIME;

        getTemperature( H_P, temp1);

        for (int i = 0; i < numberOfNods; i++)
        {
            marks[i].temp = temp1[i];
        }
        printTemperatures(currentTime);
    }
}
void Grid::print()
{
    for (int i = 0; i < numberOfNods; i++) {
        std::cout << std::setprecision(10) << Pbrzeg[i]<<" ";
        std::cout << std::setprecision(10) << vectorP[i];
        std::cout << "\n";
    }
 
    
}
void Grid::printTemperatures(int curentTime)
{
    double min = AMBIENT_TEMPERATURE;
    double max = 0.0;
    std::cout << std::setprecision(10);
    std::cout << "\nTemperatures after " << curentTime << " seconds \n";
    for (int i = 0; i < numberOfNods; i++)
    {
        if (marks[i].temp > max)
            max = marks[i].temp;
        if (marks[i].temp < min)
            min = marks[i].temp;
        std::cout << marks[i].temp << "\t";
        if ((i + 1) % 4 == 0)
           std::cout << "\n";
    }
    std::cout << "Min temp: " << min << "\tMax temp: " << max ;

}