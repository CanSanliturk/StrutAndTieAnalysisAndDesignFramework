//
// Created by can on 27.05.2020.
//

#ifndef ELASTICITY_NATURALBC_H
#define ELASTICITY_NATURALBC_H


class NaturalBC
{
public:
    NaturalBC(double xStart, double xEnd, double yStart, double yEnd, double valueX, double valueY);
    double XStart;
    double XEnd;
    double YStart;
    double YEnd;
    double ValueX;
    double ValueY;
};


#endif //ELASTICITY_NATURALBC_H
