//
// Created by can on 26.05.2020.
//

#ifndef ELASTICITY_ESSENTIALBC_H
#define ELASTICITY_ESSENTIALBC_H

class EssentialBC
{
public:
    EssentialBC(double xStart, double xEnd, double yStart, double yEnd, double valueX, double valueY);
    double XStart;
    double XEnd;
    double YStart;
    double YEnd;
    double ValueX;
    double ValueY;
};


#endif //ELASTICITY_ESSENTIALBC_H
