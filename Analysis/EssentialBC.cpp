//
// Created by can on 26.05.2020.
//
#include "Headers/EssentialBC.h"
EssentialBC::EssentialBC(double xStart, double xEnd, double yStart, double yEnd, double valueX, double valueY)
{
    this->XStart = xStart;
    this->XEnd = xEnd;
    this->YStart = yStart;
    this->YEnd = yEnd;
    this->ValueX = valueX;
    this->ValueY = valueY;

}
