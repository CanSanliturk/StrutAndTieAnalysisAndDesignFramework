#include "Headers/Gap.h"
Gap::Gap()
{
    this->XCoordInit = 0;
    this->XCoordEnd = 0;
    this->YCoordInit = 0;
    this->YCoordEnd = 0;
}
Gap::Gap(double xCoordInit, double xCoordEnd, double yCoordInit, double yCoordEnd)
{
    this->XCoordInit = xCoordInit;
    this->XCoordEnd = xCoordEnd;
    this->YCoordInit = yCoordInit;
    this->YCoordEnd = yCoordEnd;
}
