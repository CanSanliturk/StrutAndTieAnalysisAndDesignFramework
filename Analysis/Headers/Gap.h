//
// Created by can on 26.05.2020.
//

#ifndef ELASTICITY_GAP_H
#define ELASTICITY_GAP_H
class Gap
{
public:
    Gap(double xCoordInit, double xCoordEnd, double yCoordInit, double yCoordEnd);
    Gap();
    double XCoordInit;
    double XCoordEnd;
    double YCoordInit;
    double YCoordEnd;
};
#endif //ELASTICITY_GAP_H
