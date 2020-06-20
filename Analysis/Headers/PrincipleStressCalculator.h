#ifndef ELASTICITY_PRINCIPLESTRESSCALCULATOR_H
#define ELASTICITY_PRINCIPLESTRESSCALCULATOR_H

#include <vector>
#include "Element.h"
#include "Node.h"

class PrincipleStressCalculator
{
public:
    PrincipleStressCalculator(std::vector<Element> ElementList, std::vector<double> Displacements, double modulus, double poissonsRatio, double meshSize);
    std::vector<std::vector<double>> GetPrincipalStresses(std::vector<Element> ElementList, std::vector<double> Displacements, double modulus, double poissonsRatio, double meshSize);
    std::vector<std::vector<double>> PrincipleStressList; 
};
#endif 
