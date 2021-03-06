#ifndef ELASTICITY_STIFFNESSMATRIXASSEMBLER_H
#define ELASTICITY_STIFFNESSMATRIXASSEMBLER_H
#include <vector>
#include"Node.h"
#include "Element.h"

class StiffnessMatrixAssembler
{
public:
    StiffnessMatrixAssembler(std::vector<Element> ElementVector, int nDof);
    std::vector<std::vector<double>> GetGlobalStiffnessMatrix(std::vector<Element> ElementVector, int nDof);
};



#endif //ELASTICITY_STIFFNESSMATRIXASSEMBLER_H