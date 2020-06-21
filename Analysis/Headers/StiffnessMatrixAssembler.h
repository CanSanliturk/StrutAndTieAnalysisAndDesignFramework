#ifndef ELASTICITY_STIFFNESSMATRIXASSEMBLER_H
#define ELASTICITY_STIFFNESSMATRIXASSEMBLER_H
#include <vector>
#include"Gap.h"
#include"Node.h"
#include "Element.h"
#include "EssentialBC.h"

class StiffnessMatrixAssembler
{
public:
    StiffnessMatrixAssembler(std::vector<Element> ElementVector, int nDof);
    std::vector<std::vector<double>> GetGlobalStiffnessMatrix(std::vector<Element> ElementVector, int nDof);
};



#endif //ELASTICITY_STIFFNESSMATRIXASSEMBLER_H