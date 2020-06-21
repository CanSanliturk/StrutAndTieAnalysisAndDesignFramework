#ifndef ELASTICITY_FORCEVECTORASSEMBLER_H
#define ELASTICITY_FORCEVECTORASSEMBLER_H
#include <vector>
#include"Node.h"
#include "NaturalBC.h"

class ForceVectorAssembler
{
public:
    ForceVectorAssembler(std::vector<Node> NodeVector, std::vector<NaturalBC> NaturalBCVector, double meshSize, int nDof);
    std::vector<double> ForceVector;
private:
    std::vector<double> GetForceVector(std::vector<Node> NodeVector, std::vector<NaturalBC> NaturalBCVector, double meshSize, int nDof);
	bool IsEqual(double a, double b, double tol);
};



#endif //ELASTICITY_FORCEVECTORASSEMBLERASSEMBLER_H