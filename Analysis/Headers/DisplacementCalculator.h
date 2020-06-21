#ifndef ELASTICITY_DISPLACEMENTCALCULATOR_H
#define ELASTICITY_DISPLACEMENTCALCULATOR_H
#include <vector>
#include <armadillo>

class DisplacementCalculator
{
public:
    DisplacementCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> ForceVector, int nDof, int nDofRestrained);
    std::vector<double> DisplacementVector;
    std::vector<double> SupportReactions;
private:
	std::vector<double> DisplacementVectorCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> ForceVector, int nDof, int nDofRestrained);
	std::vector<double> SupportReactionCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> DispVector, int nDof, int nDofRestrained);
};



#endif //ELASTICITY_STIFFNESSMATRIXASSEMBLER_H