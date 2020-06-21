#include "Headers/DisplacementCalculator.h"
#include <iostream>
#include <vector>
#include <armadillo>

DisplacementCalculator::DisplacementCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> ForceVector, int nDof, int nDofRestrained)
{
	this->DisplacementVector = DisplacementVectorCalculator(GlobalStiffnessMatrix, ForceVector, nDof, nDofRestrained);
}

std::vector<double> DisplacementCalculator::DisplacementVectorCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> ForceVector, int nDof, int nDofRestrained)
{
	std::vector<double> dispVector;
	int nDofUnrestrained = nDof - nDofRestrained;

	arma::mat kArma(nDofUnrestrained, nDofUnrestrained);
    arma::vec fArma(nDofUnrestrained);
    kArma.fill(0);
    fArma.fill(0);

    for (int i = 0; i < nDofUnrestrained; ++i)
    {
        for (int j = 0; j < nDofUnrestrained; ++j)
        {
            kArma(i, j) = GlobalStiffnessMatrix.at(i).at(j);
        }
        fArma(i) = ForceVector.at(i);
    }

    arma::vec resDataHelper(nDofUnrestrained);
    resDataHelper = arma::solve(kArma, fArma);

    for (int i = 0; i < nDof; ++i)
    {
    	if (i < nDofUnrestrained)
    	{
    		dispVector.push_back(resDataHelper(i));
    	}
    	else
    	{
    		dispVector.push_back(0);
    	}
    }

	return dispVector;
}
