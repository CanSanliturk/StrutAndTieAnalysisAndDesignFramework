#include "Headers/DisplacementCalculator.h"
#include <iostream>
#include <vector>
#include <armadillo>

DisplacementCalculator::DisplacementCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> ForceVector, int nDof, int nDofRestrained)
{
	this->DisplacementVector = DisplacementVectorCalculator(GlobalStiffnessMatrix, ForceVector, nDof, nDofRestrained);
    this->SupportReactions = SupportReactionCalculator(GlobalStiffnessMatrix, this->DisplacementVector, nDof, nDofRestrained);
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

std::vector<double> DisplacementCalculator::SupportReactionCalculator(std::vector<std::vector<double>> GlobalStiffnessMatrix, std::vector<double> DisplacementVector, int nDof, int nDofRestrained)
{
    std::vector<double> supportReactions(nDofRestrained);
    int nDofUnrestrained = nDof - nDofRestrained;
    arma::mat kKU(nDofRestrained, nDofUnrestrained);
    arma::vec uU(nDofUnrestrained);
    kKU.fill(0);
    uU.fill(0);
    int rowCount = 0;

    for (int i = nDofUnrestrained; i < nDof; ++i)
    {
        for (int j = 0; j < nDofUnrestrained; ++j)
        {
            kKU(rowCount, j) = GlobalStiffnessMatrix.at(i).at(j);
        }
        rowCount++;
    }

    for (int i = 0; i < nDofUnrestrained; ++i)
    {
        uU(i) = DisplacementVector.at(i);
    }

    arma::mat tKKU = kKU.t();

    for (int i = 0; i < nDofRestrained; ++i)
    {
        double val = 0;

        for (int j = 0; j < nDofUnrestrained; ++j)
        {
            double valHelper = tKKU(j, i) * uU(j);
            val += valHelper;
        }

        supportReactions.at(i) = val;
    }

    return supportReactions;
}