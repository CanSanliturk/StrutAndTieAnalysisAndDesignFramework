#include "Headers/PrincipleStressCalculator.h"
#include <vector>
#include "Headers/Node.h"
#include "Headers/Element.h"
#include <armadillo>

PrincipleStressCalculator::PrincipleStressCalculator(std::vector<Element> ElementList, std::vector<double> Displacements, double modulus, double poissonsRatio, double meshSize)
{
	this->PrincipleStressList = GetPrincipalStresses(ElementList, Displacements, modulus, poissonsRatio, meshSize);
}

std::vector<std::vector<double>> PrincipleStressCalculator::GetPrincipalStresses(std::vector<Element> ElementList, std::vector<double> Displacements, double modulus, double poissonsRatio, double meshSize)
{
    std::vector<std::vector<double>> principalStressList;

    int nElm = ElementList.size();

    for (int i = 0; i < nElm; ++i)
    {
     	Element elm = ElementList.at(i);
    	Node firstNode = elm.FirstNode;
    	Node secondNode = elm.SecondNode;
    	Node thirdNode = elm.ThirdNode;
    	Node fourthNode = elm.FourthNode;

    	std::vector<Node> nodeVec;
    	nodeVec.push_back(firstNode);
    	nodeVec.push_back(secondNode);
    	nodeVec.push_back(thirdNode);
    	nodeVec.push_back(fourthNode);

    	double kModifier = elm.StiffnessModifier;

    	double midPtX = firstNode.XCoord + (meshSize / 2);
    	double midPtY = firstNode.YCoord + (meshSize / 2);
    	double h = meshSize;
    	double hSq = h * h;

	    arma::mat elasticityMat(3, 3);
	    elasticityMat.fill(0);
	    double mult = (modulus * kModifier) / ((1 + poissonsRatio) * (1 - poissonsRatio));
	    elasticityMat(0, 0) = mult * (1 - poissonsRatio);
	    elasticityMat(0, 1) = mult * poissonsRatio;
	    elasticityMat(1, 0) = mult * poissonsRatio;
	    elasticityMat(1, 1) = mult * (1 - poissonsRatio);
	    elasticityMat(2, 2) = mult * 0.5 * (1 - (2 * poissonsRatio));

	    arma::mat bMat(3, 8);
	    bMat.fill(0);
	    bMat(0, 0) = (midPtY - h) / hSq;
	    bMat(0, 2) = -1 * bMat(0, 0);
	    bMat(0, 4) = midPtY / hSq;
	    bMat(0, 6) = -1 * bMat(0, 4);
	    bMat(1, 1) = (midPtX - h) / hSq;
	    bMat(1, 3) = -1 * bMat(1, 1);
	    bMat(1, 5) = midPtX / hSq;
	    bMat(1, 7) =  -1 * bMat(1, 5);
    	bMat(2, 0) = bMat(1, 1);
    	bMat(2, 1) = bMat(0, 0);
    	bMat(2, 2) = bMat(1, 3);
    	bMat(2, 3) = bMat(0, 2);
    	bMat(2, 4) = bMat(1, 5);
    	bMat(2, 5) = bMat(0, 4);
    	bMat(2, 6) = bMat(1, 7);
    	bMat(2, 7) = bMat(0, 6);

    	arma::vec dispVec(8);
    	int counter = 0;
    	for (int j = 0; j < 4; ++j)
    	{
    		Node elmNode = nodeVec.at(j);
    		int dofX = elmNode.DofIndexX;
    		int dofY = elmNode.DofIndexY;
    		dispVec(counter) = (Displacements.at(dofX - 1));
    		counter++;
    		dispVec(counter) = (Displacements.at(dofY - 1));
    		counter++;
    	}

    	arma::vec elmMidPtStressArma(3);
    	arma::mat helperMat(3, 8);
    	helperMat = elasticityMat * bMat;
    	elmMidPtStressArma = helperMat * dispVec;
    	std::vector<double> elmMidPtStress;
    	for (int j = 0; j < 3; ++j)
    	{
    		elmMidPtStress.push_back(elmMidPtStressArma(j));
    	}

    	principalStressList.push_back(elmMidPtStress);
    }

    return principalStressList;
}
