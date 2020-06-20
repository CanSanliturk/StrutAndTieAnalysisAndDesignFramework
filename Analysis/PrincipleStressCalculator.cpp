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
    	double kModifier = elm.StiffnessModifier;

	    arma::mat elasticityMat(3, 3);
	    elasticityMat.fill(0);
	    double mult = (modulus * kModifier) / ((1 + poissonsRatio) * (1 - poissonsRatio));
	    elasticityMat(0, 0) = mult * (1 - poissonsRatio);
	    elasticityMat(0, 1) = mult * poissonsRatio;
	    elasticityMat(1, 0) = mult * poissonsRatio;
	    elasticityMat(1, 1) = mult * (1 - poissonsRatio);
	    elasticityMat(2, 2) = mult * 0.5 * (1 - (2 * poissonsRatio));

    }

	arma::vec asdasd(3);

    return principalStressList;
}
