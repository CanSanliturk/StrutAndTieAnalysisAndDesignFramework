#include "Headers/StiffnessMatrixAssembler.h"
#include <vector>
#include <armadillo>
#include "Headers/Node.h"
#include "Headers/Element.h"


StiffnessMatrixAssembler::StiffnessMatrixAssembler(std::vector<Element> ElementVector, int nDof)
{
}

std::vector<std::vector<double>> StiffnessMatrixAssembler::GetGlobalStiffnessMatrix(std::vector<Element> ElementVector, int nDof)
{
    std::vector<std::vector<double>> kGlobal;
    arma::mat asd(3, 3);
    std::vector<double> innerK(nDof);
    std::fill(innerK.begin(), innerK.end(), 0);

    for (int i = 0; i < nDof; ++i)
    {
        kGlobal.push_back(innerK);
    }

    for (int i = 0; i < ElementVector.size(); ++i)
    {
        Element elm = ElementVector.at(i);
        Node firstNode = elm.FirstNode;
        Node secondNode = elm.SecondNode;
        Node thirdNode = elm.ThirdNode;
        Node fourthNode = elm.FourthNode;

        int steeringVector[8];

        steeringVector[0] = firstNode.DofIndexX;
        steeringVector[1] = firstNode.DofIndexY;
        steeringVector[2] = secondNode.DofIndexX;
        steeringVector[3] = secondNode.DofIndexY;
        steeringVector[4] = thirdNode.DofIndexX;
        steeringVector[5] = thirdNode.DofIndexY;
        steeringVector[6] = fourthNode.DofIndexX;
        steeringVector[7] = fourthNode.DofIndexY;

        for (int j = 0; j < 8; ++j)
        {
            for (int k = 0; k < 8; ++k)
            {
                int firstIdx = steeringVector[j] - 1;
                int secondIdx = steeringVector[k] - 1;
                kGlobal.at(firstIdx).at(secondIdx) += elm.ElementMatrix[j][k];
            }
        }
    }


    return kGlobal;
}

