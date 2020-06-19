//
// Created by can on 27.05.2020.
//

#ifndef ELASTICITY_NODELISTFACTORY_H
#define ELASTICITY_NODELISTFACTORY_H
#include<vector>
#include"Gap.h"
#include"Node.h"
#include"EssentialBC.h"
#include"NaturalBC.h"

class NodeListFactory
{
public:

    NodeListFactory(std::vector<double> dimVector, std::vector<Gap> gapVector, std::vector<EssentialBC> ebcVector, std::vector<NaturalBC> nbcVector, double meshSize);
    std::vector<Node> NodeList;

private:
    std::vector<Node> GetNodeList(std::vector<double> dimVector, std::vector<Gap> gapVector, std::vector<EssentialBC> ebcVector, std::vector<NaturalBC> nbcVector, double meshSize);
    std::vector<double> _dimensionVector;
    std::vector<double> _gapVector;
    std::vector<EssentialBC> _ebcVector;
    std::vector<NaturalBC> _nbcVector;
    double _meshSize;
};


#endif //ELASTICITY_NODELISTFACTORY_H
