//
// Created by can on 26.05.2020.
//

#ifndef ELASTICITY_ELEMENTLISTFACTORY_H
#define ELASTICITY_ELEMENTLISTFACTORY_H
#include <vector>
#include"Gap.h"
#include"Node.h"
#include "Element.h"

class ElementListFactory
{
public:
    ElementListFactory(std::vector<Node> nodeList, std::vector<Gap> gapVector, double h, double lX, double lY, double e, double pois, double thickness);
    std::vector<Element> ElementList;
private:
    std::vector<Element> GetElementList(std::vector<Node> nodeList, std::vector<Gap> gapVector, double h, double lX, double lY, double e, double pois, double thickness);
};



#endif //ELASTICITY_ELEMENTLISTFACTORY_H
