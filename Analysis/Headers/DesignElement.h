#ifndef ELASTICITY_DESIGNELEMENT_H
#define ELASTICITY_DESIGNELEMENT_H
#include "Node.h"
#include "Element.h"

class DesignElement
{
public:
    DesignElement(std::vector<Element> ElementVector, double Stress);
    std::vector<Element> Elements;
    double Stress;
};

#endif //ELASTICITY_DESIHNELEMENT_H
