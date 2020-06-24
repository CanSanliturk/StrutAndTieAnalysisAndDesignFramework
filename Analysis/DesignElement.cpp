#include "Headers/DesignElement.h"
#include "Headers/Element.h"
#include "Headers/Node.h"

DesignElement::DesignElement(std::vector<Element> ElementVector, double Stress)
{
    this->Elements = ElementVector;
    this->Stress = Stress;
}