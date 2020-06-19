//
// Created by can on 26.05.2020.
//

#ifndef ELASTICITY_ELEMENT_H
#define ELASTICITY_ELEMENT_H
#include "Node.h"

class Element
{
public:
    Element();
    Element(int elmIndex, Node firstNode, Node secondNode, Node thirdNode, Node fourthNode, double elementMatrix[8][8]);
    int ElementIndex;
    Node FirstNode;
    Node SecondNode;
    Node ThirdNode;
    Node FourthNode;
    double ElementMatrix[8][8];
};

#endif //ELASTICITY_ELEMENT_H
