#include "Headers/Element.h"
#include "Headers/Node.h"



Element::Element(int elmIndex, Node firstNode, Node secondNode, Node thirdNode, Node fourthNode, double elementMatrix[8][8])
{
    this->ElementIndex = elmIndex;
    this->FirstNode = firstNode;
    this->SecondNode = secondNode;
    this->ThirdNode = thirdNode;
    this->FourthNode = fourthNode;
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            this->ElementMatrix[i][j] = elementMatrix[i][j];
        }
    }
}

Element::Element()
{
    this->ElementIndex = -1;
    this->FirstNode = Node();
    this->SecondNode = Node();
    this->ThirdNode = Node();
    this->FourthNode = Node();
}
