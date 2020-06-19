#include "Headers/Node.h"

Node::Node() // Node at gap
{
    this->NodeIndex = -5;
    this->XCoord = -5;
    this->YCoord = -5;
    this->DofIndexX = -5;
    this->DofIndexY = -5;
}

Node::Node(double xCoord, double yCoord, int nodeIndex, int dofIndexX, int dofIndexY)
{
    this->XCoord = xCoord;
    this->YCoord = yCoord;
    this->NodeIndex = nodeIndex;
    this->DofIndexX = dofIndexX;
    this->DofIndexY = dofIndexY;
}


