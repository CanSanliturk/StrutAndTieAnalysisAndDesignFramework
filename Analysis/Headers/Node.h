#ifndef ELASTICITY_NODE_H
#define ELASTICITY_NODE_H


class Node
{
public:
    Node(double xCoord, double yCoord, int nodeIndex, int dofIndexX, int dofIndexY);
    Node();
    double XCoord;
    double YCoord;
    int NodeIndex;
    int DofIndexX;
    int DofIndexY;

};


#endif //ELASTICITY_NODE_H
