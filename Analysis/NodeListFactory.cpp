#include "Headers/NodeListFactory.h"
#include <vector>
#include "Headers/Node.h"
#include "Headers/Gap.h"
#include "Headers/EssentialBC.h"
#include "Headers/NaturalBC.h"


NodeListFactory::NodeListFactory(std::vector<double> dimVector, std::vector<Gap> gapVector, std::vector<EssentialBC> ebcVector, std::vector<NaturalBC> nbcVector, double meshSize)
{
    this->NodeList = GetNodeList(dimVector, gapVector, ebcVector, nbcVector, meshSize);
}

std::vector<Node> NodeListFactory::GetNodeList(std::vector<double> dimVector, std::vector<Gap> gapVector, std::vector<EssentialBC> ebcVector, std::vector<NaturalBC> nbcVector, double meshSize)
{
    std::vector<Node> nodeList;
    std::vector<Node> restrainedNodeList;

    double tol = 0.0000000001;

    double lX = dimVector[0];
    double lY = dimVector[1];
    double h = meshSize;

    int nNodeX = lX / h;
    int nNodeY = lY / h;

    double yCoord = 0;
    int nodeIdx = 1;
    int dofIdx = 1;

    for (int i = 0; i <= nNodeY; i++)
    {
        double xCoord = 0;

        for (int j = 0; j <= nNodeX; j++)
        {
            // Check if node is inside gap or not
            bool isInsideGap = false;

            int numOfGaps = gapVector.size();
            bool xCond = false;
            bool yCond = false;

            for (int k = 0; k < numOfGaps; k++)
            {
                Gap gap = gapVector[k];

                xCond = ((xCoord > gap.XCoordInit + tol) && (xCoord < gap.XCoordEnd - tol)) || (xCoord == 0 && gap.XCoordInit == 0) || (xCoord == lX && gap.XCoordEnd == lX);
                yCond = ((yCoord > gap.YCoordInit + tol) && (yCoord < gap.YCoordEnd - tol)) || (yCoord == 0 && gap.YCoordInit == 0) || (yCoord == lY && gap.YCoordEnd == lY);
                isInsideGap = (xCond && yCond);
                if (isInsideGap)
                {
                    break;
                }
            }

            if (!isInsideGap)
            {
                bool isRestrained = false;
                bool xCondForRestraint = false;
                bool yCondForRestraint = false;
                int numberOfEBC = ebcVector.size();

                // Check if node is restrained by a boundary conditions
                for (int l = 0; l < numberOfEBC; l++)
                {
                    EssentialBC ebc = ebcVector[l];

                    xCondForRestraint = (xCoord >= ebc.XStart - tol) && (xCoord <= ebc.XEnd + tol);
                    yCondForRestraint = (yCoord >= ebc.YStart - tol) && (yCoord <= ebc.YEnd + tol);
                    isRestrained = xCondForRestraint && yCondForRestraint;
                    if (isRestrained)
                    {
                        bool isRestrainedInXDir = ebc.ValueX != -5;
                        bool isRestrainedInYDir = ebc.ValueY != -5;

                        if (isRestrainedInXDir && isRestrainedInYDir)
                        {
                            Node node(xCoord, yCoord, nodeIdx, -1, -1);
                            nodeList.push_back(node);
                            nodeIdx++;
                        }
                        else if(isRestrainedInXDir && !isRestrainedInYDir)
                        {
                            Node node(xCoord, yCoord, nodeIdx, -1, dofIdx);
                            nodeList.push_back(node);
                            nodeIdx++;
                            dofIdx++;
                        }
                        else if (!isRestrainedInXDir && isRestrainedInYDir)
                        {
                            Node node(xCoord, yCoord, nodeIdx, dofIdx, -1);
                            nodeList.push_back(node);
                            nodeIdx++;
                            dofIdx++;
                        }
                        break;
                    }
                    else // Not restrained
                    {
                        if (l == numberOfEBC - 1) // In case of all of the conditions are checked and not a single restraining condition is found
                        {
                            Node node(xCoord, yCoord, nodeIdx, dofIdx, dofIdx + 1);
                            nodeList.push_back(node);
                            nodeIdx++;
                            dofIdx += 2;
                        }
                        else  // Not restrained for that BC but may be restrained by the others
                        {
                            continue;
                        }
                    }
                }

            }
            else
            {
                Node nullNode(-5, -5, -5, -5, -5);
                nodeList.push_back(nullNode);
            }

            xCoord += h;
        }
        yCoord += h;
    }

    for (int i = 0; i < nodeList.size(); i++)
    {
        Node node = nodeList[i];

        int dofIdxInX = node.DofIndexX;
        int dofIdxInY = node.DofIndexY;

        if (dofIdxInX == -1 || dofIdxInY == -1)
        {
            int newIdxX;
            int newIdxY;

            if (dofIdxInX == -1)
            {
                newIdxX = dofIdx;
                dofIdx++;
            }
            else
            {
                newIdxX = node.DofIndexX;
            }

            if (dofIdxInY == -1)
            {
                newIdxY = dofIdx;
                dofIdx++;
            }
            else
            {
                newIdxY = node.DofIndexY;
            }

            Node newNode(node.XCoord, node.YCoord, node.NodeIndex, newIdxX, newIdxY);
            nodeList[i] = newNode;
        }

    }

    return nodeList;
}

