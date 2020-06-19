#include "Headers/ElementListFactory.h"
#include <vector>
#include "Headers/Node.h"
#include "Headers/Gap.h"
#include "Headers/Element.h"

ElementListFactory::ElementListFactory(std::vector<Node> nodeList, std::vector<Gap> gapVector, double h, double lX, double lY, double e, double pois, double thickness)
{
    this->ElementList = GetElementList(nodeList, gapVector, h, lX, lY, e, pois, thickness);
}

std::vector<Element> ElementListFactory::GetElementList(std::vector<Node> nodeList, std::vector<Gap> gapVector, double h, double lX, double lY, double e, double pois, double thickness)
{
    std::vector<Element> elmList;

    int numberOfElmXDir = lX / h;
    int numberOfElmYDir = lY / h;

    int numberOfNodeInXDir = numberOfElmXDir + 1;

    int elmIndex = 1;

    double mult =  thickness * e * h / (1 - (pois * pois));
    double v = pois;

    double kElm[8][8];

    kElm[0][0] = mult * (3 - v) / (6);
    kElm[0][1] = mult * (1 + v) / (8);
    kElm[0][2] = mult * (-3 - v) / (12);
    kElm[0][3] = mult * (-1 + 3 * v) / (8);
    kElm[0][4] = mult * (-3 + v) / (12);
    kElm[0][5] = mult * (-1 - v) / (8);
    kElm[0][6] = mult * (v) / (6);
    kElm[0][7] = mult * (1 - 3 * v) / (8);

    kElm[1][0] = kElm[0][1];
    kElm[1][1] = mult * (3 - v) / (6);
    kElm[1][2] = mult * (1 - 3 * v) / (8);
    kElm[1][3] = mult * (v) / (6);
    kElm[1][4] = mult * (-1 + v) / (8);
    kElm[1][5] = mult * (-3 + v) / (12);
    kElm[1][6] = mult * (-1 + 3 * v) / (8);
    kElm[1][7] = mult * (-3 + v) / (12);

    kElm[2][0] = kElm[0][2];
    kElm[2][1] = kElm[1][2];
    kElm[2][2] = mult * (3 - v) / (6);
    kElm[2][3] = mult * (-1 - v) / (8);
    kElm[2][4] = mult * (v) / (6);
    kElm[2][5] = mult * (-1 + 3 * v) / (8);
    kElm[2][6] = mult * (-3 + v) / (12);
    kElm[2][7] = mult * (1 + v) / (8);

    kElm[3][0] = kElm[0][3];
    kElm[3][1] = kElm[1][3];
    kElm[3][2] = kElm[2][3];
    kElm[3][3] = mult * (3 - v) / (6);
    kElm[3][4] = mult * (1 - 3 * v) / (8);
    kElm[3][5] = mult * (-3 - v) / (12);
    kElm[3][6] = mult * (1 + v) / (8);
    kElm[3][7] = mult * (-3 + v) / (12);

    kElm[4][0] = kElm[0][4];
    kElm[4][1] = kElm[1][4];
    kElm[4][2] = kElm[2][4];
    kElm[4][3] = kElm[3][4];
    kElm[4][4] = mult * (3 - v) / (6);
    kElm[4][5] = mult * (1 + v) / (8);
    kElm[4][6] = mult * (-3 - v) / (12);
    kElm[4][7] = mult * (-1 + 3 * v) / (8);

    kElm[5][0] = kElm[0][5];
    kElm[5][1] = kElm[1][5];
    kElm[5][2] = kElm[2][5];
    kElm[5][3] = kElm[3][5];
    kElm[5][4] = kElm[4][5];
    kElm[5][5] = mult * (3 - v) / (6);
    kElm[5][6] = mult * (1 - 3 * v) / (8);
    kElm[5][7] = mult * (v) / (6);

    kElm[6][0] = kElm[0][6];
    kElm[6][1] = kElm[1][6];
    kElm[6][2] = kElm[2][6];
    kElm[6][3] = kElm[3][6];
    kElm[6][4] = kElm[4][6];
    kElm[6][5] = kElm[5][6];
    kElm[6][6] = mult * (3 - v) / (6);
    kElm[6][7] = mult * (-1 - v) / (8);

    kElm[7][0] = kElm[0][7];
    kElm[7][1] = kElm[1][7];
    kElm[7][2] = kElm[2][7];
    kElm[7][3] = kElm[3][7];
    kElm[7][4] = kElm[4][7];
    kElm[7][5] = kElm[5][7];
    kElm[7][6] = kElm[6][7];
    kElm[7][7] = mult * (3 - v) / (6);


    for (int i = 0; i < numberOfElmYDir; i++)
    {
        for (int j = 0; j < numberOfElmXDir; j++)
        {
            // Specify nodes
            Node firstNode = nodeList[j + (i * numberOfNodeInXDir)];
            Node secondNode = nodeList[j + (i * numberOfNodeInXDir) + 1];
            Node thirdNode = nodeList[j + ((i + 1) * (numberOfNodeInXDir)) + 1];
            Node fourthNode = nodeList[j + ((i + 1) * (numberOfNodeInXDir))];

            // Check if any of node is in gap or not
            bool isFirstNodeInGap = (firstNode.XCoord == -5) && (firstNode.YCoord == -5);
            bool isSecondNodeInGap = (secondNode.XCoord == -5) && (secondNode.YCoord == -5);
            bool isThirdNodeInGap = (thirdNode.XCoord == -5) && (thirdNode.YCoord == -5);
            bool isFourthNodeInGap = (fourthNode.XCoord == -5) && (fourthNode.YCoord == -5);

            bool isAnyNodeInGap = isFirstNodeInGap || isSecondNodeInGap || isThirdNodeInGap || isFourthNodeInGap;

            if (isAnyNodeInGap)
            {
                continue;
            }

            // Check if element is in gap or not
            double centerX = (firstNode.XCoord + thirdNode.XCoord) / 2;
            double centerY = (firstNode.YCoord + thirdNode.YCoord) / 2;
            bool isInsideGap = false;

            for (int k = 0; k < gapVector.size(); k++)
            {
                Gap gap = gapVector[k];

                bool xCond = (centerX >= gap.XCoordInit) && (centerX <= gap.XCoordEnd);
                bool yCond = (centerY >= gap.YCoordInit) && (centerY <= gap.YCoordEnd);

                isInsideGap = xCond && yCond;

                if (isInsideGap)
                {
                    break;
                }

            }

            if (!isInsideGap)
            {
                // Create element and push it to element vector
                Element elm(elmIndex, firstNode, secondNode, thirdNode, fourthNode, kElm);
                elmList.push_back(elm);
                elmIndex++;
            }
        }
    }
    return elmList;
}
