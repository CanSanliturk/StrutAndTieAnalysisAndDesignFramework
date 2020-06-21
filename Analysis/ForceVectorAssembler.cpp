#include "Headers/ForceVectorAssembler.h"
#include <vector>
#include "Headers/Node.h"
#include "Headers/NaturalBC.h"


ForceVectorAssembler::ForceVectorAssembler(std::vector<Node> NodeVector, std::vector<NaturalBC> NaturalBCVector, double meshSize, int nDof)
{
    this->ForceVector = GetForceVector(NodeVector, NaturalBCVector, meshSize, nDof);
}

std::vector<double> ForceVectorAssembler::GetForceVector(std::vector<Node> nodeVec, std::vector<NaturalBC> NaturalBCVector, double meshSize, int nDof)
{
    std::vector<double> forceVector(nDof);
    std::fill(forceVector.begin(), forceVector.end(), 0);
    double tol = 0.000001; // Set tolerance value to check equality

    for (auto nbc : NaturalBCVector)
    {
        double xSt = nbc.XStart;
        double xEnd = nbc.XEnd;
        double ySt = nbc.YStart;
        double yEnd = nbc.YEnd;
        double xVal = nbc.ValueX;
        double yVal = nbc.ValueY;

        // Check the type of the load(whether it is point load, distributed load in x-dir or distributed load in y-dir)
        bool isPtLoad = IsEqual(xSt, xEnd, tol) && IsEqual(ySt, yEnd, tol);
        bool isDistXDir = IsEqual(ySt, yEnd, tol) && (!IsEqual(xSt, xEnd, tol));
        bool isDistYDir = IsEqual(xSt, xEnd, tol) && (!IsEqual(ySt, yEnd, tol));

        if (isPtLoad)
        {
            for (int j = 0; j < nodeVec.size(); ++j)
            {
                Node loadedNode = nodeVec.at(j);
                double nodalX = loadedNode.XCoord;
                double nodalY = loadedNode.YCoord;
                double dofIdxXDir = loadedNode.DofIndexX;
                double dofIdxYDir = loadedNode.DofIndexY;

                if (IsEqual(xSt, nodalX, tol) && IsEqual(ySt, nodalY, tol))
                {
                    forceVector.at(dofIdxXDir) += xVal;
                    forceVector.at(dofIdxYDir) += yVal;
                    break;
                }
            }
        }

        else if (isDistXDir) // if it is distributed in x-dir, then it is loaded in y-dir
        {
            for (int j = 0; j < nodeVec.size(); ++j)
            {
                Node loadedNode = nodeVec.at(j);
                double nodalX = loadedNode.XCoord;
                double nodalY = loadedNode.YCoord;
                double dofIdxYDir = loadedNode.DofIndexY;

                if (IsEqual(nodalX, xSt, tol))
                {
                    double loadToBeAdded = yVal * meshSize / 2;
                    forceVector.at(dofIdxYDir) += loadToBeAdded;
                }
                else if (IsEqual(nodalX, xSt, tol))
                {
                    double loadToBeAdded = yVal * meshSize / 2;
                    forceVector.at(dofIdxYDir) += loadToBeAdded;
                }
                else if ((nodalX > xSt) && (nodalX < xEnd))
                {
                    double loadToBeAdded = yVal * meshSize;
                    forceVector.at(dofIdxYDir) += loadToBeAdded;
                }
            }
        }
        else if (isDistYDir) // if it is distributed in y-dir, then it is loaded in x-dir
        {
            for (int j = 0; j < nodeVec.size(); ++j)
            {
                Node loadedNode = nodeVec.at(j);
                double nodalX = loadedNode.XCoord;
                double nodalY = loadedNode.YCoord;
                double dofIdxXDir = loadedNode.DofIndexX;

                if (IsEqual(nodalY, ySt, tol))
                {
                    double loadToBeAdded = xVal * meshSize / 2;
                    forceVector.at(dofIdxXDir) += loadToBeAdded;
                }
                else if (IsEqual(nodalY, ySt, tol))
                {
                    double loadToBeAdded = xVal * meshSize / 2;
                    forceVector.at(dofIdxXDir) += loadToBeAdded;
                }
                else if ((nodalX > xSt) && (nodalX < xEnd))
                {
                    double loadToBeAdded = yVal * meshSize;
                    forceVector.at(dofIdxXDir) += loadToBeAdded;
                }
            }
        }
    }


    return forceVector;
}

bool ForceVectorAssembler::IsEqual(double a, double b, double tol)
{
    double upperBoundary = b + tol;
    double lowerBoundary = b - tol;
    bool isEqual = (a > lowerBoundary) && (a < upperBoundary);
    return isEqual;
}