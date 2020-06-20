#include <iostream>
#include <vector>
#include <chrono>
#include "Analysis/Headers/Node.h"
#include "Analysis/Headers/Gap.h"
#include "Analysis/Headers/NaturalBC.h"
#include "Analysis/Headers/EssentialBC.h"
#include "Analysis/Headers/NodeListFactory.h"
#include "Analysis/Headers/ElementListFactory.h"
#include "Analysis/Headers/PrincipleStressCalculator.h"
#include <armadillo>

using namespace std;

bool IsEqual(double a, double b)
{
    double tol = 0.0000000001;
    double upperBoundary = b + tol;
    double lowerBoundary = b - tol;
    bool isEqual = (a > lowerBoundary) && (a < upperBoundary);
    return isEqual;
}

int main()
{

    /// INPUT CARD ///
    // Length : m, Force : N
    // Surface outer dimensions
    double lX = 6; // in meters
    double lY = 3; // in meters
    double thickness = 0.25; // in meters
    vector<double> dimVector{ lX, lY };

    // Material properties
    double e = 33000000000; // Elasticity modululus in Pa
    double v = 0.3; // Poisson's ratio
    double rho = 2500; // Density of material in kg/m3 (if mass is not gonna be encountered, simply send it as "0")

    // Info of mesh
    double meshSize = 0.1; // in meters
    string meshType = "Quad";  // It is either "Quad" or "Triangular". Triangular mesh is not prepared yet (2020.05.21)

    // Info of gap(s)
    //Gap firstGap(2.2, 3.8, 0, 1.5);
 	Gap nullGap(-1, -1, -1, -1);
    vector<Gap> gapVector{ nullGap };

    // TODO : For now, if there is no gap, Gap nullGap(-1, -1, -1, -1). Define it and add to gapVector. But there is a need of improvement for that case.

    // Essential bc's on primary variable (Displacements in case of elasticity problem)
    //EssentialBC nullEBC(-1, -1, -1, -1, -1, -1);
    EssentialBC ebcFirst(0, 0, 0, 0, 1, 1);
    EssentialBC ebcSecond(0, 0, 3, 3, 1, 1);
    vector<EssentialBC> ebcVector{ ebcFirst, ebcSecond };

    //// Natural bc's on secondary variable (Forces in case of elasticity problem)
    /*NaturalBC firstNBC(0, 0, 0, 0, 0, 221667);
    NaturalBC secondNBC(6, 6, 0, 0, 0, 128333);*/
    NaturalBC firstNBC(6, 6, 3, 3, 0, -100000);
    NaturalBC secondNBC(6, 6, 0, 0, 0, -100000);
    vector<NaturalBC> nbcVector{ firstNBC };

    /// SOLVER PART /// This part is going to be moved to a seperate class named as solver

    cout<<"Beginning of analysis"<<endl;
    auto timenow =
            chrono::system_clock::to_time_t(chrono::system_clock::now());

    cout << ctime(&timenow) << endl;

    NodeListFactory nLF(dimVector, gapVector, ebcVector, nbcVector,  meshSize);
    vector<Node> nodeVec = nLF.NodeList;
    cout<<"Nodes are created"<<endl;

    ElementListFactory elmLF(nodeVec, gapVector, meshSize, lX, lY, e, v, thickness);
    vector<Element> elmVec = elmLF.ElementList;
    int numElm = elmVec.size();
    cout<<"Number of elements = "<<numElm<<endl;
    Element elmArr[numElm];
    for (int i = 0; i < numElm; ++i)
    {
        Element helperElm = elmVec.at(i);
        Element tempElm(helperElm.ElementIndex, helperElm.FirstNode, helperElm.SecondNode, helperElm.ThirdNode, helperElm.FourthNode, helperElm.ElementMatrix);
        elmArr[i] = tempElm;
    }
    cout<<"Meshes are created"<<endl;

    // Normally, i use StiffnessMatrixAssembler for the calcudouble elementMatrix[8][8]lations below but i get segmentation error in runtime
    // (not building). This is a temporary solution. I will work with pointers and arrays instead of vectors because
    // the error could not be solved. Really wonder why...
    int nDof = 0;
    for (int i = 0; i < nodeVec.size(); ++i)
    {
        Node counterNode = nodeVec.at(i);
        int idxX = counterNode.DofIndexX;
        int idxY = counterNode.DofIndexY;

        if(idxX > 0)
        {
            nDof++;
        }

        if(idxY > 0)
        {
            nDof++;
        }
    }

    int nDofRestrained = 0;
    int nDofUnrestrained = 0;
    for (int i = 0; i < ebcVector.size(); ++i)
    {
        EssentialBC essentialBc = ebcVector.at(i);
        bool isPoint = IsEqual(essentialBc.XStart, essentialBc.XEnd) && IsEqual(essentialBc.YStart, essentialBc.YEnd);
        bool isXDir = (!IsEqual(essentialBc.XStart, essentialBc.XEnd)) && IsEqual(essentialBc.YStart, essentialBc.YEnd);
        bool isYDir = IsEqual(essentialBc.XStart, essentialBc.XEnd) && (!IsEqual(essentialBc.YStart, essentialBc.YEnd));

        if (isPoint)
        {
            nDofRestrained++;
        }
        else if(isXDir)
        {
            int restrainedNodeNumber = (essentialBc.XEnd - essentialBc.XStart) / meshSize;
            nDofRestrained += restrainedNodeNumber;
        }
        else if (isYDir)
        {
            int restrainedNodeNumber = (essentialBc.YEnd - essentialBc.XEnd) / meshSize;
            nDofRestrained += restrainedNodeNumber;
        }
    }

    nDofUnrestrained = nDof - nDofRestrained;

    vector<vector<double>> kGlobal;
    vector<double> innerK(nDof);
    std::fill(innerK.begin(), innerK.end(), 0);

    for (int i = 0; i < nDof; ++i)
    {
        kGlobal.push_back(innerK);
    }

    for (int i = 0; i < numElm; ++i)
    {
        Element elm = elmArr[i];
        Node firstNode = elm.FirstNode;
        Node secondNode = elm.SecondNode;
        Node thirdNode = elm.ThirdNode;
        Node fourthNode = elm.FourthNode;

        int steeringVector[8];

        steeringVector[0] = firstNode.DofIndexX;
        steeringVector[1] = firstNode.DofIndexY;
        steeringVector[2] = secondNode.DofIndexX;
        steeringVector[3] = secondNode.DofIndexY;
        steeringVector[4] = thirdNode.DofIndexX;
        steeringVector[5] = thirdNode.DofIndexY;
        steeringVector[6] = fourthNode.DofIndexX;
        steeringVector[7] = fourthNode.DofIndexY;

        for (int j = 0; j < 8; ++j)
        {
            for (int k = 0; k < 8; ++k)
            {
                int firstIdx = steeringVector[j] - 1;
                int secondIdx = steeringVector[k] - 1;
                kGlobal.at(firstIdx).at(secondIdx) += elm.ElementMatrix[j][k];
            }
        }
    }

    // Same story of stiffness matrix is valid for force vector, too...
    // For now, distribute distributed loads equally to adjacent nodes. It will not cause the program to work poorly unless
    // mesh size is very large. If mesh size is very large, well it should not be any way :D.

    vector<double> forceVector(nDof);
    fill(forceVector.begin(), forceVector.end(), 0);

    for (auto nbc : nbcVector)
    {
        double xSt = nbc.XStart;
        double xEnd = nbc.XEnd;
        double ySt = nbc.YStart;
        double yEnd = nbc.YEnd;
        double xVal = nbc.ValueX;
        double yVal = nbc.ValueY;

        // Check the type of the load(whether it is point load, distributed load in x-dir or distributed load in y-dir)
        bool isPtLoad = IsEqual(xSt, xEnd) && IsEqual(ySt, yEnd);
        bool isDistXDir = IsEqual(ySt, yEnd) && (!IsEqual(xSt, xEnd));
        bool isDistYDir = IsEqual(xSt, xEnd) && (!IsEqual(ySt, yEnd));

        if (isPtLoad)
        {
            for (int j = 0; j < nodeVec.size(); ++j)
            {
                Node loadedNode = nodeVec.at(j);
                double nodalX = loadedNode.XCoord;
                double nodalY = loadedNode.YCoord;
                double dofIdxXDir = loadedNode.DofIndexX;
                double dofIdxYDir = loadedNode.DofIndexY;

                if (IsEqual(xSt, nodalX) && IsEqual(ySt, nodalY))
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

                if (IsEqual(nodalX, xSt))
                {
                    double loadToBeAdded = yVal * meshSize / 2;
                    forceVector.at(dofIdxYDir) += loadToBeAdded;
                }
                else if (IsEqual(nodalX, xSt))
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

                if (IsEqual(nodalY, ySt))
                {
                    double loadToBeAdded = xVal * meshSize / 2;
                    forceVector.at(dofIdxXDir) += loadToBeAdded;
                }
                else if (IsEqual(nodalY, ySt))
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

    cout<<"Assembly of matrices and vectors are completed"<<endl;

    arma::mat kArma(nDofUnrestrained, nDofUnrestrained);
    arma::mat kGlobalArma(nDof, nDof);
    arma::vec fGlobalArma(nDof);
    arma::vec fArma(nDofUnrestrained);
    kArma.fill(0);
    kGlobalArma.fill(0);	
    fArma.fill(0);
    fGlobalArma.fill(0);

    for (int i = 0; i < nDofUnrestrained; ++i)
    {
        for (int j = 0; j < nDofUnrestrained; ++j)
        {
            kArma(i, j) = kGlobal.at(i).at(j);
        }
        fArma(i) = forceVector.at(i);
    }

    for (int i = 0; i < nDof; ++i)
    {
    	for (int j = 0; j < nDof; ++j)
    	{
    		kGlobalArma(i, j) = kGlobal.at(i).at(j);
    	}
    }

    cout<<"Beginning of solver"<<endl;
    
    arma::vec resDataHelper(nDofUnrestrained);
    resDataHelper = arma::solve(kArma, fArma);
    arma::vec resData(nDof);
    std::vector<double> dispVector;

    for (int i = 0; i < nDof; ++i)
    {
    	if (i < nDofUnrestrained)
    	{
    		resData(i) = resDataHelper(i);
    		dispVector.push_back(resDataHelper(i));
    	}
    	else
    	{
    		resData(i) = 0;
    		dispVector.push_back(0);
    	}
    }

    cout<<"End of solver"<<endl;
    cout<<"Displacements are calculated"<<endl;

    PrincipleStressCalculator pSC(elmVec, dispVector, e, v, meshSize);
    cout<<"Chk1"<<endl;
    std::vector<std::vector<double>> principleStressVector(elmVec.size());
    cout<<"Chk2"<<endl;
    principleStressVector = pSC.PrincipleStressList;
    cout<<"Chk3"<<endl;
    int stressListSize = principleStressVector.size();
    cout<<"stressListSize="<<stressListSize<<endl;

    for (int i = 0; i < stressListSize; ++i)
    {
    	std::vector<double> elmStress = principleStressVector.at(i);
    	double sigmaXX = elmStress.at(0) * 0.000001;
    	double sigmaYY = elmStress.at(1) * 0.000001;
    	double sigmaXY = elmStress.at(2) * 0.000001;

    	cout<<"--------------------------------------------------------------"<<endl;
    	cout<<"Element No: "<<i + 1<<endl;
    	cout<<"Stress in XX-Direction = "<<sigmaXX<<" MPa"<<endl;
    	cout<<"Stress in YY-Direction = "<<sigmaYY<<" MPa"<<endl;
    	cout<<"Stress in XY-Direction = "<<sigmaXY<<" MPa"<<endl;
    }


    auto timenow2 =
            chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << ctime(&timenow2) << endl;
    cout<< "Elapsed Time = " << timenow2 - timenow << " seconds"<< endl;

    cout<<kGlobal.size()<<endl;

    return 0;
}

