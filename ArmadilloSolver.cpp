#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include "Analysis/Headers/Node.h"
#include "Analysis/Headers/Gap.h"
#include "Analysis/Headers/NaturalBC.h"
#include "Analysis/Headers/EssentialBC.h"
#include "Analysis/Headers/NodeListFactory.h"
#include "Analysis/Headers/ElementListFactory.h"
#include "Analysis/Headers/StiffnessMatrixAssembler.h"
#include "Analysis/Headers/ForceVectorAssembler.h"
#include "Analysis/Headers/DisplacementCalculator.h"
#include "Analysis/Headers/PrincipleStressCalculator.h"
#include <armadillo>

using namespace std;

bool IsEqual(double a, double b, double tol)
{
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
    double thickness = 1; // in meters

	// Control node to monitor displacement    
    double controlPointX = 3;
    double controlPointY = 0;
    bool isControlDisplacementInXDirection = false; // Name explains itself

    vector<double> dimVector{ lX, lY };
    vector<double> controlPointCoord{ controlPointX, controlPointY};


    // Material properties
    double e = 36000000000; // Elasticity modululus in Pa
    double v = 0.3; // Poisson's ratio
    double rho = 0; // Density of material in kg/m3 (if mass is not gonna be encountered, simply send it as "0")

    // Info of mesh
    double meshSize = 0.5; // in meters
    string meshType = "Quad";  // It is either "Quad" or "Triangular". Triangular mesh is not prepared yet (2020.05.21)

    // Info of gap(s)
    //Gap firstGap(2.2, 3.8, 0, 1.5);
 	Gap nullGap(-1, -1, -1, -1);
    vector<Gap> gapVector{ nullGap };

    // TODO : For now, if there is no gap, Gap nullGap(-1, -1, -1, -1). Define it and add to gapVector. But there is a need of improvement for that case.

    // Essential bc's on primary variable (Displacements in case of elasticity problem)
    //EssentialBC nullEBC(-1, -1, -1, -1, -1, -1);
    EssentialBC ebcFirst(0, 0, 0, 0, 1, 1);
    EssentialBC ebcSecond(6, 6, 0, 0, 1, 1);
    vector<EssentialBC> ebcVector{ ebcFirst, ebcSecond };

    //// Natural bc's on secondary variable (Forces in case of elasticity problem)
    /*NaturalBC firstNBC(0, 0, 0, 0, 0, 221667);
    NaturalBC secondNBC(6, 6, 0, 0, 0, 128333);*/
    NaturalBC firstNBC(3, 3, 3, 3, 0, -350000);
    vector<NaturalBC> nbcVector{ firstNBC };
    double tol = 0.001; // Set tolerance value to check equality

    /// SOLVER PART /// This part is going to be moved to a seperate class named as solver
    double controlDisplacement = 0;
    cout<<"Beginning of analysis"<<endl;
    //auto timenow =
    //        chrono::system_clock::to_time_t(chrono::system_clock::now());

    NodeListFactory nLF(dimVector, gapVector, ebcVector, nbcVector,  meshSize);
    vector<Node> nodeVec = nLF.NodeList;
    cout<<"Nodes are created"<<endl;

    ElementListFactory elmLF(nodeVec, gapVector, meshSize, lX, lY, e, v, thickness);
    vector<Element> elmVec = elmLF.ElementList;
    cout<<"Meshes are created"<<endl;

    int nDof = 0;
    int controlDof = 0;

    for (int i = 0; i < nodeVec.size(); ++i)
    {
        Node counterNode = nodeVec.at(i);
        int idxX = counterNode.DofIndexX;
        int idxY = counterNode.DofIndexY;

        if(idxX > 0)
        {
            nDof++; // Number of total degree of freedoms. (Nodes in gaps have dof index less than zero)
        }

        if(idxY > 0)
        {
            nDof++;
        }

        // Set degree of freedom that represents control displacement
    	if (IsEqual(counterNode.XCoord, controlPointX, tol) && IsEqual(counterNode.YCoord, controlPointY, tol))
    	{
			if (isControlDisplacementInXDirection)
			{
				controlDof = counterNode.DofIndexX;
			}
			else
			{
				controlDof = counterNode.DofIndexY;
			}
    	}
    }

    int nDofRestrained = 0;
    int nDofUnrestrained = 0;
    for (int i = 0; i < ebcVector.size(); ++i)
    {
        EssentialBC essentialBc = ebcVector.at(i);
        bool isPoint = IsEqual(essentialBc.XStart, essentialBc.XEnd, tol) && IsEqual(essentialBc.YStart, essentialBc.YEnd, tol);
        bool isXDir = (!IsEqual(essentialBc.XStart, essentialBc.XEnd, tol)) && IsEqual(essentialBc.YStart, essentialBc.YEnd, tol);
        bool isYDir = IsEqual(essentialBc.XStart, essentialBc.XEnd, tol) && (!IsEqual(essentialBc.YStart, essentialBc.YEnd, tol));

        int nRest = 1;
        if (essentialBc.ValueX > (-1 * tol) && essentialBc.ValueY > (-1 * tol))
        {
            nRest = 2;
        }

        int addVal = 0;

        if (isPoint)
        {
            addVal++;
        }
        else if(isXDir)
        {
            int restrainedNodeNumber = (essentialBc.XEnd - essentialBc.XStart) / meshSize;
            addVal += restrainedNodeNumber;
        }
        else if (isYDir)
        {
            int restrainedNodeNumber = (essentialBc.YEnd - essentialBc.XEnd) / meshSize;
            addVal += restrainedNodeNumber;
        }

        addVal *= nRest;

        nDofRestrained += addVal;

    }

    nDofUnrestrained = nDof - nDofRestrained;
    StiffnessMatrixAssembler sMA(elmVec, nDof);
    vector<vector<double>> kGlobal = sMA.GetGlobalStiffnessMatrix(elmVec, nDof);
    ForceVectorAssembler fVA(nodeVec, nbcVector, meshSize, nDof);
    std::vector<double> fGlobal = fVA.ForceVector;

    cout<<"Assembly of matrices and vectors are completed"<<endl;

    DisplacementCalculator dispCal(kGlobal, fGlobal, nDof, nDofRestrained);
    std::vector<double> dispVector = dispCal.DisplacementVector;
    std::vector<double> supportReactions = dispCal.SupportReactions;
    
    controlDisplacement = dispVector.at(controlDof - 1);

    ofstream StiffnessMatrixFile;
    ofstream ForceVectorFile;
    ofstream DisplacementFile;
    ofstream SupportReactionsFile;

    StiffnessMatrixFile.open("Outputs/AnalysisOutputs/StiffnessMatrix");
    ForceVectorFile.open("Outputs/AnalysisOutputs/ForceVectorFile");
    DisplacementFile.open("Outputs/AnalysisOutputs/DisplacementFile");
    SupportReactionsFile.open("Outputs/AnalysisOutputs/SupportReactions");

    for (int i = 0; i < nDof; ++i)
    {
        
        for (int j = 0; j < nDof; ++j)
        {            
            StiffnessMatrixFile << kGlobal.at(i).at(j);
            StiffnessMatrixFile << " ";
        }
        
        StiffnessMatrixFile << "\n";
        ForceVectorFile << fGlobal.at(i);
        ForceVectorFile << "\n";
        DisplacementFile << dispVector.at(i);
        DisplacementFile << "\n";

    }

    StiffnessMatrixFile.close();
    ForceVectorFile.close();
    DisplacementFile.close();
    
    for (int i = 0; i < nDofRestrained; ++i)
    {
        SupportReactionsFile << supportReactions.at(i);
        SupportReactionsFile << "\n";
    }
    
    SupportReactionsFile.close();

    cout<<"Displacements are calculated"<<endl;

    PrincipleStressCalculator pSC(elmVec, dispVector, e, v, meshSize);
    std::vector<std::vector<double>> principleStressVector(elmVec.size());
    principleStressVector = pSC.PrincipleStressList;
    int stressListSize = principleStressVector.size();
    
    cout<< "Principle stresses are calculated"<<endl;
    
	ofstream PrincipleStressFile;
	PrincipleStressFile.open("Outputs/AnalysisOutputs/PrincipleStresses");    

    for (int i = 0; i < stressListSize; ++i)
    {
    	std::vector<double> elmStress = principleStressVector.at(i);
    	double sigmaXX = elmStress.at(0) * 0.000001;
    	double sigmaYY = elmStress.at(1) * 0.000001;
    	double sigmaXY = elmStress.at(2) * 0.000001;
    	PrincipleStressFile <<"Element No: ";
    	PrincipleStressFile << (i + 1);
    	PrincipleStressFile <<"\n";
    	PrincipleStressFile <<"Stress in XX-Direction = ";
    	PrincipleStressFile <<sigmaXX;
    	PrincipleStressFile <<"\n";
    	PrincipleStressFile <<"Stress in YY-Direction = ";
    	PrincipleStressFile <<sigmaYY;
    	PrincipleStressFile <<"\n";
    	PrincipleStressFile <<"Stress in XY-Direction = ";;
    	PrincipleStressFile <<sigmaXY;
    	PrincipleStressFile <<"\n";
    }
    PrincipleStressFile.close();

    double sigmaMinAll = 0;
    double sigmaMaxAll = 0;

    std::vector<std::vector<double>> compressiveTensileStresses;

    ofstream MainStressFile;
    MainStressFile.open("Outputs/AnalysisOutputs/MainStressFile");

    for (int i = 0; i < elmVec.size(); ++i)
    {
    	Element elm = elmVec.at(i);
    	std::vector<double> elmStressVec = principleStressVector.at(i);
    	double sigmaXX = elmStressVec.at(0);
    	double sigmaYY = elmStressVec.at(1);
    	double sigmaXY = elmStressVec.at(2);

    	double lambdaOne = 0;
    	double lambdaTwo = 0;

    	double b = -sigmaXX - sigmaYY;
    	double c = (sigmaXX * sigmaYY) - (sigmaXY * sigmaXY);
    	double sqrtTerm = sqrt((b * b) - (4 * c));

    	double lambdaHelperOne = ((-1 * b) + (sqrtTerm)) / 2; 
    	double lambdaHelperTwo = ((-1 * b) - (sqrtTerm)) / 2;

    	double sigmaMin = 0;
    	double sigmaMax = 0;
    	
    	std::vector<double> elmSigma;

    	if (lambdaHelperOne > lambdaHelperTwo)
    	{
    		sigmaMin = lambdaHelperTwo;
    		sigmaMax = lambdaHelperOne;
    	}
    	else
    	{
    		sigmaMin = lambdaHelperOne;
    		sigmaMax = lambdaHelperTwo;
    	}

    	if (sigmaMin > 0)
    	{
    		sigmaMin = 0;
    	}
    	if (sigmaMax < 0)
    	{
    		sigmaMax = 0;
    	}

    	elmSigma.push_back(sigmaMin);
    	elmSigma.push_back(sigmaMax);

    	if (sigmaMin < sigmaMinAll)
    	{
    		sigmaMinAll = sigmaMin;
    	}
    	if (sigmaMax > sigmaMaxAll)
    	{
    		sigmaMaxAll = sigmaMax;
    	}
    	compressiveTensileStresses.push_back(vector<double>{ sigmaMin, sigmaMax });
    	MainStressFile << "Element Index ";
    	MainStressFile << i + 1;
    	MainStressFile << "\n";
    	MainStressFile << "Compressive Stress = ";
    	MainStressFile << sigmaMin * 0.000001;
    	MainStressFile << "MPa \n";
    	MainStressFile << "Tensile  Stress = ";
    	MainStressFile << sigmaMax * 0.000001;
    	MainStressFile << "MPa \n";
	}

    MainStressFile.close();


    // Element modification part
    std::vector<Element> ModifiedElmVec = elmVec;
    
    for (int i = 0; i < ModifiedElmVec.size(); ++i)
    {
    	Element elm = ModifiedElmVec.at(i);
    	std::vector<double> elmPrincipleStresses = compressiveTensileStresses.at(i);
    	
    	double elmCompStress = elmPrincipleStresses.at(0);
    	double elmTensStress = elmPrincipleStresses.at(1);

    	double alpha = elmCompStress / sigmaMinAll;
    	double beta = elmTensStress / sigmaMaxAll;

    	double modificationFactor = 1;

    	if (alpha > beta)
    	{
    		modificationFactor = alpha;
    	}
    	else
    	{
    		modificationFactor = beta;
    	}

		double mult =  thickness * e * meshSize / (1 - (v * v)) * modificationFactor;

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

    	Element modifiedElement(elm.ElementIndex, elm.FirstNode, elm.SecondNode, elm.ThirdNode, elm.FourthNode, kElm, 1);
    	ModifiedElmVec.at(i) = modifiedElement;
    }

    StiffnessMatrixAssembler modSMA(ModifiedElmVec, nDof);
    kGlobal = modSMA.GetGlobalStiffnessMatrix(ModifiedElmVec, nDof);
    ForceVectorAssembler modFVA(nodeVec, nbcVector, meshSize, nDof);
    fGlobal = modFVA.ForceVector;
    DisplacementCalculator modDispCal(kGlobal, fGlobal, nDof, nDofRestrained);
    dispVector = modDispCal.DisplacementVector;
    PrincipleStressCalculator modPSC(elmVec, dispVector, e, v, meshSize);
    principleStressVector = modPSC.PrincipleStressList;
    
    double monitoredDisp = dispVector.at(controlDof - 1);
    double displacementModificationFactor = monitoredDisp / controlDisplacement;



    //auto timenow2 =
    //        chrono::system_clock::to_time_t(chrono::system_clock::now());
    //cout<< "Elapsed Time = " << timenow2 - timenow << " seconds"<< endl;
    return 0;
}
