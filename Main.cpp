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

vector<vector<double>> MainStressCalculator(vector<Element> elmVec, vector<vector<double>> principleStressVector)
{
    std::vector<std::vector<double>> compressiveTensileStresses;

    double sigmaMinAll = 0;
    double sigmaMaxAll = 0;
    double sigmaMinAverage = 0;
    double sigmaMaxAverage = 0;

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

    	sigmaMinAverage += sigmaMin;
    	sigmaMaxAverage += sigmaMax;

	}

	sigmaMinAverage = sigmaMinAverage / elmVec.size();
	sigmaMaxAverage = sigmaMaxAverage / elmVec.size();

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

    	if (sigmaMin < sigmaMinAll && sigmaMin > (50 * sigmaMinAverage)) // To avoid stress localization effect
    	{
    		sigmaMinAll = sigmaMin;
    	}
    	if (sigmaMax > sigmaMaxAll && sigmaMax < (50 * sigmaMaxAverage))
    	{
    		sigmaMaxAll = sigmaMax;
    	}
    	compressiveTensileStresses.push_back(vector<double>{ sigmaMin, sigmaMax });
	}

    return compressiveTensileStresses;
}

vector<double> SigmaMinMaxCalculator(vector<Element> elmVec, vector<vector<double>> principleStressVector)
{
	vector<double> sigmaMinMax;
    double sigmaMinAll = 0;
    double sigmaMaxAll = 0;

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
	}
	sigmaMinMax.push_back(sigmaMinAll);
	sigmaMinMax.push_back(sigmaMaxAll);
    return sigmaMinMax;
}

vector<Element> ElementModificator(vector<Element> elmVec, vector<double> dispVector, double e, double v, double meshSize, double thickness, vector<double> fGlobal, int nDof, int nDofRestrained, vector<Node> nodeVec, int controlDof, double controlDisplacement, vector<NaturalBC> nbcVector)
{
	PrincipleStressCalculator pSC(elmVec, dispVector, e, v, meshSize);
    std::vector<std::vector<double>> principleStressVector(elmVec.size());
    principleStressVector = pSC.PrincipleStressList;
    
    vector<vector<double>> compressiveTensileStresses = MainStressCalculator(elmVec, principleStressVector);

    // Element modification part
    vector<Element> ModifiedElmVec = elmVec;   
    vector<double> sigmaMinMax = SigmaMinMaxCalculator(elmVec, principleStressVector);

    double sigmaMinAll = sigmaMinMax.at(0);
    double sigmaMaxAll = sigmaMinMax.at(1);

    for (int i = 0; i < elmVec.size(); ++i)
    {
    	Element elm = elmVec.at(i);
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

    	double kElm[8][8] = { { 0 } };

    	for (int j = 0; j < 8; ++j)
    	{
    		for (int k = 0; k < 8; ++k)
    		{
    			kElm[j][k] = elm.ElementMatrix[j][k] * modificationFactor;
    		}
    	}

    	Element modifiedElement(elm.ElementIndex, elm.FirstNode, elm.SecondNode, elm.ThirdNode, elm.FourthNode, kElm, modificationFactor);
    	ModifiedElmVec.at(i) = modifiedElement;
    }

    StiffnessMatrixAssembler modSMA(ModifiedElmVec, nDof);
    std::vector<std::vector<double>> kGlobal = modSMA.GetGlobalStiffnessMatrix(ModifiedElmVec, nDof);
    ForceVectorAssembler modFVA(nodeVec, nbcVector, meshSize, nDof);
    fGlobal = modFVA.ForceVector;
    DisplacementCalculator modDispCal(kGlobal, fGlobal, nDof, nDofRestrained);
    dispVector = modDispCal.DisplacementVector;
    
    // Multiply elasticity matrix with displacementModificationFactor
    double monitoredDisp = dispVector.at(controlDof - 1);
    double displacementModificationFactor = monitoredDisp / controlDisplacement;

    vector<Element> newModifiedElementVec;

    for (int i = 0; i < ModifiedElmVec.size(); ++i)
    {
    	Element elm = ModifiedElmVec.at(i);

    	double kElm[8][8] = { { 0 } };

    	for (int j = 0; j < 8; ++j)
    	{
    		for (int k = 0; k < 8; ++k)
    		{
    			kElm[j][k] = elm.ElementMatrix[j][k] * displacementModificationFactor;
    		}
    	}

    	Element modifiedElement(elm.ElementIndex, elm.FirstNode, elm.SecondNode, elm.ThirdNode, elm.FourthNode, kElm, elm.StiffnessModifier * displacementModificationFactor);
    	newModifiedElementVec.push_back(modifiedElement);
    }
	
    return newModifiedElementVec;
}

double CalculateShearReinforcementRatio(double dia, double nRow, double spacing, double thickness)
{
    double rho = 0;
    rho = (nRow * 3.141593 * dia * dia * 0.25) / (spacing * thickness);
    return rho;
}

int main()
{
    /// INPUT CARD ///
    // Length : m, Force : N
    // Surface outer dimensions
    double lX = 6; // in meters
    double lY = 3; // in meters
    double thickness = 0.5; // in meters

	// Control node to monitor displacement    
    double controlPointX = 3;
    double controlPointY = 0;
    bool isControlDisplacementInXDirection = false; // Name explains itself

    vector<double> dimVector{ lX, lY };
    vector<double> controlPointCoord{ controlPointX, controlPointY};

    // Material properties
    double e = 36000000000; // Elasticity modululus in Pa
    double v = 0.0; // Poisson's ratio
    double rho = 0.3; // Density of material in kg/m3 (if mass is not gonna be encountered, simply send it as "0")
    double fC = 20000000; // Compressive strength of concrete in Pa. Pay attention that it is the design strength
    double fY = 420000000; // Yield strength of steel in Pa. Pay attention that it is the design strength
    double longitudinalBarDiameter = 0.014; // Tension bars diameter in meters
    double transverseBarDiameter = 0.012; // Transverse bars diameters in meters

    // Info of mesh
    double meshSize = 0.5; // in meters
    string meshType = "Quad";  // It is either "Quad" or "Triangular". Triangular mesh is not prepared yet (2020.05.21)

    // Info of gap(s)
    //Gap firstGap(2.2, 3.8, 0, 1.5);
 	Gap nullGap(-1, -1, -1, -1);
    vector<Gap> gapVector{ nullGap };

    // Essential bc's on primary variable (Displacements in case of elasticity problem)
    //EssentialBC nullEBC(-1, -1, -1, -1, -1, -1);
    EssentialBC ebcFirst(0, 0, 0, 0, 1, 1);
    EssentialBC ebcSecond(6, 6, 0, 0, 1, 1);
    vector<EssentialBC> ebcVector{ ebcFirst, ebcSecond };

    //// Natural bc's on secondary variable (Forces in case of elasticity problem)
    NaturalBC firstNBC(3, 3, 3, 3, 0, -1400000);
    vector<NaturalBC> nbcVector{ firstNBC };
    double tol = 0.001; // Set tolerance value to check equality

    /// SOLVER PART /// 
    double controlDisplacement = 0;
    cout<<"Beginning of analysis"<<endl;
    auto timenow =
            chrono::system_clock::to_time_t(chrono::system_clock::now());

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
    SupportReactionsFile.open("Outputs/AnalysisOutputs/SupportReactionsFile");
    for (int i = 0; i < supportReactions.size(); ++i)
    {
        double supportReaction = supportReactions.at(i);
        SupportReactionsFile << supportReaction;
        SupportReactionsFile << "\n";
    }
    SupportReactionsFile.close();
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
        DisplacementFile <<"Displacement at Dof ";
        DisplacementFile << i + 1;
        DisplacementFile <<" = ";
        DisplacementFile << dispVector.at(i);
        DisplacementFile << " m";
        DisplacementFile << "\n";

    }

    StiffnessMatrixFile.close();
    ForceVectorFile.close();
    DisplacementFile.close();
    
    cout<<"Displacements are calculated"<<endl;

    PrincipleStressCalculator pSC(elmVec, dispVector, e, v, meshSize);
    std::vector<std::vector<double>> principleStressVector(elmVec.size());
    principleStressVector = pSC.PrincipleStressList;
    
    cout<< "Principle stresses are calculated"<<endl;
	vector<vector<double>> MainStresses = MainStressCalculator(elmVec, principleStressVector);

	ofstream MainStressesFile; 
	MainStressesFile.open("Outputs/AnalysisOutputs/MainStressesFile");
	for (int i = 0; i < MainStresses.size(); ++i)
	{
		vector<double> stressCouple = MainStresses.at(i);
		//MainStressesFile<<"Element Index ";
		//MainStressesFile<<(i + 1);
		//MainStressesFile<<"\n";
		//MainStressesFile<<"Compressive Stress = ";
		MainStressesFile<<stressCouple.at(0) * 0.000001;
		//MainStressesFile<<" MPa";
		MainStressesFile<<" ";		
		//MainStressesFile<<"\n";
		//MainStressesFile<<"Tensile Stress = ";
		MainStressesFile<<stressCouple.at(1) * 0.000001;
		//MainStressesFile<<" MPa";
		MainStressesFile<<"\n";
	}
	MainStressesFile.close();

	// Element modification part

    vector<Element> ModifiedElmVec = ElementModificator(elmVec, dispVector, e, v, meshSize, thickness, fGlobal, nDof, nDofRestrained, nodeVec, controlDof, controlDisplacement, nbcVector);
    StiffnessMatrixAssembler modifiedSMA(ModifiedElmVec, nDof);
    vector<vector<double>> modifiedKGlobal = modifiedSMA.GetGlobalStiffnessMatrix(ModifiedElmVec, nDof);
    ForceVectorAssembler modifiedFVA(nodeVec, nbcVector, meshSize, nDof);
    vector<double> modifiedFGlobal = modifiedFVA.ForceVector;
    DisplacementCalculator modifiedDispCal(modifiedKGlobal, modifiedFGlobal, nDof, nDofRestrained);
    std::vector<double> modifiedDispVector = modifiedDispCal.DisplacementVector;
    PrincipleStressCalculator modifiedPSC(ModifiedElmVec, modifiedDispVector, e, v, meshSize);
    vector<vector<double>> modifiedPrincipleStressVector = modifiedPSC.PrincipleStressList;
    vector<vector<double>> modifiedCompressiveTensileStresses = MainStressCalculator(ModifiedElmVec, modifiedPrincipleStressVector);

    double monitoredDisp = modifiedDispVector.at(controlDof - 1);
    double convRatio = monitoredDisp / controlDisplacement;
    int counter = 0;
    
    vector<Element> controlElmVector = ModifiedElmVec;
	for (int i = 0; i < 10; ++i)
	{
		cout<<"Elements are being modified";    	
		cout<<"\n";
		ModifiedElmVec =  ElementModificator(controlElmVector, dispVector, e, v, meshSize, thickness, fGlobal, nDof, nDofRestrained, nodeVec, controlDof, controlDisplacement, nbcVector);
    	StiffnessMatrixAssembler loopdSMA(ModifiedElmVec, nDof);
    	modifiedKGlobal = loopdSMA.GetGlobalStiffnessMatrix(ModifiedElmVec, nDof);
    	ForceVectorAssembler loopFVA(nodeVec, nbcVector, meshSize, nDof);
    	modifiedFGlobal = loopFVA.ForceVector;
    	DisplacementCalculator loopDispCalc(modifiedKGlobal, modifiedFGlobal, nDof, nDofRestrained);
    	modifiedDispVector = loopDispCalc.DisplacementVector;
    	PrincipleStressCalculator loopPSC(ModifiedElmVec, modifiedDispVector, e, v, meshSize);
    	modifiedPrincipleStressVector = loopPSC.PrincipleStressList;
    	modifiedCompressiveTensileStresses = MainStressCalculator(ModifiedElmVec, modifiedPrincipleStressVector);
    	monitoredDisp = modifiedDispVector.at(controlDof - 1);
		convRatio = monitoredDisp / controlDisplacement;
    	controlElmVector = ModifiedElmVec;

    	if (IsEqual(convRatio, 1, 0.1))
    	{
    		break;
    	}
    }

    cout<<"Element modification according to principal stresses"<<endl;
    cout<<"and displacement correction is completed"<<endl;

    ofstream ModifiedMainStressFile;
    ModifiedMainStressFile.open("Outputs/AnalysisOutputs/ModifiedMainStressFile");

    double nom = 0;
    double denom = 0;

    double maxCompressiveStress;

    for (int i = 0; i < modifiedCompressiveTensileStresses.size(); ++i)
    {
    	vector<double> stressCouple = modifiedCompressiveTensileStresses.at(i);
    	double sigmaMin = stressCouple.at(0);
    	double sigmaMax = stressCouple.at(1);    	
    	ModifiedMainStressFile << sigmaMin * 0.000001;
    	ModifiedMainStressFile << " ";
    	ModifiedMainStressFile << sigmaMax * 0.000001;
    	ModifiedMainStressFile << "\n";

        if (sigmaMax > tol)
        {
            nom += sigmaMax;
            denom++;
        }

        if (sigmaMin < maxCompressiveStress)
        {
            maxCompressiveStress = sigmaMin;
        }
    }
    cout<<"Analysis is completed"<<endl;
    cout<<"See Outputs file for Analysis Outputs"<<endl;
    cout<<"Stable system is obtained. Design procedure starts..."<<endl;

    ofstream ShearDesignOutput;
    ofstream TensileDesignOutput;
    TensileDesignOutput.open("Outputs/DesignOutputs/TensileDesignOutputs");    
    ShearDesignOutput.open("Outputs/DesignOutputs/ShearDesignOutputs");

    // Tensile Design
    if (IsEqual(denom, 0.0, tol))
    {
        cout<<"There is no element in tension. Only shear reinforcements are designed"<<endl;
        TensileDesignOutput<<"There is no elements in tension. Only shear reinforcements are designed"<<endl;
    }
    else
    {
        TensileDesignOutput<<"Note that Yield Stress is multiplied by 0.75 in design";
        TensileDesignOutput<<"\n";
        TensileDesignOutput<<"-------------------------------------------------------";
        TensileDesignOutput<<"\n";
        TensileDesignOutput<<"\n";

        double diaMM = longitudinalBarDiameter * 1000;
        double barArea = 3.141593 * longitudinalBarDiameter * longitudinalBarDiameter / 4;
        double averageTensileStress = nom / denom;
        for (int i = 0; i < ModifiedElmVec.size(); ++i)
        {
            Element elm = ModifiedElmVec.at(i);
            vector<double> stressCouple = modifiedCompressiveTensileStresses.at(i);
            double tensileStress = stressCouple.at(1);
            
            if (tensileStress > averageTensileStress)
            {
                TensileDesignOutput<<"Tensile Design For Element#";
                TensileDesignOutput<<elm.ElementIndex;
                TensileDesignOutput<<"\n";
                double requiredArea = (tensileStress / (fY * 0.75)) * (thickness * meshSize);
                TensileDesignOutput<<"Required Area = ";
                TensileDesignOutput<<requiredArea * 1000000;
                TensileDesignOutput<<" mm^2";
                TensileDesignOutput<<"\n";
                TensileDesignOutput<<"Use ";
                TensileDesignOutput<< ceil(requiredArea / barArea);
                TensileDesignOutput<<"Ø";
                TensileDesignOutput<<diaMM;
                TensileDesignOutput<<"\n";
                TensileDesignOutput<<"Provided Area = ";
                TensileDesignOutput<< ceil(requiredArea / barArea) * barArea * 1000000;
                TensileDesignOutput<< " mm^2";
                TensileDesignOutput<<"\n";
                TensileDesignOutput<<"A_provided > A_required -> OK ✔";
                TensileDesignOutput<<"\n";
                TensileDesignOutput<<"\n";

            }
        }
    }

    // Shear Design
    ShearDesignOutput<<"Stress check...";
    ShearDesignOutput<<"\n";
    ShearDesignOutput<<"ØFn = Ø * 0.85 * fc where Ø = 075";
    ShearDesignOutput<<"\n";
    ShearDesignOutput<<"Maximum Compressive Stress for an element = ";
    ShearDesignOutput<<-1 * maxCompressiveStress * 0.000001;
    ShearDesignOutput<<" MPa";
    ShearDesignOutput<<"\n";
    ShearDesignOutput<<"ØFn = ";
    ShearDesignOutput<<0.75 * 0.85 * fC * 0.000001;
    ShearDesignOutput<<" MPa";
    ShearDesignOutput<<"\n";
    ShearDesignOutput<<"Fu = ";
    ShearDesignOutput<<-1 * maxCompressiveStress * 0.000001;
    ShearDesignOutput<<" MPa";
    ShearDesignOutput<<"\n";
    
    if (-1 * maxCompressiveStress > (0.75 * 0.85 * fC))
    {
        ShearDesignOutput<<"ØFn < Fu";
        ShearDesignOutput<<"\n";
        ShearDesignOutput<<"Compressive stress is above limits. You may increase thickness or use concrete of higher strength...";
    }
    else
    {
        double transverseDiaMM = transverseBarDiameter * 1000;
        ShearDesignOutput<<"ØFn => Fu";
        ShearDesignOutput<<"\n";
        ShearDesignOutput<<"Compressive stress is within limits. Design can be conducted.";
        ShearDesignOutput<<"\n";
        ShearDesignOutput<<"\n";
        ShearDesignOutput<<"Vertical Reinforcement Design";
        ShearDesignOutput<<"\n";
        int nRow = ceil(thickness / 0.1) - 1;
        double rho200 = CalculateShearReinforcementRatio(transverseBarDiameter, nRow, 0.2, thickness);
        double rho100 = CalculateShearReinforcementRatio(transverseBarDiameter, nRow, 0.1, thickness);
        double rho50 = CalculateShearReinforcementRatio(transverseBarDiameter, nRow, 0.05, thickness);

        if (rho200 > 0.0025)
        {
            ShearDesignOutput<<"Use ";
            ShearDesignOutput<<nRow;
            ShearDesignOutput<<" Ø";
            ShearDesignOutput<<transverseDiaMM;
            ShearDesignOutput<<"/200 mm";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"(p_v = ";
            ShearDesignOutput<<rho200;
            ShearDesignOutput<<" > 0.0025  -> OK ✔";
        }
        else if (rho100 > 0.0025)
        {
            ShearDesignOutput<<"Use ";
            ShearDesignOutput<<nRow;
            ShearDesignOutput<<" Ø";
            ShearDesignOutput<<transverseDiaMM;
            ShearDesignOutput<<"/100 mm";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"(p_v = ";
            ShearDesignOutput<<rho100;
            ShearDesignOutput<<" > 0.0025  -> OK ✔";
        }
        else if (rho50 > 0.0025)
        {
            ShearDesignOutput<<"Use ";
            ShearDesignOutput<<nRow;
            ShearDesignOutput<<" Ø";
            ShearDesignOutput<<transverseDiaMM;
            ShearDesignOutput<<"/50 mm";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"(p_v = ";
            ShearDesignOutput<<rho50;
            ShearDesignOutput<<" > 0.0025  -> OK ✔";
        }
        else
        {
            ShearDesignOutput<<"Could not find a configuration for vertical shear reinforcement.";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"Please increase transverse reinforcement bar diameter";
        }
        
        ShearDesignOutput<<"\n";
        ShearDesignOutput<<"\n";
        ShearDesignOutput<<"Horizontal Reinforcement Design";
        ShearDesignOutput<<"\n";
        double rho200Hor = CalculateShearReinforcementRatio(transverseBarDiameter, nRow, 0.2, thickness);
        double rho100Hor = CalculateShearReinforcementRatio(transverseBarDiameter, nRow, 0.1, thickness);
        double rho50Hor= CalculateShearReinforcementRatio(transverseBarDiameter, nRow, 0.05, thickness);

        if (rho200Hor > 0.0025)
        {
            ShearDesignOutput<<"Use ";
            ShearDesignOutput<<nRow;
            ShearDesignOutput<<" Ø";
            ShearDesignOutput<<transverseDiaMM;
            ShearDesignOutput<<"/200 mm";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"(p_v = ";
            ShearDesignOutput<<rho200Hor;
            ShearDesignOutput<<" > 0.0025  -> OK ✔";
        }
        else if (rho100Hor > 0.0025)
        {
            ShearDesignOutput<<"Use ";
            ShearDesignOutput<<nRow;
            ShearDesignOutput<<" Ø";
            ShearDesignOutput<<transverseDiaMM;
            ShearDesignOutput<<"/100 mm";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"(p_v = ";
            ShearDesignOutput<<rho100Hor;
            ShearDesignOutput<<" > 0.0025  -> OK ✔";
        }
        else if (rho50Hor > 0.0025)
        {
            ShearDesignOutput<<"Use ";
            ShearDesignOutput<<nRow;
            ShearDesignOutput<<" Ø";
            ShearDesignOutput<<transverseDiaMM;
            ShearDesignOutput<<"/50 mm";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"(p_v = ";
            ShearDesignOutput<<rho50Hor;
            ShearDesignOutput<<" > 0.0025  -> OK ✔";
        }
        else
        {
            ShearDesignOutput<<"Could not find a configuration for horizontal shear reinforcement.";
            ShearDesignOutput<<"\n";
            ShearDesignOutput<<"Please increase transverse reinforcement bar diameter";  
        }
    }    

    cout<<"Design is completed"<<endl;
    cout<<"See Outputs file for Design Outputs"<<endl;
    cout<<"\n";
    TensileDesignOutput.close();
    ShearDesignOutput.close();
    auto timenow2 =
            chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout<< "Elapsed Time = " << timenow2 - timenow << " seconds"<< endl;
    return 0;
}
