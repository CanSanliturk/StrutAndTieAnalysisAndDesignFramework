#include <iostream>
#include <vector>
#include <chrono>
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
    double tol = 0.000001; // Set tolerance value to check equality

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
    cout<<"Meshes are created"<<endl;

    int nDof = 0;
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
    }

    int nDofRestrained = 0;
    int nDofUnrestrained = 0;
    for (int i = 0; i < ebcVector.size(); ++i)
    {
        EssentialBC essentialBc = ebcVector.at(i);
        bool isPoint = IsEqual(essentialBc.XStart, essentialBc.XEnd, tol) && IsEqual(essentialBc.YStart, essentialBc.YEnd, tol);
        bool isXDir = (!IsEqual(essentialBc.XStart, essentialBc.XEnd, tol)) && IsEqual(essentialBc.YStart, essentialBc.YEnd, tol);
        bool isYDir = IsEqual(essentialBc.XStart, essentialBc.XEnd, tol) && (!IsEqual(essentialBc.YStart, essentialBc.YEnd, tol));

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
    StiffnessMatrixAssembler sMA(elmVec, nDof);
    vector<vector<double>> kGlobal = sMA.GetGlobalStiffnessMatrix(elmVec, nDof);
    ForceVectorAssembler fVA(nodeVec, nbcVector, meshSize, nDof);
    std::vector<double> fGlobal = fVA.ForceVector;

    cout<<"Assembly of matrices and vectors are completed"<<endl;

    DisplacementCalculator dispCal(kGlobal, fGlobal, nDof, nDofRestrained);
    std::vector<double> dispVector = dispCal.DisplacementVector;
    
    cout<<"Displacements are calculated"<<endl;

    PrincipleStressCalculator pSC(elmVec, dispVector, e, v, meshSize);
    std::vector<std::vector<double>> principleStressVector(elmVec.size());
    principleStressVector = pSC.PrincipleStressList;
    int stressListSize = principleStressVector.size();

    
    for (int i = 0; i < stressListSize; ++i)
    {
    	std::vector<double> elmStress = principleStressVector.at(i);
    	double sigmaXX = elmStress.at(0) * 0.000001;
    	double sigmaYY = elmStress.at(1) * 0.000001;
    	double sigmaXY = elmStress.at(2) * 0.000001;

    	cout<<"--------------------------------------------------------------"<<endl;
    	cout<<"Element No: "<< (i + 1) <<endl;
    	cout<<"Stress in XX-Direction = "<<sigmaXX<<" MPa"<<endl;
    	cout<<"Stress in YY-Direction = "<<sigmaYY<<" MPa"<<endl;
    	cout<<"Stress in XY-Direction = "<<sigmaXY<<" MPa"<<endl;
    }
	

    auto timenow2 =
            chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << ctime(&timenow2) << endl;
    cout<< "Elapsed Time = " << timenow2 - timenow << " seconds"<< endl;

    
    
    return 0;
}
