/*===========
C++ Libraries
=============*/ 
#include <iostream>
#include <vector>
#include <cmath> 
#include <iomanip>
#include <thread> // Added for optional small delays

/*=====
Headers
========*/ 

//Stops The Loop Solver Iterations by Pressing ESC
#include "EssentialHeaders/KeyboardHandler.hpp"

//Inputs
#include "EssentialHeaders/SolverInputs.hpp"


/*=====
Classes
========*/

//Geometry
#include "Classes1D/Geometry/Geometry.hpp"

//Thermophysical_Properties
#include "Classes1D/Thermophysical_Properties/ThermophysicalProperties.hpp"

//SolverSettings
#include "Classes1D/SolverSettings/SolutionParameters.hpp"

//Mesh1D
#include "Classes1D/Mesh1D/Mesh1D.hpp"

//Field1D
#include "Classes1D/Field1D/Field1D.hpp"

//InitialConditions
#include "Classes1D/InitialConditions/InitialConditions.hpp"

//Boundaries
#include "Classes1D/Boundaries/Boundary.hpp"

//BoundaryConditions
#include "Classes1D/BoundaryConditions/BoundaryConditions.hpp"

//AnalyticalSolution
#include "Classes1D/AnalyticalSolution1D/AnalyticalDiffusion1D.hpp"

//FileWriter1D
#include "Classes1D/FileWriter1D/FileWriter1D.hpp"

//Plotter
#include "Classes1D/Plotter/ResidualPlotter.hpp"
#include "Classes1D/Plotter/OutputPlotter1D.hpp"

//Matrices
#include "Classes1D/Matrices/CoefficientMatrix.hpp"
#include "Classes1D/Matrices/TridiagonalMatrix.hpp"
#include "Classes1D/Matrices/RHS.hpp"

//Linear Solvers
#include "Classes1D/LinearSolvers/GaussianElimination.hpp"
#include "Classes1D/LinearSolvers/ThomasAlgorithm.hpp"

//Convergence
#include "Classes1D/Convergence/RelativeResidual.hpp"

//Printer
#include "Classes1D/Printer/StatusPrinter.hpp"


enum class SolverScheme 
{
    Laasonen,
    CrankNicolson
};

//constexpr SolverScheme ActiveScheme = SolverScheme::Laasonen;
constexpr SolverScheme ActiveScheme = SolverScheme::CrankNicolson;
//============
//Start Of main
//=============
int main()
{

const std::string schemeName = (ActiveScheme == SolverScheme::Laasonen) 
                                  ? "Laasonen" 
                                  : "CrankNicolson";

std::cout << "Selected Scheme: " << schemeName << std::endl;

//=======
//Geometry
//========= 
Geometry Height;
Height.SetValue(SolverInputs::Height_Value);


//======================
//ThermophysicalProperties
//========================
ThermophysicalProperties nu;
nu.SetValue(SolverInputs::nu_Value);


//===============
//Solver Settings
//===============

  SolutionParameters<std::size_t> NumTimeLevels;
  SolutionParameters<double> dt;
  SolutionParameters<double> Tolerance;
 


  NumTimeLevels.SetValue(SolverInputs::NumTimeStep_Value);
  dt.SetValue(SolverInputs::dt_Value);
  Tolerance.SetValue(SolverInputs::Tolerance_Value); //Tolerance is used for Convergence
  

//======
//Mesh1D
//====== 
Mesh1D mesh1D(Height);
mesh1D.SetNumOfNodes(SolverInputs::NumOfNodes_Value);
//Printer1D::PrintMesh1D(mesh1D);

//===============
//Diffusion Number
//================= 
SolutionParameters<double> DiffNumber;
DiffNumber.SetValue((nu.GetValue() * dt.GetValue())/(mesh1D.GetdeltaY()*mesh1D.GetdeltaY()));
  
    
//=====
//Field1D
//======= 
Field1D U_0(mesh1D);
Field1D U_n(mesh1D);
Field1D U_nPlus1(mesh1D);

//Printer1D::PrintField1D(U_n);

//=====================================
//Apply Initial Conditions To Field1D U_0
//======================================
InitialConditions::ApplyInitialCondition(U_0, SolverInputs::Uniform_Initial_Velocity);

//====================
// Boundaries For Field1D
// ======================
Boundary TopWall;
Boundary BottomWall;

//Boundary Setters
TopWall.SetValue(SolverInputs::TopWall_Velocity_Values);
BottomWall.SetValue(SolverInputs::BottomWall_Velocity_Values);

//=======================================
//Apply Boundary Conditions To Field1D U_0
//=========================================

BoundaryConditions::ApplyBoundaryCondition(U_0,TopWall,BoundaryLocation::Top);
BoundaryConditions::ApplyBoundaryCondition(U_0,BottomWall,BoundaryLocation::Bottom);


//Printer1D::Print_Mesh1D_Field1D(mesh1D,U_0);

// =================
// AnalyticalSolution
// ===================
AnalyticalDiffusion1D U_Analytical_1D(mesh1D,Height,TopWall,BottomWall);

// ======================================================
// Writing Initial Velocity and Analytical Solution To file
// ========================================================

std::string FilePath = "OutPutData/" + schemeName;
FileWriter1D::WriteField1D(U_0, mesh1D, FieldType::Initial_Numerical, FilePath);
FileWriter1D::WriteField1D(U_Analytical_1D, mesh1D, FieldType::Steady_State_Analytical, FilePath);


// =====================================
// Initialize Plotter and Data Containers
// ======================================
ResidualPlotter plotter("CFD Convergence: Diffusion 1D (" + schemeName + ")");
std::vector<double> stepHistory;
std::vector<double> residualHistory;

// =============================================================================
// Create Coefficient  Matrix A : 2D Matrix
// When you create this 2D matrix, you should use the Gaussian elimination method.
// ===============================================================================
CoefficientMatrix A(Scheme::Laasonen,mesh1D,DiffNumber);
//CoefficientMatrix A(Scheme::CrankNicolson,mesh1D,DiffNumber);


//===========================================================
//Apply Dirichlet boundary condition To Coefficient  Matrix A 
//===========================================================
A.SetBoundaryConditions();

A.Print();



//=========================================================
// Create Tridiagonal  Matrix .
// When you create this class, you should use Thomas' algorithm.
// Tridiagonal Matrix Holds three 1D Vectors:
//  1.LowerDiagonal;
//  2.MainDiagonal;
//  3.UpperDiagonal;
// =============================================================
//TridiagonalMatrix A3(SchemeName::Laasonen,mesh1D,DiffNumber);
TridiagonalMatrix A3(SchemeName::CrankNicolson,mesh1D,DiffNumber);



//=======================================================================
//Apply Dirichlet boundary condition To TridiagonalMatrix A3 
//=======================================================================
A3.SetBoundaryConditions();


// ==========
// Create RHS
// ==========

if (ActiveScheme == SolverScheme::Laasonen)
{
// ==============
// Lassonen Method
// ===============
// EQ : -r * U[i - 1][n + 1] + (1.0 + 2.0 * r) * U[i][n + 1] - r * U[i + 1][n + 1] = RHS(U_n)
// RHS = U[i][n]
// Create RHS for Laasonen Method
// ========================================
RHS RHS1D(U_0);
}

else if (ActiveScheme == SolverScheme::CrankNicolson)
{
// =================
// CrankNicolson Method
// ===========================================================================================================
//  EQ : -(r / 2) * U[i - 1]^(n + 1) + (1 + r) * U[i]^(n + 1) - (r / 2) * U[i + 1]^(n + 1) = RHS(U_n,r)
//  RHS(U_n,r) = (r / 2) * U[i - 1]^n + (1 - r) * U[i]^n + (r / 2) * U[i + 1]^n
//  Create RHS for CrankNicolson Method  
// ===========================================================================================================
RHS RHS1D(U_0,DiffNumber);
}


//=======
//StartUp
//=======
U_n.Swap(U_0);


// ==========================
// Solver Loop 
// ===========================

for (std::size_t TimeLevel = 1; TimeLevel < NumTimeLevels.GetValue(); TimeLevel++)
{
//This Condition Stops The Loop Solver Iterations by Pressing ESC
if (isEscPressed())
{
    std::cout << "\nESC pressed. Exiting program now..." << std::endl;
    return 0;
}

// TODO: Temporarily copying matrix A to prevent in-place modification by Gaussian Elimination.
// While this copy operation adds O(N^2) overhead per time step, it keeps the current CFD solver 
// implementation straightforward. In future optimization phases, this should be replaced with 
// LU Decomposition (A = L*U). By decomposing A only once outside the time loop, we can reduce 
// the per-step solver cost from O(N^3) (Gaussian Elimination) to O(N^2) (Forward/Backward Substitution), 
// eliminating both the N^3 solver overhead and the need for copying matrix A at each step.


// =====================
// Gaussian_Elimination
// ====================
// CoefficientMatrix A_copy = A; // LU Decomposition is better insted of this statement
// RHS RHS1D(U_n);
// GaussianElimination::Solve_nPlus1(A_copy,RHS1D,U_nPlus1);


// ================
// Thomas Algorithm
// =================
TridiagonalMatrix A3_Copy = A3;
RHS RHS1D(U_n);
ThomasAlgorithm::Solve_nPlus1(A3_Copy, RHS1D, U_nPlus1);


//===================================
// Apply Boundary Conditions To Field1D
// =====================================
// BoundaryConditions::ApplyBoundaryCondition(U_nPlus1,mesh1D,TopWall);
// BoundaryConditions::ApplyBoundaryCondition(U_nPlus1,mesh1D,BottomWall);


//Claculate The Relative Residual
double Residual = RelativeResidual::ReturnResidual(U_n,U_nPlus1);


// Print Step Status such as Time // Step // Relative Residual
StatusPrinter::PrintStepStatus(TimeLevel,Residual,dt);

//U <--- Unew : U=Unew Copy The Values But Swap does not Copy
U_n.Swap(U_nPlus1);



// Record data for plotting
stepHistory.push_back(static_cast<double>(TimeLevel));
residualHistory.push_back(Residual);

// Update plot every 20 steps
if (TimeLevel % 20 == 0) 
{
    plotter.updatePlot(stepHistory, residualHistory);
}


//Codition For Convergence To stady state Solution
if(Residual < Tolerance.GetValue())
{
StatusPrinter::PrintConvergenceStatus(TimeLevel,Residual,dt);
FileWriter1D::WriteField1D(U_n,mesh1D,FieldType::Steady_State_Numerical,FilePath);
break;
}


}//end of Solver Loop

// =========================================================================
// Common Post-Processing
// =========================================================================
OutputPlotter1D outputPlotter("Solution Comparison (" + schemeName + ")");
outputPlotter.plotComparison(
    FilePath + "/Initial_Numerical1D.dat", "Initial Condition",
    FilePath + "/Steady_State_Numerical1D.dat", "Steady State Numerical Solution",
    FilePath + "/Steady_State_Analytical1D.dat", "Steady State Analytical Solution"
);
  std::cout << "Calculations Completed Successfully" <<std::endl;
  std::cin.get();
  return 0;
  
}//End Of main
