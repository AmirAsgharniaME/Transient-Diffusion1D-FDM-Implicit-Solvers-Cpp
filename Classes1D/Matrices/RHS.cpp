#include "Matrices/RHS.hpp"
//#include <algorithm>

RHS::RHS(
        const Field1D& Field1D_Obj,
        double Default_Initial_Value_)
        :RHSValues(Field1D_Obj.Size(), Default_Initial_Value_)
    
{
    // ==============
    // Lassonen Method
    // ===============
    // EQ : -r * U[i - 1][n + 1] + (1.0 + 2.0 * r) * U[i][n + 1] - r * U[i + 1][n + 1] = RHS(U_n)
    // RHS = U[i][n]
    // Create RHS for Laasonen Method
    // ========================================
    const std::size_t NumOfNodes = Field1D_Obj.Size();

     // Iterate through the interior nodes
    for (std::size_t i = 1; i < NumOfNodes-1; i++)
    {
            RHSValues[i] = Field1D_Obj.GetValue(i);
    }
    //Apply Dirichlet boundary condition To RHS
     RHSValues[0] = Field1D_Obj.GetValue(0);
     RHSValues[NumOfNodes-1] = Field1D_Obj.GetValue(NumOfNodes-1);
    }
    
RHS::RHS(
        const Field1D& Field1D_Obj,
        const SolutionParameters<double>& DiffNumer_Obj,
        double Default_Initial_Value_ 
    )
    :RHSValues(Field1D_Obj.Size(), Default_Initial_Value_),
     DiffNumer(DiffNumer_Obj.GetValue())
    {

    // =================
    // CrankNicolson Method
    // ===========================================================================================================
    //  EQ : -(r / 2) * U[i - 1]^(n + 1) + (1 + r) * U[i]^(n + 1) - (r / 2) * U[i + 1]^(n + 1) = RHS(U_n,r)
    //  RHS(U_n,r) = (r / 2) * U[i - 1]^n + (1 - r) * U[i]^n + (r / 2) * U[i + 1]^n
    //  Create RHS for CrankNicolson Method  
    // ===========================================================================================================
  
    // Create a local copy of RHS to act as the source stencil (n-th time level).
    // This is crucial to prevent 'Stencil Pollution', where updated values from the 
    // current time-step (n+1) would otherwise corrupt the stencil for subsequent spatial nodes.
    //Field1D U_old_local = RHS;

    const std::size_t NumOfNodes = Field1D_Obj.Size();

    // Iterate through the interior nodes, excluding boundary nodes.
    for (size_t i = 1; i < NumOfNodes-1; i++)
    {

        // Reading strictly from U_old_local preserves the n-th time level data.
        double rhs =(DiffNumer / 2.0) *   Field1D_Obj.GetValue(i - 1)
                    + (1.0 - DiffNumer) * Field1D_Obj.GetValue(i)
                    + (DiffNumer / 2.0) * Field1D_Obj.GetValue(i + 1);

        // Update the RHS vector with the new calculated implicit source term.
        RHSValues[i]=rhs;
    }

    //Apply Dirichlet boundary condition To RHS
    RHSValues[0]=Field1D_Obj.GetValue(0);
    RHSValues[Field1D_Obj.Size()-1]=Field1D_Obj.GetValue(Field1D_Obj.Size()-1);

    }


void RHS::SetValue(std::size_t PIndex, double Value_)
{
    RHSValues[PIndex] = Value_;
}

double RHS::GetValue(std::size_t PIndex) const
{

    return RHSValues[PIndex];
}

std::size_t RHS::GetSize() const
{
return RHSValues.size();
}