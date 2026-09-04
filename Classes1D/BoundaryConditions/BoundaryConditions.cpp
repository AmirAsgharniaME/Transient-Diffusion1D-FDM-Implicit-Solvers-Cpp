#include "BoundaryConditions/BoundaryConditions.hpp" // class BoundaryConditions;  forward declaration

void BoundaryConditions::ApplyBoundaryCondition(Field1D& Field1D_Obj,const Mesh1D& Mesh1D_Obj,const Boundary& Boundary_Obj)
{
    const double Boundary_Value = Boundary_Obj.GetValue();
    double  BoundaryPosition = Boundary_Obj.GetPosition();
    
    Field1D_Obj.SetValue(Mesh1D_Obj.find_Position_Index(BoundaryPosition),Boundary_Value);
    
}
