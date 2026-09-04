#include "Mesh1D/Mesh1D.hpp"

#include <cstddef>
#include <iostream>
#include <stdexcept>



Mesh1D::Mesh1D(
    const Geometry& HeightObj,
    std::size_t NumOfNodes_
)
    : NumOfNodes(
          (NumOfNodes_ >= 2)
              ? NumOfNodes_
              : throw std::invalid_argument("The number of mesh nodes must be at least 2.")
      ),
      NumOfIntervals(NumOfNodes - 1),
      ChannelHeight(HeightObj.GetValue()),
      deltaY(ChannelHeight / static_cast<double>(NumOfIntervals)),
      YPositions(NumOfNodes)
{
    DiscretizeYDirection();
}

void Mesh1D::DiscretizeYDirection()
{
    for (std::size_t Index = 0; Index < YPositions.size(); ++Index)
    {
        YPositions[Index] = static_cast<double>(Index) * deltaY;
    }
}

void Mesh1D::PrintYPositions() const
{
    for (std::size_t Index = 0; Index < YPositions.size(); ++Index)
    {
        std::cout << YPositions[Index] << " ";
    }

    std::cout << '\n';
}

std::size_t Mesh1D::GetNumOfNodes() const
{
    return NumOfNodes;
}

void Mesh1D::SetNumOfNodes(std::size_t NumOfNodes_)
{
    if (NumOfNodes_ < 2)
    {
        throw std::invalid_argument("The number of mesh nodes must be at least 2.");
    }

    NumOfNodes = NumOfNodes_;
    NumOfIntervals = NumOfNodes - 1;
    deltaY = ChannelHeight / static_cast<double>(NumOfIntervals);
    YPositions.resize(NumOfNodes);
    DiscretizeYDirection();
}


std::size_t Mesh1D::GetNumOfIntervals() const
{
    return NumOfIntervals;
}

double Mesh1D::GetdeltaY() const
{
    return deltaY;
}

double Mesh1D::GetYPositions(std::size_t Index) const
{
    return YPositions.at(Index);
}



std::size_t Mesh1D::find_Position_Index(double target_Position) const
{
    const double tolerance = 1e-9;

    auto it = std::find_if(YPositions.begin(), YPositions.end(), [target_Position, tolerance](double pos) {
        return std::abs(pos - target_Position) < tolerance;
    });

    if (it == YPositions.end()) 
    {
        return static_cast<std::size_t>(-1);
    }

    return static_cast<std::size_t>(it - YPositions.begin());
}
