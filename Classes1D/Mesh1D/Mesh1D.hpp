#ifndef Mesh1D_HPP
#define Mesh1D_HPP

#include <cstddef>
#include <vector>
#include <algorithm> // std::find
#include <cmath>
#include "Geometry/Geometry.hpp"

class Mesh1D
{


public:
    explicit Mesh1D(
        const Geometry& HeightObj,
        std::size_t NumOfNodes_ = 10
    );

    void PrintYPositions() const;

    std::size_t GetNumOfNodes() const;
    void SetNumOfNodes(std::size_t NumOfNodes_);

    std::size_t GetNumOfIntervals() const;
    double GetdeltaY() const;
    double GetYPositions(std::size_t Index) const;
    std::size_t find_Position_Index(double target_Position) const;


    private:
    std::size_t NumOfNodes;
    std::size_t NumOfIntervals;
    double ChannelHeight;
    double deltaY;
    std::vector<double> YPositions;

    void DiscretizeYDirection();
};

#endif // Mesh1D_HPP
