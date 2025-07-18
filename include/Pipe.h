#pragma once
#include <cmath>

const float Pi = 3.1415926;

struct Pipe {
    int Id;
    int StartNode;
    int EndNode;
    float DeltaX;
    float DeltaY;
    float DeltaZ;
    float Diameter;         // 外径
    float WallThickness;
    float ElasticModulus;
    float PoissionRatio;
    float Density;

    // 计算派生属性
    float I_x;     // 转动惯量
    float J_x, J_y; // 惯性矩
    float ShearModulus;
    float Mass;
    float Length;
    float Aera;

    void computeProperties() {
        ShearModulus = ElasticModulus / 2 / (1 + PoissionRatio);
        Length = std::sqrt(DeltaX*DeltaX + DeltaY*DeltaY + DeltaZ*DeltaZ);
        Aera = Pi / 4 * (std::pow(Diameter, 2) - std::pow((Diameter - 2 * WallThickness), 2));
        Mass = Aera * Density * Length;
        I_x = Pi / 32 * (std::pow(Diameter, 4) - std::pow((Diameter - 2 * WallThickness), 4)) * Density * Length;
        J_x = Pi / 32 * (std::pow(Diameter, 4) - std::pow((Diameter - 2 * WallThickness), 4));
        J_y = J_x / 2;
    }
};
