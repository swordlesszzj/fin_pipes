#include "PipeSystem.h"
#include "MatrixUtils.h"
#include <fstream>
#include <iostream>
#include "include/json.hpp"


using json = nlohmann::json;

static double K1[12][12]{}, M1[12][12]{}, C[12][12]{}, c[3][3]{}, CT[12][12]{};
static double MultiResult[12][12]{};

PipeSystem::PipeSystem() : StiffnessMatrix(0, 0), MassMatrix(0, 0) { }

void PipeSystem::resetTempMatrices() {
    memset(K1, 0, sizeof(K1));
    memset(M1, 0, sizeof(M1));
    memset(C, 0, sizeof(C));
    memset(c, 0, sizeof(c));
    memset(CT, 0, sizeof(CT));
    memset(MultiResult, 0, sizeof(MultiResult));
}

void PipeSystem::loadFromFile(const string& filename) {
    pipes.clear();
    id.clear();
    number.clear();

    ifstream file(filename);
    if (!file) {
        throw runtime_error("无法打开文件: " + filename);
    }

    json j;
    file >> j;

    const auto& pipe_system = j["pipe_system"];
    const auto& pipe_jsons = pipe_system["pipes"];

    int Count_pipes2 = 0;
    for (const auto& pipe_json : pipe_jsons) {
        const auto& elements = pipe_json["elements"];
        for (const auto& elem : elements) {
            Count_pipes2++;
            Pipe p;
            p.StartNode = elem.value("from_node", -1);
            p.EndNode = elem.value("to_node", -1);
            p.DeltaX = elem.value("delta_x", 0.0) / 1000;
            p.DeltaY = elem.value("delta_y", 0.0) / 1000;
            p.DeltaZ = elem.value("delta_z", 0.0) / 1000;
            p.Diameter = elem.value("diameter", 0.0) / 1000;
            p.WallThickness = elem.value("wall_thickness", 0.0) / 1000;
            p.ElasticModulus = elem.value("elastic_modulus_cold", 0.0);
            if (p.ElasticModulus == 0) p.ElasticModulus = 200000000000;
            p.PoissionRatio = elem.value("poisson_ratio", 0.0);
            if (p.PoissionRatio == 0) p.PoissionRatio = 0.3;

            p.Density = 0.0;
            if (elem.contains("fluid_density") && elem["fluid_density"].is_array() && !elem["fluid_density"].empty()) {
                p.Density = elem["fluid_density"][0].get<double>();
            }

            p.computeProperties();
            pipes.push_back(p);

            // 记录节点编号顺序
            if (number.find(p.StartNode) == number.end()) {
                number[p.StartNode] = id.size() + 1;
                id.push_back(p.StartNode);
            }
            if (number.find(p.EndNode) == number.end()) {
                number[p.EndNode] = id.size() + 1;
                id.push_back(p.EndNode);
            }
        }
    }

    // 初始化全局矩阵大小
    StiffnessMatrix = Matrix(id.size() * 6, id.size() * 6);
    MassMatrix = Matrix(id.size() * 6, id.size() * 6);
}

void PipeSystem::assembleMatrices() {
    resetTempMatrices();

    for (const Pipe& p1 : pipes) {
        // 计算局部刚度矩阵 K1 和质量矩阵 M1
        // 集中质量矩阵
        for (int i = 0; i <= 8; i++) {
            if (i <= 2) M1[i][i] = p1.Mass / 2;
            if (i >= 6) M1[i][i] = p1.Mass / 2;
        }
        // 刚度矩阵
        K1[0][0] = p1.ElasticModulus * p1.Aera / p1.Length; K1[6][0] = -p1.ElasticModulus * p1.Aera / p1.Length;
        K1[1][1] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1[5][1] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1[7][1] = -12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1[11][1] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1[2][2] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1[4][2] = -6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1[8][2] = -12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1[10][2] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1[3][3] = p1.J_x * p1.ShearModulus / p1.Length; K1[9][3] = -p1.J_x * p1.ShearModulus / p1.Length;
        K1[4][4] = 4 * p1.ElasticModulus * p1.J_y / p1.Length; K1[8][4] = 6 * p1.ElasticModulus * p1.J_y / p1.Length; K1[10][4] = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1[5][5] = 4 * p1.ElasticModulus * p1.J_y / p1.Length; K1[7][5] = -6 * p1.ElasticModulus * p1.J_y / p1.Length; K1[11][5] = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1[6][6] = p1.ElasticModulus * p1.Aera / p1.Length;
        K1[7][7] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3); K1[11][7] = -6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1[8][8] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3); K1[10][8] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1[9][9] = p1.J_x * p1.ShearModulus / p1.Length;
        K1[10][10] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1[11][11] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        symmetry(K1, 12);

        // 计算旋转矩阵 c
        float h = std::sqrt(p1.DeltaX * p1.DeltaX + p1.DeltaZ * p1.DeltaZ);
        if (h == 0) {
            c[0][0] = p1.DeltaX / p1.Length; c[0][1] = p1.DeltaY / p1.Length; c[0][2] = p1.DeltaZ / p1.Length;
            c[1][0] = 0; c[1][1] = h / p1.Length; c[1][2] = 1;
            c[2][0] = 1; c[2][1] = 0; c[2][2] = 0;
        }
        else {
            c[0][0] = p1.DeltaX / p1.Length; c[0][1] = p1.DeltaY / p1.Length; c[0][2] = p1.DeltaZ / p1.Length;
            c[1][0] = -p1.DeltaX * p1.DeltaY / p1.Length / h; c[1][1] = h / p1.Length; c[1][2] = -p1.DeltaZ * p1.DeltaY / p1.Length / h;
            c[2][0] = -p1.DeltaZ / h; c[2][1] = 0; c[2][2] = p1.DeltaX / h;
        }

        // 构造大尺寸旋转矩阵 C 和转置 CT
        for (int i = 0; i <= 3; i++) {
            for (int j = 0; j <= 2; j++) {
                for (int k = 0; k <= 2; k++) {
                    C[4 * i + j][4 * i + k] = c[j][k];
                }
            }
        }
        for (int i = 0; i <= 11; i++) {
            for (int j = 0; j <= 11; j++) {
                CT[i][j] = C[j][i];
            }
        }

        int id_start = number[p1.StartNode] - 1;
        int id_end = number[p1.EndNode] - 1;

        // 计算旋转后的矩阵
        Multi(CT, M1, C, MultiResult);
        for (int i = 0; i <= 5; i++) {
            for (int j = 0; j <= 5; j++) {
                MassMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult[i][j];
            }
        }
        for (int i = 6; i <= 11; i++) {
            for (int j = 0; j <= 5; j++) {
                MassMatrix(i + 6 * id_end - 6, j + 6 * id_start) += MultiResult[i][j];
            }
        }
        for (int i = 0; i <= 5; i++) {
            for (int j = 6; j <= 11; j++) {
                MassMatrix(i + 6 * id_start, j + 6 * id_end - 6) += MultiResult[i][j];
            }
        }
        for (int i = 6; i <= 11; i++) {
            for (int j = 6; j <= 11; j++) {
                MassMatrix(i + 6 * id_end - 6, j + 6 * id_end - 6) += MultiResult[i][j];
            }
        }

        Multi(CT, K1, C, MultiResult);
        for (int i = 0; i <= 5; i++) {
            for (int j = 0; j <= 5; j++) {
                StiffnessMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult[i][j];
            }
        }
        for (int i = 6; i <= 11; i++) {
            for (int j = 0; j <= 5; j++) {
                StiffnessMatrix(i + 6 * id_end - 6, j + 6 * id_start) += MultiResult[i][j];
            }
        }
        for (int i = 0; i <= 5; i++) {
            for (int j = 6; j <= 11; j++) {
                StiffnessMatrix(i + 6 * id_start, j + 6 * id_end - 6) += MultiResult[i][j];
            }
        }
        for (int i = 6; i <= 11; i++) {
            for (int j = 6; j <= 11; j++) {
                StiffnessMatrix(i + 6 * id_end - 6, j + 6 * id_end - 6) += MultiResult[i][j];
            }
        }
    }
}
