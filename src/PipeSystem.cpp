#include "PipeSystem.h"
#include "MatrixUtils.h"
#include <fstream>
#include <iostream>
#include "json.hpp"


using json = nlohmann::json;

static double K1[12][12]{}, M1[12][12]{}, C[12][12]{}, c[3][3]{}, CT[12][12]{};
static double MultiResult[12][12]{};

// 解析单个 DOF 约束，添加 static 限定
static DofConstraint parse_dof(const json& j) {
    return {
        j.value("node", 0),
        j.value("connecting_node", 0),
        j.value("direction_x", 0.0),
        j.value("direction_y", 0.0),
        j.value("direction_z", 0.0),
        j.value("friction", 0.0),
        j.value("gap", 0.0),
        j.value("stiffness", 0.0),
        j.value("support_guid", std::string{}),
        j.value("support_tag", std::string{}),
        j.value("type", 0)
    };
}

// 解析完整的 restraint 数组，添加 static 限定
static std::vector<Restraint> parse_restraints(const json& root) {
    std::vector<Restraint> restraints;
    if (!root.contains("boundary_condition") || !root["boundary_condition"].contains("restraints"))
        return restraints;

    const auto& j_restraints = root["boundary_condition"]["restraints"];
    for (const auto& r : j_restraints) {
        Restraint restraint;
        for (const auto& dof : r["dofs"]) {
            restraint.dofs.push_back(parse_dof(dof));
        }
        restraints.push_back(restraint);
    }
    return restraints;
    
}


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

    restraints = parse_restraints(j);

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
        for (size_t i = 0; i < restraints.size(); ++i) {
        const auto& r = restraints[i];
        for (size_t j = 0; j < r.dofs.size(); ++j) {
            const auto& dof = r.dofs[j];
            if (dof.type == 0) continue;

            // 定义DOF类型到索引的映射
            int dof_index = 0;
            switch (dof.type) {
            case 1: dof_index = 6; break;
            case 2: dof_index = 5; break;
            case 3: dof_index = 4; break;
            case 10: dof_index = 3; break;
            case 11: dof_index = 2; break;
            case 12: dof_index = 1; break;
            default: throw std::runtime_error("nonlinear restraint");
            }

            int row_idx = number[dof.node] * 6 - dof_index;

            // 添加对角线项
            StiffnessMatrix(row_idx, row_idx) += dof.stiffness;

            // 如果有连接节点，添加耦合项
            if (dof.connecting_node != 0) {
                int col_idx = number[dof.connecting_node] * 6 - dof_index;

                StiffnessMatrix(col_idx, col_idx) += dof.stiffness;

                // 添加反对角线项（负值）
                StiffnessMatrix(row_idx, col_idx) -= dof.stiffness;
                StiffnessMatrix(col_idx, row_idx) -= dof.stiffness;

            }
        }
    }
}
