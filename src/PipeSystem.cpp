// #include "PipeSystem.h"
// #include "MatrixUtils.h"
// #include <fstream>
// #include <iostream>
// #include "json.hpp"

// void orthogonalize(double c[3][3]);

// using json = nlohmann::json;

// static double K1[12][12]{}, M1[12][12]{}, C[12][12]{}, c[3][3]{}, CT[12][12]{};
// static double MultiResult[12][12]{};

// // 解析单个 DOF 约束，添加 static 限定
// static DofConstraint parse_dof(const json &j)
// {
//     return {
//         j.value("node", 0),
//         j.value("connecting_node", 0),
//         j.value("direction_x", 0.0),
//         j.value("direction_y", 0.0),
//         j.value("direction_z", 0.0),
//         j.value("friction", 0.0),
//         j.value("gap", 0.0),
//         j.value("stiffness", 0.0),
//         j.value("support_guid", std::string{}),
//         j.value("support_tag", std::string{}),
//         j.value("type", 0)};
// }

// // 解析完整的 restraint 数组，添加 static 限定
// static std::vector<Restraint> parse_restraints(const json &root)
// {
//     std::vector<Restraint> restraints;
//     if (!root.contains("boundary_condition") || !root["boundary_condition"].contains("restraints"))
//         return restraints;

//     const auto &j_restraints = root["boundary_condition"]["restraints"];
//     for (const auto &r : j_restraints)
//     {
//         Restraint restraint;
//         for (const auto &dof : r["dofs"])
//         {
//             restraint.dofs.push_back(parse_dof(dof));
//         }
//         restraints.push_back(restraint);
//     }
//     return restraints;
// }

// PipeSystem::PipeSystem() : StiffnessMatrix(0, 0), MassMatrix(0, 0) {}

// void PipeSystem::resetTempMatrices()
// {
//     memset(K1, 0, sizeof(K1));
//     memset(M1, 0, sizeof(M1));
//     memset(C, 0, sizeof(C));
//     memset(c, 0, sizeof(c));
//     memset(CT, 0, sizeof(CT));
//     memset(MultiResult, 0, sizeof(MultiResult));
// }

// void PipeSystem::loadFromFile(const string &filename)
// {
//     pipes.clear();
//     id.clear();
//     number.clear();

//     ifstream file(filename);
//     if (!file)
//     {
//         throw runtime_error("无法打开文件: " + filename);
//     }

//     json j;
//     file >> j;

//     const auto &pipe_system = j["pipe_system"];
//     const auto &pipe_jsons = pipe_system["pipes"];

//     int Count_pipes2 = 0;
//     for (const auto &pipe_json : pipe_jsons)
//     {
//         const auto &elements = pipe_json["elements"];
//         for (const auto &elem : elements)
//         {
//             Count_pipes2++;
//             Pipe p;
//             p.StartNode = elem.value("from_node", -1);
//             p.EndNode = elem.value("to_node", -1);
//             p.DeltaX = elem.value("delta_x", 0.0) / 1000;
//             p.DeltaY = elem.value("delta_y", 0.0) / 1000;
//             p.DeltaZ = elem.value("delta_z", 0.0) / 1000;
//             p.Diameter = elem.value("diameter", 0.0) / 1000;
//             p.WallThickness = elem.value("wall_thickness", 0.0) / 1000;
//             p.ElasticModulus = elem.value("elastic_modulus_cold", 0.0);
//             if (p.ElasticModulus == 0)
//                 p.ElasticModulus = 200000000000;
//             p.PoissionRatio = elem.value("poisson_ratio", 0.0);
//             if (p.PoissionRatio == 0)
//                 p.PoissionRatio = 0.3;

//             p.Density = 0.0;
//             if (elem.contains("fluid_density") && elem["fluid_density"].is_array() && !elem["fluid_density"].empty())
//             {
//                 p.Density = elem["fluid_density"][0].get<double>();
//             }

//             p.computeProperties();
//             pipes.push_back(p);

//             // 记录节点编号顺序
//             if (number.find(p.StartNode) == number.end())
//             {
//                 number[p.StartNode] = id.size() + 1;
//                 id.push_back(p.StartNode);
//             }
//             if (number.find(p.EndNode) == number.end())
//             {
//                 number[p.EndNode] = id.size() + 1;
//                 id.push_back(p.EndNode);
//             }
//         }
//     }

//     restraints = parse_restraints(j);

//     // 初始化全局矩阵大小
//     StiffnessMatrix = Matrix(id.size() * 6, id.size() * 6);
//     MassMatrix = Matrix(id.size() * 6, id.size() * 6);
// }

// void PipeSystem::assembleMatrices()
// {
//     resetTempMatrices();

//     for (const Pipe &p1 : pipes)
//     {
//         // 计算局部刚度矩阵 K1 和质量矩阵 M1

//                 //一致质量矩阵
//         M1[0][0] = p1.Mass / 3; M1[6][0] = p1.Mass / 6;
//         M1[1][1] = p1.Mass * 13 / 35; M1[5][1] = p1.Mass*p1.Length * 11 / 210; M1[7][1] = p1.Mass * 9 / 70; M1[11][1] = -p1.Mass * p1.Length * 13 / 420;
//         M1[2][2] = p1.Mass * 13 / 35; M1[4][2] = -p1.Mass * p1.Length * 11 / 210; M1[8][2] = p1.Mass * 9 / 70; M1[10][2] = p1.Mass * p1.Length * 13 / 420;
//         M1[3][3] = p1.I_x / 3; M1[9][3] = p1.I_x / 6;
//         M1[4][4] = p1.Mass / 105 * pow(p1.Length, 2); M1[8][4] = -13 * p1.Mass / 420 * p1.Length; M1[10][4] = -p1.Mass / 140 * pow(p1.Length, 2);
//         M1[5][5] = p1.Mass / 105 * pow(p1.Length, 2); M1[7][5] = 13 * p1.Mass / 420 * p1.Length; M1[11][5] = -p1.Mass / 140 * pow(p1.Length, 2);
//         M1[6][6] = p1.Mass / 3;
//         M1[7][7] = p1.Mass * 13 / 35; M1[11][7] = -p1.Mass * p1.Length * 11 / 210;
//         M1[8][8] = p1.Mass * 13 / 35; M1[10][8] = p1.Mass * p1.Length * 11 / 210;
//         M1[9][9] = p1.I_x / 3;
//         M1[10][10] = p1.Mass / 105 * pow(p1.Length, 2);
//         M1[11][11] = p1.Mass / 105 * pow(p1.Length, 2);
//         symmetry(M1,12);
        

//         //集中质量矩阵
//         // for (int i = 0; i <= 8; i++)
//         // {
//         //     if (i <= 2)
//         //         M1[i][i] = p1.Mass / 2;
//         //     if (i >= 6)
//         //         M1[i][i] = p1.Mass / 2;
//         // }
//         // 刚度矩阵
//         K1[0][0] = p1.ElasticModulus * p1.Aera / p1.Length;
//         K1[6][0] = -p1.ElasticModulus * p1.Aera / p1.Length;
//         K1[1][1] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
//         K1[5][1] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
//         K1[7][1] = -12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
//         K1[11][1] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
//         K1[2][2] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
//         K1[4][2] = -6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
//         K1[8][2] = -12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
//         K1[10][2] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
//         K1[3][3] = p1.J_x * p1.ShearModulus / p1.Length;
//         K1[9][3] = -p1.J_x * p1.ShearModulus / p1.Length;
//         K1[4][4] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
//         K1[8][4] = 6 * p1.ElasticModulus * p1.J_y / p1.Length / p1.Length;
//         K1[10][4] = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
//         K1[5][5] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
//         K1[7][5] = -6 * p1.ElasticModulus * p1.J_y / p1.Length / p1.Length;
//         K1[11][5] = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
//         K1[6][6] = p1.ElasticModulus * p1.Aera / p1.Length;
//         K1[7][7] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
//         K1[11][7] = -6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
//         K1[8][8] = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
//         K1[10][8] = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
//         K1[9][9] = p1.J_x * p1.ShearModulus / p1.Length;
//         K1[10][10] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
//         K1[11][11] = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
//         symmetry(K1, 12);

//                 // 计算旋转矩阵 c
//         float h = std::sqrt(p1.DeltaX * p1.DeltaX + p1.DeltaZ * p1.DeltaZ);
//         if (fabs(h/p1.Length) < 1e-8)
//         {
//             c[0][0] = p1.DeltaX / p1.Length;
//             c[0][1] = p1.DeltaY / p1.Length;
//             c[0][2] = p1.DeltaZ / p1.Length;
//             c[1][0] = 0;
//             c[1][1] = h / p1.Length;
//             c[1][2] = 1;
//             c[2][0] = 1;
//             c[2][1] = 0;
//             c[2][2] = 0;
//         }
//         else
//         {
//             c[0][0] = p1.DeltaX / p1.Length;
//             c[0][1] = p1.DeltaY / p1.Length;
//             c[0][2] = p1.DeltaZ / p1.Length;
//             c[1][0] = -p1.DeltaX * p1.DeltaY / p1.Length / h;
//             c[1][1] = h / p1.Length;
//             c[1][2] = -p1.DeltaZ * p1.DeltaY / p1.Length / h;
//             c[2][0] = -p1.DeltaZ / h;
//             c[2][1] = 0;
//             c[2][2] = p1.DeltaX / h;
//         }
//         orthogonalize(c);
//         double cc[3][3];
//         memset(cc,0,sizeof(cc));
//         for(int i=0;i<=2;i++){
//             for(int j=0;j<=2;j++){
//                 for(int k=0;k<=2;k++){
//                     cc[i][j]=c[i][k]*c[j][k]+cc[i][j];
//                 }
//                 if(i==j){
//                     if(fabs(cc[i][j]-1)>=1e-8) 
//                     cout<<"不正交";
//                 }
//                 else if(fabs(cc[i][j])>=1e-8)
//                 cout<<"不正交";
//             }
//         }
    

//         // 构造大尺寸旋转矩阵 C 和转置 CT
//         for (int i = 0; i <= 3; i++)
//         {
//             for (int j = 0; j <= 2; j++)
//             {
//                 for (int k = 0; k <= 2; k++)
//                 {
//                     C[3 * i + j][3 * i + k] = c[j][k];
//                 }
//             }
//         }
//         for (int i = 0; i <= 11; i++)
//         {
//             for (int j = 0; j <= 11; j++)
//             {
//                 CT[i][j] = C[j][i];
//             }
//         }

//         int id_start = number[p1.StartNode] - 1;
//         int id_end = number[p1.EndNode] - 1;
//                double cc1[12][12];
//         memset(cc1,0,sizeof(cc1));
//         for(int i=0;i<=11;i++){
//             for(int j=0;j<=11;j++){
//                 for(int k=0;k<=11;k++){
//                     cc1[i][j]=C[i][k]*C[j][k]+cc1[i][j];
//                 }
//                 if(i==j){
//                     if(fabs(cc1[i][j]-1)>=1e-5){
//                         cout<<fabs(cc1[i][j]-1)<<endl;
//                         cout<<"不正交";
//                     }
//                 }
//                 else if(fabs(cc1[i][j])>=1e-5){
//                     cout<<fabs(cc1[i][j])<<endl;
//                     cout<<"不正交";
//                 }
//             }
//         }

//         // 计算旋转后的矩阵
//         Multi(CT, M1, C, MultiResult);
//         for (int i = 0; i <= 5; i++)
//         {
//             for (int j = 0; j <= 5; j++)
//             {
//                 MassMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult[i][j];
//             }
//         }
//         for (int i = 6; i <= 11; i++)
//         {
//             for (int j = 0; j <= 5; j++)
//             {
//                 MassMatrix(i + 6 * id_end - 6, j + 6 * id_start) += MultiResult[i][j];
//             }
//         }
//         for (int i = 0; i <= 5; i++)
//         {
//             for (int j = 6; j <= 11; j++)
//             {
//                 MassMatrix(i + 6 * id_start, j + 6 * id_end - 6) += MultiResult[i][j];
//             }
//         }
//         for (int i = 6; i <= 11; i++)
//         {
//             for (int j = 6; j <= 11; j++)
//             {
//                 MassMatrix(i + 6 * id_end - 6, j + 6 * id_end - 6) += MultiResult[i][j];
//             }
//         }

//         Multi(CT, K1, C, MultiResult);
//         Matrix A(12,12);
//         for(int i=0;i<=11;i++){
//             for(int j=0;j<=11;j++){
//                 A(i,j)=K1[i][j];
//             }
//         }
//         Matrix L=A.choleskyDecomposition(A);
        


//         for (int i = 0; i <= 5; i++)
//         {
//             for (int j = 0; j <= 5; j++)
//             {
//                 StiffnessMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult[i][j];
//             }
//         }
//         for (int i = 6; i <= 11; i++)
//         {
//             for (int j = 0; j <= 5; j++)
//             {
//                 StiffnessMatrix(i + 6 * id_end - 6, j + 6 * id_start) += MultiResult[i][j];
//             }
//         }
//         for (int i = 0; i <= 5; i++)
//         {
//             for (int j = 6; j <= 11; j++)
//             {
//                 StiffnessMatrix(i + 6 * id_start, j + 6 * id_end - 6) += MultiResult[i][j];
//             }
//         }
//         for (int i = 6; i <= 11; i++)
//         {
//             for (int j = 6; j <= 11; j++)
//             {
//                 StiffnessMatrix(i + 6 * id_end - 6, j + 6 * id_end - 6) += MultiResult[i][j];
//             }
//         }
//     }
//     for (size_t i = 0; i < restraints.size(); ++i)
//     {
//         const auto &r = restraints[i];
//         for (size_t j = 0; j < r.dofs.size(); ++j)
//         {
//             const auto &dof = r.dofs[j];
//             if (dof.type == 0)
//                 continue;

//             // 定义DOF类型到索引的映射
//             int dof_index = 0;
//             switch (dof.type)
//             {
//             case 2:
//                 dof_index = 6;
//                 break;
//             case 3:
//                 dof_index = 5;
//                 break;
//             case 4:
//                 dof_index = 4;
//                 break;
//             case 11:
//                 dof_index = 3;
//                 break;
//             case 12:
//                 dof_index = 2;
//                 break;
//             case 13:
//                 dof_index = 1;
//                 break;
//             case 14:
//                 dof_index = 5;
//                 break;
//             case 1:
//                 dof_index = 6;
//                 break;
//             default:
//                 std::printf("nonlinear restraint");
//                 continue;
//             }

//             int row_idx = number[dof.node] * 6 - dof_index;
//             if (dof.type != 1)
//             {
//                 // 添加对角线项
//                 StiffnessMatrix(row_idx, row_idx) += dof.stiffness;

//                 // 如果有连接节点，添加耦合项
//                 if (dof.connecting_node != 0)
//                 {
//                     int col_idx = number[dof.connecting_node] * 6 - dof_index;

//                     StiffnessMatrix(col_idx, col_idx) += dof.stiffness;

//                     // 添加反对角线项（负值）
//                     StiffnessMatrix(row_idx, col_idx) -= dof.stiffness;
//                     StiffnessMatrix(col_idx, row_idx) -= dof.stiffness;
//                 }
//             }
//             else
//             {
//                 for (int i = 0; i <= 5; i++)
//                 {
//                     StiffnessMatrix(row_idx + i, row_idx + i) += dof.stiffness;
//                     if (dof.connecting_node != 0)
//                     {
//                         int col_idx = number[dof.connecting_node] * 6 + i;
//                         StiffnessMatrix(col_idx, col_idx) += dof.stiffness;
//                         StiffnessMatrix(row_idx + i, col_idx) -= dof.stiffness;
//                         StiffnessMatrix(col_idx, row_idx + i) -= dof.stiffness;
//                     }
//                 }
//             }
//         }
//     }
// }

// // 对3x3矩阵c进行正交化
// void orthogonalize(double c[3][3]) {
//     // Gram-Schmidt正交化
//     // 取三列向量
//     double v0[3] = {c[0][0], c[1][0], c[2][0]};
//     double v1[3] = {c[0][1], c[1][1], c[2][1]};
//     double v2[3] = {c[0][2], c[1][2], c[2][2]};

//     // 正交化v1
//     double dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
//     for(int i=0;i<3;++i) v1[i] -= dot01 * v0[i];

//     // 归一化v0,v1
//     double norm0 = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
//     double norm1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
//     for(int i=0;i<3;++i) {
//         v0[i] /= norm0;
//         v1[i] /= norm1;
//     }

//     // v2正交化
//     double dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2];
//     double dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
//     for(int i=0;i<3;++i) v2[i] -= dot02 * v0[i] + dot12 * v1[i];
//     double norm2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
//     for(int i=0;i<3;++i) v2[i] /= norm2;

//     // 回写
//     for(int i=0;i<3;++i) {
//         c[i][0] = v0[i];
//         c[i][1] = v1[i];
//         c[i][2] = v2[i];
//     }
// }
#include "PipeSystem.h"
#include "MatrixUtils.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "json.hpp"
#include "Matrix.h"

void orthogonalize(Matrix &c);

using json = nlohmann::json;

// 全部用Matrix类替换
static Matrix K1(12, 12), M1(12, 12), C(12, 12), c(3, 3), CT(12, 12), MultiResult(12, 12);

// 解析单个 DOF 约束，添加 static 限定
static DofConstraint parse_dof(const json &j)
{
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
        j.value("type", 0)};
}

// 解析完整的 restraint 数组，添加 static 限定
static std::vector<Restraint> parse_restraints(const json &root)
{
    std::vector<Restraint> restraints;
    if (!root.contains("boundary_condition") || !root["boundary_condition"].contains("restraints"))
        return restraints;

    const auto &j_restraints = root["boundary_condition"]["restraints"];
    for (const auto &r : j_restraints)
    {
        Restraint restraint;
        for (const auto &dof : r["dofs"])
        {
            restraint.dofs.push_back(parse_dof(dof));
        }
        restraints.push_back(restraint);
    }
    return restraints;
}

PipeSystem::PipeSystem() : StiffnessMatrix(0, 0), MassMatrix(0, 0) {}

void PipeSystem::resetTempMatrices()
{
    K1 = Matrix(12, 12);
    M1 = Matrix(12, 12);
    C = Matrix(12, 12);
    CT = Matrix(12, 12);
    MultiResult = Matrix(12, 12);
    c = Matrix(3, 3);
}

void PipeSystem::loadFromFile(const string &filename)
{
    pipes.clear();
    id.clear();
    number.clear();

    ifstream file(filename);
    if (!file)
    {
        throw runtime_error("无法打开文件: " + filename);
    }

    json j;
    file >> j;

    const auto &pipe_system = j["pipe_system"];
    const auto &pipe_jsons = pipe_system["pipes"];

    int Count_pipes2 = 0;
    for (const auto &pipe_json : pipe_jsons)
    {
        const auto &elements = pipe_json["elements"];
        for (const auto &elem : elements)
        {
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
            if (p.ElasticModulus == 0)
                p.ElasticModulus = 200000000000;
            p.PoissionRatio = elem.value("poisson_ratio", 0.0);
            if (p.PoissionRatio == 0)
                p.PoissionRatio = 0.3;

            
            if (elem.contains("fluid_density") && elem["fluid_density"].is_array() && !elem["fluid_density"].empty())
            {
                p.Density = elem["fluid_density"][0].get<double>();
            }
            if (p.Density == 0)
                p.Density = 7850;


            p.computeProperties();
            pipes.push_back(p);

            // 记录节点编号顺序
            if (number.find(p.StartNode) == number.end())
            {
                number[p.StartNode] = id.size() + 1;
                id.push_back(p.StartNode);
            }
            if (number.find(p.EndNode) == number.end())
            {
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

void PipeSystem::assembleMatrices()
{
    resetTempMatrices();

    for (const Pipe &p1 : pipes)
    {
        // 计算局部刚度矩阵 K1 和质量矩阵 M1
       // 一致质量矩阵
        for (int i = 0; i <= 8; i++)
        {
            if (i <= 2)
                M1(i,i) = p1.Mass / 2;
            if (i >= 6)
                M1(i,i) = p1.Mass / 2;
        }
        // M1(0,0) = p1.Mass / 3; M1(6,0) = p1.Mass / 6;
        // M1(1,1) = p1.Mass * 13 / 35; M1(5,1) = p1.Mass*p1.Length * 11 / 210; M1(7,1) = p1.Mass * 9 / 70; M1(11,1) = -p1.Mass * p1.Length * 13 / 420;
        // M1(2,2) = p1.Mass * 13 / 35; M1(4,2) = -p1.Mass * p1.Length * 11 / 210; M1(8,2) = p1.Mass * 9 / 70; M1(10,2) = p1.Mass * p1.Length * 13 / 420;
        // M1(3,3) = p1.I_x / 3; M1(9,3) = p1.I_x / 6;
        // M1(4,4) = p1.Mass / 105 * pow(p1.Length, 2); M1(8,4) = -13 * p1.Mass / 420 * p1.Length; M1(10,4) = -p1.Mass / 140 * pow(p1.Length, 2);
        // M1(5,5) = p1.Mass / 105 * pow(p1.Length, 2); M1(7,5) = 13 * p1.Mass / 420 * p1.Length; M1(11,5) = -p1.Mass / 140 * pow(p1.Length, 2);
        // M1(6,6) = p1.Mass / 3;
        // M1(7,7) = p1.Mass * 13 / 35; M1(11,7) = -p1.Mass * p1.Length * 11 / 210;
        // M1(8,8) = p1.Mass * 13 / 35; M1(10,8) = p1.Mass * p1.Length * 11 / 210;
        // M1(9,9) = p1.I_x / 3;
        // M1(10,10) = p1.Mass / 105 * pow(p1.Length, 2);
        // M1(11,11) = p1.Mass / 105 * pow(p1.Length, 2);
        // for (int i = 0; i <= 10; i++) {
        //      for (int j = i + 1; j <= 11; j++) {
        //        M1(i, j) = M1(j, i);
        //     }
        // }
        
        // 刚度矩阵
        K1(0,0) = p1.ElasticModulus * p1.Aera / p1.Length;
        K1(6,0) = -p1.ElasticModulus * p1.Aera / p1.Length;
        K1(1,1) = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1(5,1) = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1(7,1) = -12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1(11,1) = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1(2,2) = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1(4,2) = -6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1(8,2) = -12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1(10,2) = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1(3,3) = p1.J_x * p1.ShearModulus / p1.Length;
        K1(9,3) = -p1.J_x * p1.ShearModulus / p1.Length;
        K1(4,4) = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1(8,4) = 6 * p1.ElasticModulus * p1.J_y / p1.Length / p1.Length;
        K1(10,4) = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1(5,5) = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1(7,5) = -6 * p1.ElasticModulus * p1.J_y / p1.Length / p1.Length;
        K1(11,5) = 2 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1(6,6) = p1.ElasticModulus * p1.Aera / p1.Length;
        K1(7,7) = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1(11,7) = -6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1(8,8) = 12 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 3);
        K1(10,8) = 6 * p1.ElasticModulus * p1.J_y / std::pow(p1.Length, 2);
        K1(9,9) = p1.J_x * p1.ShearModulus / p1.Length;
        K1(10,10) = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        K1(11,11) = 4 * p1.ElasticModulus * p1.J_y / p1.Length;
        for (int i = 0; i <= 10; i++) {
             for (int j = i + 1; j <= 11; j++) {
               K1(i, j) = K1(j, i);
            }
        }
        //cout<<K1;
        // 计算旋转矩阵 c
        float h = std::sqrt(p1.DeltaX * p1.DeltaX + p1.DeltaZ * p1.DeltaZ);
        if (fabs(h/p1.Length) < 1e-8)
        {
            c(0,0) = 0;
            c(0,1) = p1.DeltaY/p1.Length;
            c(0,2) = 0;
            c(1,0) = 0;
            c(1,1) = 0;
            c(1,2) = 1;
            c(2,0) = p1.DeltaY/p1.Length;
            c(2,1) = 0;
            c(2,2) = 0;
        }
        else
        {
            c(0,0) = p1.DeltaX / p1.Length;
            c(0,1) = p1.DeltaY / p1.Length;
            c(0,2) = p1.DeltaZ / p1.Length;
            c(1,0) = -p1.DeltaX * p1.DeltaY / p1.Length / h;
            c(1,1) = h / p1.Length;
            c(1,2) = -p1.DeltaZ * p1.DeltaY / p1.Length / h;
            c(2,0) = -p1.DeltaZ / h;
            c(2,1) = 0;
            c(2,2) = p1.DeltaX / h;
        }
        orthogonalize(c);

        // 构造大尺寸旋转矩阵 C 和转置 CT
        for (int i = 0; i <= 3; i++)
        {
            for (int j = 0; j <= 2; j++)
            {
                for (int k = 0; k <= 2; k++)
                {
                    C(3 * i + j, 3 * i + k) = c(j, k);
                }
            }
        }
        for (int i = 0; i <= 11; i++)
        {
            for (int j = 0; j <= 11; j++)
            {
                CT(i, j) = C(j, i);
            }
        }

        int id_start = number[p1.StartNode] - 1;
        int id_end = number[p1.EndNode] - 1;

        // 计算旋转后的矩阵
        MultiResult = CT * M1 * C;
        for (int i = 0; i <= 5; i++)
        {
            for (int j = 0; j <= 5; j++)
            {
                MassMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult(i, j);
            }
        }
        for (int i = 6; i <= 11; i++)
        {
            for (int j = 0; j <= 5; j++)
            {
                MassMatrix(i + 6 * id_end - 6, j + 6 * id_start) += MultiResult(i, j);
            }
        }
        for (int i = 0; i <= 5; i++)
        {
            for (int j = 6; j <= 11; j++)
            {
                MassMatrix(i + 6 * id_start, j + 6 * id_end - 6) += MultiResult(i, j);
            }
        }
        for (int i = 6; i <= 11; i++)
        {
            for (int j = 6; j <= 11; j++)
            {
                MassMatrix(i + 6 * id_end - 6, j + 6 * id_end - 6) += MultiResult(i, j);
            }
        }

        MultiResult = CT * K1 * C;

        for (int i = 0; i <= 5; i++)
        {
            for (int j = 0; j <= 5; j++)
            {
                StiffnessMatrix(i + 6 * id_start, j + 6 * id_start) += MultiResult(i, j);
            }
        }
        for (int i = 6; i <= 11; i++)
        {
            for (int j = 0; j <= 5; j++)
            {
                StiffnessMatrix(i + 6 * id_end - 6, j + 6 * id_start) += MultiResult(i, j);
            }
        }
        for (int i = 0; i <= 5; i++)
        {
            for (int j = 6; j <= 11; j++)
            {
                StiffnessMatrix(i + 6 * id_start, j + 6 * id_end - 6) += MultiResult(i, j);
            }
        }
        for (int i = 6; i <= 11; i++)
        {
            for (int j = 6; j <= 11; j++)
            {
                StiffnessMatrix(i + 6 * id_end - 6, j + 6 * id_end - 6) += MultiResult(i, j);
            }
        }
    }
    for (size_t i = 0; i < restraints.size(); ++i)
    {
        const auto &r = restraints[i];
        for (size_t j = 0; j < r.dofs.size(); ++j)
        {
            const auto &dof = r.dofs[j];
            if (dof.type == 0)
                continue;

            // 定义DOF类型到索引的映射
            int dof_index = 0;
            switch (dof.type)
            {
            case 2:
                dof_index = 6;
                break;
            case 3:
                dof_index = 5;
                break;
            case 4:
                dof_index = 4;
                break;
            case 11:
                dof_index = 3;
                break;
            case 12:
                dof_index = 2;
                break;
            case 13:
                dof_index = 1;
                break;
            case 14:
                dof_index = 5;
                break;
            case 1:
                dof_index = 6;
                break;
            default:
                std::printf("nonlinear restraint");
                continue;
            }

            int row_idx = number[dof.node] * 6 - dof_index;
            if (dof.type != 1)
            {
                // 添加对角线项
                StiffnessMatrix(row_idx, row_idx) += dof.stiffness;

                // 如果有连接节点，添加耦合项
                if (dof.connecting_node != 0)
                {
                    int col_idx = number[dof.connecting_node] * 6 - dof_index;

                    StiffnessMatrix(col_idx, col_idx) += dof.stiffness;

                    // 添加反对角线项（负值）
                    StiffnessMatrix(row_idx, col_idx) -= dof.stiffness;
                    StiffnessMatrix(col_idx, row_idx) -= dof.stiffness;
                }
            }
            else
            {
                for (int i = 0; i <= 5; i++)
                {
                    StiffnessMatrix(row_idx + i, row_idx + i) += dof.stiffness;
                    if (dof.connecting_node != 0)
                    {
                        int col_idx = number[dof.connecting_node] * 6 + i;
                        StiffnessMatrix(col_idx, col_idx) += dof.stiffness;
                        StiffnessMatrix(row_idx + i, col_idx) -= dof.stiffness;
                        StiffnessMatrix(col_idx, row_idx + i) -= dof.stiffness;
                    }
                }
            }
        }
    }
int n = id.size() * 6;
int Num = 0;

//第一步：收集需要保留和移除的索引
std::vector<int> keep_indices;
std::vector<int> remove_indices;
keep_indices.reserve(n);  // 预分配内存以提高性能
remove_indices.reserve(n);

for (int i = 0; i < n; ++i) {
    bool should_keep = false;
    
    // 检查整行中是否有任何元素大于阈值
    for (int j = 0; j < n; ++j) {
        if (std::fabs(MassMatrix(i, j)) > 1e-28) {
            should_keep = true;
            break;  // 找到一个大于阈值的元素就足够
        }
    }
    
    if (should_keep) {
        keep_indices.push_back(i);
        Num++;
    } else {
        remove_indices.push_back(i);
    }
}

int num_keep = keep_indices.size();
int num_remove = remove_indices.size();

// 第二步：创建索引映射表
std::vector<int> index_map(n);
for (int i = 0; i < num_keep; ++i) {
    index_map[keep_indices[i]] = i;
}
for (int i = 0; i < num_remove; ++i) {
    index_map[remove_indices[i]] = num_keep + i;
}

// 第三步：直接重新排列矩阵元素
Matrix new_MassMatrix(num_keep, num_keep);  // 只需要保留部分
Matrix new_StiffnessMatrix(n, n);

// 并行化循环以提高大矩阵处理速度
//#pragma omp parallel for collapse(2)
for (int i = 0; i < n; ++i) {
    int new_i = index_map[i];
    for (int j = 0; j < n; ++j) {
        int new_j = index_map[j];
        new_StiffnessMatrix(new_i, new_j) = StiffnessMatrix(i, j);
        
        // 只处理MassMatrix的保留部分
        if (new_i < num_keep && new_j < num_keep) {
            new_MassMatrix(new_i, new_j) = MassMatrix(i, j);
        }
        else{if(fabs(MassMatrix(i,j)>=1e-15)) cout<<MassMatrix(i,j);}
    }
}

// 第四步：交换矩阵内容
MassMatrix = std::move(new_MassMatrix);
StiffnessMatrix = std::move(new_StiffnessMatrix);

// 提取子矩阵
int remove_size = n - Num;
Matrix Kcc(remove_size, remove_size);
Matrix Kca(remove_size, Num);
Matrix Kaa(Num, Num);

for (int i = 0; i < remove_size; ++i) {
    int global_i = i + Num;
    for (int j = 0; j < remove_size; ++j) {
        int global_j = j + Num;
        Kcc(i, j) = StiffnessMatrix(global_i, global_j);
    }
    
    for (int j = 0; j < Num; ++j) {
        Kca(i, j) = StiffnessMatrix(global_i, j);
    }
}

for (int i = 0; i < Num; ++i) {
    for (int j = 0; j < Num; ++j) {
        Kaa(i, j) = StiffnessMatrix(i, j);
    }
}

// 舒尔补计算
Matrix L = Kcc.choleskyDecomposition(Kcc);
Matrix Y = L.Solve_Matrix(L,Kca);
StiffnessMatrix = Kaa - Y.transpose() * Y;
}

// 对3x3矩阵c进行正交化
void orthogonalize(Matrix &c) {
    // Gram-Schmidt正交化
    // 取三列向量
    std::vector<double> v0(3), v1(3), v2(3);
    for(int i=0;i<3;++i) {
        v0[i] = c(i,0);
        v1[i] = c(i,1);
        v2[i] = c(i,2);
    }

    // 正交化v1
    double dot01 = v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
    for(int i=0;i<3;++i) v1[i] -= dot01 * v0[i];

    // 归一化v0,v1
    double norm0 = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);
    double norm1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    for(int i=0;i<3;++i) {
        v0[i] /= norm0;
        v1[i] /= norm1;
    }

    // v2正交化
    double dot02 = v0[0]*v2[0] + v0[1]*v2[1] + v0[2]*v2[2];
    double dot12 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    for(int i=0;i<3;++i) v2[i] -= dot02 * v0[i] + dot12 * v1[i];
    double norm2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    for(int i=0;i<3;++i) v2[i] /= norm2;

    // 回写
    for(int i=0;i<3;++i) {
        c(i,0) = v0[i];
        c(i,1) = v1[i];
        c(i,2) = v2[i];
    }
}

// 存储矩阵到文件
bool PipeSystem::saveMatrices(const Matrix& massMatrix, const Matrix& stiffnessMatrix, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "无法打开文件进行写入: " << filename << std::endl;
        return false;
    }
    
    try {
        // 写入质量矩阵
        int rows = massMatrix.rows;
        int cols = massMatrix.cols;
        file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double value = massMatrix(i, j);
                file.write(reinterpret_cast<const char*>(&value), sizeof(value));
            }
        }
        
        // 写入刚度矩阵
        rows = stiffnessMatrix.rows;
        cols = stiffnessMatrix.cols;
        file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
        
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double value = stiffnessMatrix(i, j);
                file.write(reinterpret_cast<const char*>(&value), sizeof(value));
            }
        }
        
        file.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "写入文件时出错: " << e.what() << std::endl;
        file.close();
        return false;
    }
}

// 从文件读取矩阵
bool PipeSystem::loadMatrices(Matrix& massMatrix, Matrix& stiffnessMatrix, const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "无法打开文件进行读取: " << filename << std::endl;
        return false;
    }
    
    try {
        // 读取质量矩阵
        int rows, cols;
        file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
        file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
        
        massMatrix = Matrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double value;
                file.read(reinterpret_cast<char*>(&value), sizeof(value));
                massMatrix(i, j) = value;
            }
        }
        
        // 读取刚度矩阵
        file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
        file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
        
        stiffnessMatrix = Matrix(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double value;
                file.read(reinterpret_cast<char*>(&value), sizeof(value));
                stiffnessMatrix(i, j) = value;
            }
        }
        
        file.close();
        return true;
    } catch (const std::exception& e) {
        std::cerr << "读取文件时出错: " << e.what() << std::endl;
        file.close();
        return false;
    }
}
