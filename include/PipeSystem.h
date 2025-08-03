#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include "Pipe.h"
#include "Matrix.h"

using namespace std;

class PipeSystem {
public:
    vector<Pipe> pipes;
    unordered_map<int, int> number; // <节点编号, 顺序索引>
    vector<int> id;  // 节点编号列表，按顺序存储
    vector<Restraint> restraints;

    Matrix StiffnessMatrix;
    Matrix MassMatrix;

    PipeSystem();
    void loadFromFile(const string& filename);
    void assembleMatrices();

private:
    void resetTempMatrices();
};
