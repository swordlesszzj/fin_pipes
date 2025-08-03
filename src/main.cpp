#include <iostream>
#include "PipeSystem.h"
#include "Matrix.h"

int main() {
    PipeSystem pipeSystem;

    try {
        pipeSystem.loadFromFile("data/Default.txt");
        pipeSystem.assembleMatrices();

        Matrix temp(1, 1);
        auto result = temp.subspaceIteration(pipeSystem.StiffnessMatrix, pipeSystem.MassMatrix, 8);
        auto eigenvalues = result.first;
        auto eigenvectors = result.second;

        std::cout << "特征值: ";
        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            std::cout << eigenvalues[i];
            if (i < eigenvalues.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
