#include <iostream>
#include "PipeSystem.h"
#include "Matrix.h"
#include "cblas.h"

int main()
{
    PipeSystem pipeSystem;

    try
    {
        pipeSystem.loadFromFile("data/simple2.txt");
        pipeSystem.assembleMatrices();
        // pipeSystem.saveMatrices(pipeSystem.MassMatrix,pipeSystem.StiffnessMatrix," matrix");
        // pipeSystem.loadMatrices(pipeSystem.MassMatrix,pipeSystem.StiffnessMatrix," matrix");
        Matrix temp(1, 1);
        // temp.choleskyDecomposition(pipeSystem.MassMatrix);

        auto result = temp.subspaceIteration(pipeSystem.StiffnessMatrix, pipeSystem.MassMatrix, 3);
        auto eigenvalues = result.first;
        auto eigenvectors = result.second;

        std::cout << "特征值: ";
        for (size_t i = 0; i < eigenvalues.size(); ++i)
        {
            std::cout << eigenvalues[i];
            if (i < eigenvalues.size() - 1)
            {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "错误: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
