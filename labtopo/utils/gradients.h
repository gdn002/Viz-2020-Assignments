/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <labtopo/labtopomoduledefine.h>
#include <labtopo/utils/tnt/jama_eig.h>

namespace inviwo {
namespace util {

struct EigenResult
{
    mat2 eigenvectors;
    vec2 eigenvaluesRe;
    vec2 eigenvaluesIm;
};

EigenResult eigenAnalysis(mat2 matrix)
{
    EigenResult result;

    // Copy into tnt matrix.
    // WARNING: glm is column-major, tnt is row-major.
    TNT::Array2D<float> tntMat(2, 2);
    for (int c = 0; c < 2; ++c)
        for (int r = 0; r < 2; ++r)
            tntMat[c][r] = matrix[r][c];

    // Work magic.
    JAMA::Eigenvalue<float> eigs(tntMat);

    // Get values. Copy into glm vectors.
    TNT::Array1D<float> tntRe, tntIm;
    TNT::Array2D<float> tntEigVec;
    eigs.getRealEigenvalues(tntRe);
    eigs.getImagEigenvalues(tntIm);
    eigs.getV(tntEigVec);

    for (int dim = 0; dim < 2; ++dim)
    {
        result.eigenvaluesRe[dim] = tntRe[dim];
        result.eigenvaluesIm[dim] = tntIm[dim];

        for (int r = 0; r < 2; ++r)
            result.eigenvectors[dim][r] = tntEigVec[r][dim];
    }
    
    return result;
}

}// namespace
}