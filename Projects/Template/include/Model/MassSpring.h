#pragma once

#include "CRLHelper/VecMatDef.h"

class MassSpring {
   public:
    /// Fixed parameters
    F rest_length = 0.9;
    Vector2F endpoint0 = Vector2F(-1, 0);
    Vector2F endpoint1 = Vector2F(1, 0);

    /// Degrees of freedom
    Vector2F y = Vector2F(0, 0);

   public:
    void makeConfigMenu();

   public:
    void compute_energy(F &value) const;
    void compute_gradient(VectorXF &gradient) const;
    void compute_hessian(SparseMatrixF &hessian) const;
};
