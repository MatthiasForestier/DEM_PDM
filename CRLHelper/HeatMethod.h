#pragma once

#include "CRLHelper/VecMatDef.h"

namespace HeatMethod {

bool solveHeat(const MatrixXF &V, const MatrixXI &E, const VectorXI &B, VectorXF &sol);

bool mesh2D(const MatrixXF &V_in, const MatrixXI &E_in, MatrixXF &V_out, MatrixXI &T_out, VectorXI &B_out,
            F far_field_radius = 1.0, F max_area = -1.0);

}  // namespace HeatMethod
