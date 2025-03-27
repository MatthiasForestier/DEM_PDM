#pragma once

#include <Eigen/Core>
#include "CRLHelper/VecMatDef.h"  // Ensure this defines MatrixXF

float computeConditionNumber(const MatrixXF &denseMatrix);