#ifndef LIS_RCI_NATIVE_PROTOCOL_H
#define LIS_RCI_NATIVE_PROTOCOL_H
#include "lis_rci_status.h"
enum lis_rci_opaque_operation {
 LIS_RCI_OP_A=101, LIS_RCI_OP_B, LIS_RCI_OP_PRECON,
 LIS_RCI_OP_SHIFTED_SOLVE, LIS_RCI_OP_ZERO, LIS_RCI_OP_AXPBY
};
#endif
