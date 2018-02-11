#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "beachtime.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(get_numeric_margins, 2),
    REGISTER(get_numeric_default_margins, 2),
    REGISTER(get_numeric_sparse_row_margins, 1),

    REGISTER(get_numeric_col_margins_const, 1),
    REGISTER(get_numeric_simple_row_margins, 1),

    REGISTER(get_numeric_arma_margins, 2),
    REGISTER(get_numeric_eigen_margins, 2),

    REGISTER(get_numeric_margins_random, 3),
    REGISTER(get_numeric_sparse_row_margins_random, 2),

    REGISTER(edit_numeric_matrix, 1),
    REGISTER(edit_numeric_arma_matrix, 1),

    REGISTER(dense_matrix_multiplication, 2),
    REGISTER(sparse_matrix_multiplication, 2),
    {NULL, NULL, 0}
};

void attribute_visible R_init_beachtime(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

