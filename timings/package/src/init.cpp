#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "beachtime.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(get_numeric_margins, 2),
    REGISTER(get_numeric_default_margins, 2),
    REGISTER(get_numeric_random_margins, 3),
    REGISTER(get_numeric_arma_margins, 2),
    REGISTER(get_numeric_simple_margins, 1),
    {NULL, NULL, 0}
};

void attribute_visible R_init_beachtime(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}

