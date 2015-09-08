#include <cloog/ppl/cloog.h>

static int cloog_ppl_count = 0;

/**
 * Allocate and initialize full state.
 */
CloogState *cloog_state_malloc(void)
{
	CloogState *state;

	if (cloog_ppl_count == 0) {
		if (ppl_initialize() < 0)
			return NULL;
	}
	cloog_ppl_count++;

	state = cloog_core_state_malloc();
	state->backend = NULL;
	return state;
}

/**
 * Free state and backend independent parts.
 */
void cloog_state_free(CloogState *state)
{
	cloog_core_state_free(state);

	if (--cloog_ppl_count == 0)
	    ppl_finalize();
}
