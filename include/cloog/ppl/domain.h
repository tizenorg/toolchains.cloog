#ifndef CLOOG_PPL_DOMAIN_H
#define CLOOG_PPL_DOMAIN_H

#include <ppl_c.h>

#if defined(__cplusplus)
extern "C" 
  {
#endif 


struct cloogdomain {
	CloogState *state;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	int nb_par;
	int refs;
};

struct cloogscattering {
	struct cloogdomain dom;
};


#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
