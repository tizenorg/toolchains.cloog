#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <cloog/ppl/cloog.h>

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n)*sizeof(type))


CloogDomain *cloog_domain_from_powerset(CloogState *state,
	ppl_Pointset_Powerset_C_Polyhedron_t ps, int nb_par)
{
	CloogDomain *domain;

	assert(ps);

	domain = ALLOC(CloogDomain);
	if (!domain)
		cloog_die("memory overflow.\n");

	domain->state = state;
	domain->ps = ps;
	domain->nb_par = nb_par;
	domain->refs = 1;

	return domain;
}


CloogScattering *cloog_scattering_from_powerset(CloogState *state,
	ppl_Pointset_Powerset_C_Polyhedron_t ps, int nb_par)
{
	return (CloogScattering *)cloog_domain_from_powerset(state, ps, nb_par);
}


/**
 * Returns true if each scattering dimension is or may be defined in terms
 * of the original iterators.
 */
int cloog_scattering_fully_specified(CloogScattering *scattering,
				      CloogDomain *domain)
{
	return 1;
}


/*
 * Return the number of constraints in the constraint system.
 */
static int ppl_Constraint_System_count(ppl_const_Constraint_System_t cs)
{
	int n = 0;
	ppl_Constraint_System_const_iterator_t cit, cit_end;

	if (ppl_new_Constraint_System_const_iterator(&cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Constraint_System_const_iterator(&cit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_begin(cs, cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_end(cs, cit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Constraint_System_const_iterator_equal_test(cit, cit_end)) {
		++n;

		if (ppl_Constraint_System_const_iterator_increment(cit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Constraint_System_const_iterator(cit_end);
	ppl_delete_Constraint_System_const_iterator(cit);

	return n;
}


/*
 * Create a CloogMatrix with the constraints of pol in PolyLib format.
 */
static CloogMatrix *cloog_matrix_from_polyhedron(ppl_const_Polyhedron_t pol)
{
	int row = 0;
	CloogMatrix *M;
	ppl_dimension_type dim;
	ppl_Coefficient_t coeff;
	ppl_const_Constraint_System_t cs;
	ppl_Constraint_System_const_iterator_t cit, cit_end;

	if (ppl_Polyhedron_get_minimized_constraints(pol, &cs) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint_System_const_iterator(&cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Constraint_System_const_iterator(&cit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_begin(cs, cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_end(cs, cit_end) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	ppl_Polyhedron_space_dimension(pol, &dim);

	M = cloog_matrix_alloc(ppl_Constraint_System_count(cs), 1 + dim + 1);

	while (!ppl_Constraint_System_const_iterator_equal_test(cit, cit_end)) {
		int col;
		int eq;
		ppl_const_Constraint_t con;

		if (ppl_Constraint_System_const_iterator_dereference(cit, &con) < 0)
			cloog_die("PPL error.\n");

		eq = ppl_Constraint_type(con) == PPL_CONSTRAINT_TYPE_EQUAL;

		cloog_int_set_si(M->p[row][0], eq ? 0 : 1);

		for (col = 0; col < dim; ++col) {
			if (ppl_Constraint_coefficient(con, col, coeff) < 0)
				cloog_die("PPL error.\n");
			if (ppl_Coefficient_to_mpz_t(coeff, M->p[row][1 + col]) < 0)
				cloog_die("PPL error.\n");
		}
		if (ppl_Constraint_inhomogeneous_term(con, coeff) < 0)
			cloog_die("PPL error.\n");
		if (ppl_Coefficient_to_mpz_t(coeff, M->p[row][1 + dim]) < 0)
			cloog_die("PPL error.\n");

		++row;

		if (ppl_Constraint_System_const_iterator_increment(cit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Coefficient(coeff);
	ppl_delete_Constraint_System_const_iterator(cit_end);
	ppl_delete_Constraint_System_const_iterator(cit);

	return M;
}


/*
 * Create a CloogMatrix with the constraints of ps, which is assumed
 * to contain exactly one disjunct, in PolyLib format.
 */
static CloogMatrix *cloog_matrix_from_powerset(ppl_Pointset_Powerset_C_Polyhedron_t ps)
{
	CloogMatrix *M;
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
		cloog_die("PPL error.\n");

	M = cloog_matrix_from_polyhedron(pol);

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);

	return M;
}


/**
 * Extract contraints from domain.
 */
CloogConstraintSet *cloog_domain_constraints(CloogDomain *domain)
{
	CloogMatrix *M;
	ppl_Pointset_Powerset_C_Polyhedron_t ps = domain->ps;
	size_t n;
	ppl_dimension_type dim;

	if (ppl_Pointset_Powerset_C_Polyhedron_omega_reduce(ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_size(ps, &n) < 0)
		cloog_die("PPL error.\n");
	assert(n <= 1);

	ppl_Pointset_Powerset_C_Polyhedron_space_dimension(ps, &dim);

	if (n == 0)
		M = cloog_matrix_alloc(0, 1 + dim + 1);
	else
		M = cloog_matrix_from_powerset(ps);

	return cloog_constraint_set_from_cloog_matrix(M);
}


void cloog_domain_print_constraints(FILE *foo, CloogDomain *domain,
					int print_number)
{
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;
	ppl_Pointset_Powerset_C_Polyhedron_t ps = domain->ps;
	size_t n;

	if (ppl_Pointset_Powerset_C_Polyhedron_omega_reduce(ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_size(ps, &n) < 0)
		cloog_die("PPL error.\n");

	if (print_number)
		fprintf(foo, "%d\n", (int)n);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(ps, pit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		ppl_const_Polyhedron_t pol;
		CloogMatrix *M;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");
		M = cloog_matrix_from_polyhedron(pol);
		cloog_matrix_print(foo, M);
		cloog_matrix_free(M);

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);
}

void cloog_scattering_print_constraints(FILE *foo, CloogScattering *scattering)
{
	cloog_domain_print_constraints(foo, &(scattering->dom), 1);
}

void cloog_domain_free(CloogDomain *domain)
{
	if (!domain)
		return;
	if (--domain->refs > 0)
		return;

	ppl_delete_Pointset_Powerset_C_Polyhedron(domain->ps);

	free(domain);
}


void cloog_scattering_free(CloogScattering *scatt)
{
	cloog_domain_free(&scatt->dom);
}


CloogDomain *cloog_domain_copy(CloogDomain *domain)
{
	if (!domain)
		return domain;

	domain->refs++;
	return domain;
}


CloogDomain *cloog_domain_duplicate(CloogDomain *domain)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");

	return cloog_domain_from_powerset(domain->state, ps, domain->nb_par);
}


/*
 * Create a new polyhedron that contains the convex hull of the polyhedra
 * is ps.
 */
static ppl_Polyhedron_t convex_hull(ppl_const_Pointset_Powerset_C_Polyhedron_t ps)
{
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;
	ppl_const_Polyhedron_t pol;
	ppl_Polyhedron_t hull;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(ps, pit_end) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_C_Polyhedron_from_C_Polyhedron(&hull, pol) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");
		if (ppl_Polyhedron_poly_hull_assign(hull, pol) < 0)
			cloog_die("PPL error.\n");
		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);

	return hull;
}


/**
 * cloog_domain_convex function:
 * Computes the convex hull of domain.
 */ 
CloogDomain *cloog_domain_convex(CloogDomain *domain)
{
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	size_t n;

	if (ppl_Pointset_Powerset_C_Polyhedron_omega_reduce(domain->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_size(domain->ps, &n) < 0)
		cloog_die("PPL error.\n");

	if (n <= 1)
		return cloog_domain_copy(domain);

	pol = convex_hull(domain->ps);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Polyhedron(pol);

	return cloog_domain_from_powerset(domain->state, ps, domain->nb_par);
}


/**
 * cloog_domain_simple_convex:
 * Given a list (union) of polyhedra, this function returns a "simple"
 * convex hull of this union.  For now, we just return the actual
 * convex hull.
 */
CloogDomain *cloog_domain_simple_convex(CloogDomain *domain)
{
	return cloog_domain_convex(domain);
}


/**
 * cloog_domain_simplify function:
 * Given two polyhedral domains (dom1) and (dom2),
 * this function finds the largest domain set (or the smallest list
 * of non-redundant constraints), that when intersected with polyhedral
 * domain (dom2) equals (dom1)intersect(dom2). The output is a new CloogDomain
 * structure with a polyhedral domain with the "redundant" constraints removed.
 * NB: the second domain is required not to be a union.
 */ 
CloogDomain *cloog_domain_simplify(CloogDomain *dom1, CloogDomain *dom2)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, dom1->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_simplify_using_context_assign(ps, dom2->ps) < 0)
		cloog_die("PPL error.\n");

	return cloog_domain_from_powerset(dom1->state, ps, dom1->nb_par);
}


/**
 * cloog_domain_union function:
 * This function returns a new polyhedral domain which is the union of
 * two polyhedral domains (dom1) U (dom2).
 */
CloogDomain *cloog_domain_union(CloogDomain *dom1, CloogDomain *dom2)
{
	CloogDomain *dom;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, dom1->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_upper_bound_assign(ps, dom2->ps) < 0)
		cloog_die("PPL error.\n");

	dom = cloog_domain_from_powerset(dom1->state, ps, dom1->nb_par);
	cloog_domain_free(dom1);
	cloog_domain_free(dom2);
	return dom;
}



/**
 * cloog_domain_intersection function:
 * This function returns a new polyhedral domain which is the intersection of
 * two polyhedral domains (dom1) \cap (dom2).
 */ 
CloogDomain *cloog_domain_intersection(CloogDomain *dom1, CloogDomain *dom2)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, dom1->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_intersection_assign(ps, dom2->ps) < 0)
		cloog_die("PPL error.\n");

	return cloog_domain_from_powerset(dom1->state, ps, dom1->nb_par);
}


/*
 * Given an equality "a x + b = 0" and a powerset U,
 * construct a new powerset V equal to
 *
 * (U \cap { a x + b >= 1 }) \cup (U \cap { a x + b <= -1 })
 */
static ppl_Pointset_Powerset_C_Polyhedron_t subtract_eq(
	ppl_const_Pointset_Powerset_C_Polyhedron_t ps, ppl_const_Constraint_t eq)
{
	ppl_Pointset_Powerset_C_Polyhedron_t res1, res2;
	ppl_Linear_Expression_t le, le_eq;
	ppl_Constraint_t con;
	ppl_Coefficient_t coeff;
	cloog_int_t v;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&res1, ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&res2, ps) < 0)
		cloog_die("PPL error.\n");

	cloog_int_init(v);
	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	cloog_int_set_si(v, -1);
	ppl_assign_Coefficient_from_mpz_t(coeff, v);

	if (ppl_new_Linear_Expression(&le) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Linear_Expression_from_Constraint(&le_eq, eq) < 0)
		cloog_die("PPL error.\n");
	if (ppl_subtract_Linear_Expression_from_Linear_Expression(le, le_eq) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_add_constraint(res1, con) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Constraint(con);

	if (ppl_assign_Linear_Expression_from_Linear_Expression(le, le_eq) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_add_constraint(res2, con) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Constraint(con);

	ppl_delete_Linear_Expression(le);
	ppl_delete_Linear_Expression(le_eq);

	ppl_delete_Coefficient(coeff);
	cloog_int_clear(v);

	if (ppl_Pointset_Powerset_C_Polyhedron_upper_bound_assign(res1, res2) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Pointset_Powerset_C_Polyhedron(res2);

	return res1;
}


/*
 * Given an inequality "a x + b >= 0" and a powerset U,
 * construct a new powerset V equal to U \cap { a x + b <= -1 }
 */
static ppl_Pointset_Powerset_C_Polyhedron_t subtract_ineq(
	ppl_const_Pointset_Powerset_C_Polyhedron_t ps, ppl_const_Constraint_t ineq)
{
	ppl_Pointset_Powerset_C_Polyhedron_t res;
	ppl_Linear_Expression_t le, le_ineq;
	ppl_Constraint_t con;
	ppl_Coefficient_t coeff;
	cloog_int_t v;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&res, ps) < 0)
		cloog_die("PPL error.\n");

	cloog_int_init(v);
	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	cloog_int_set_si(v, -1);
	ppl_assign_Coefficient_from_mpz_t(coeff, v);

	if (ppl_new_Linear_Expression(&le) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Linear_Expression_from_Constraint(&le_ineq, ineq) < 0)
		cloog_die("PPL error.\n");
	if (ppl_subtract_Linear_Expression_from_Linear_Expression(le, le_ineq) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_add_constraint(res, con) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Constraint(con);
	ppl_delete_Linear_Expression(le);
	ppl_delete_Linear_Expression(le_ineq);

	ppl_delete_Coefficient(coeff);
	cloog_int_clear(v);

	return res;
}


/*
 * Update the powerset ps to ps \ pol by intersecting the original ps
 * with the negation of each of the constraints of pol and collecting
 * the results.
 */
static ppl_Pointset_Powerset_C_Polyhedron_t subtract(
	ppl_Pointset_Powerset_C_Polyhedron_t ps, ppl_const_Polyhedron_t pol)
{
	ppl_dimension_type dim;
	ppl_Pointset_Powerset_C_Polyhedron_t res;
	ppl_const_Constraint_System_t cs;
	ppl_Constraint_System_const_iterator_t cit, cit_end;

	ppl_Pointset_Powerset_C_Polyhedron_space_dimension(ps, &dim);
	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_space_dimension(&res, dim, 1) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Polyhedron_get_minimized_constraints(pol, &cs) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint_System_const_iterator(&cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Constraint_System_const_iterator(&cit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_begin(cs, cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_end(cs, cit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Constraint_System_const_iterator_equal_test(cit, cit_end)) {
		ppl_const_Constraint_t con;
		ppl_Pointset_Powerset_C_Polyhedron_t d;
		int eq;

		if (ppl_Constraint_System_const_iterator_dereference(cit, &con) < 0)
			cloog_die("PPL error.\n");

		eq = ppl_Constraint_type(con) == PPL_CONSTRAINT_TYPE_EQUAL;
		if (eq)
			d = subtract_eq(ps, con);
		else
			d = subtract_ineq(ps, con);

		if (ppl_Pointset_Powerset_C_Polyhedron_upper_bound_assign(res, d) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Pointset_Powerset_C_Polyhedron(d);

		if (ppl_Constraint_System_const_iterator_increment(cit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Constraint_System_const_iterator(cit_end);
	ppl_delete_Constraint_System_const_iterator(cit);

	ppl_delete_Pointset_Powerset_C_Polyhedron(ps);

	return res;
}


/**
 * cloog_domain_difference function:
 * Returns the set difference domain \ minus.
 */ 
CloogDomain *cloog_domain_difference(CloogDomain *domain, CloogDomain *minus)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(minus->ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(minus->ps, pit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		ppl_const_Polyhedron_t pol;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");

		ps = subtract(ps, pol);

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);

	return cloog_domain_from_powerset(domain->state, ps, domain->nb_par);
}


/*
 * Update powerset ps by inserting n unconstrained dimensions before
 * dimension pos.
 */
static void powerset_insert_dims(ppl_Pointset_Powerset_C_Polyhedron_t ps,
	int pos, int n)
{
	ppl_dimension_type dim;
	ppl_dimension_type *ds;
	int i;

	ppl_Pointset_Powerset_C_Polyhedron_space_dimension(ps, &dim);
	if (ppl_Pointset_Powerset_C_Polyhedron_add_space_dimensions_and_embed(ps, n) < 0)
		cloog_die("PPL error.\n");

	ds = ALLOCN(ppl_dimension_type, dim + n);
	if (!ds)
		cloog_die("memory overflow.\n");

	for (i = 0; i < pos; ++i)
		ds[i] = i;
	for (i = 0; i < dim - pos; ++i)
		ds[pos + i] = pos + n + i;
	for (i = 0; i < n; ++i)
		ds[dim + i] = pos + i;

	if (ppl_Pointset_Powerset_C_Polyhedron_map_space_dimensions(ps, ds, dim + n) < 0)
		cloog_die("PPL error.\n");

	free(ds);
}


/*
 * Given two disjoint powersets ps1 and ps2 in a space of dimension d,
 * check whether for any common value of the parameters and dimensions
 * preceding pos in both powersets, the values of dimension pos in ps1 are
 * smaller or larger than those in ps2.
 *
 * Return
 *	 1 if bset1 follows bset2
 *	-1 if bset1 precedes bset2
 *	 0 if bset1 and bset2 are incomparable
 */
static int powerset_compare_at(ppl_const_Pointset_Powerset_C_Polyhedron_t ps1,
	ppl_const_Pointset_Powerset_C_Polyhedron_t ps2, int pos, int d)
{
	int cmp;
	int bounded, dummy;
	ppl_Pointset_Powerset_C_Polyhedron_t ps1c, ps2c;
	ppl_Linear_Expression_t le;
	ppl_dimension_type sdim;
	ppl_Coefficient_t coeff_n;
	ppl_Coefficient_t coeff_d;
	cloog_int_t v;

	cloog_int_init(v);
	if (ppl_new_Coefficient(&coeff_d) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Coefficient(&coeff_n) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps1c, ps1) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps2c, ps2) < 0)
		cloog_die("PPL error.\n");
	powerset_insert_dims(ps1c, d, d - pos);
	powerset_insert_dims(ps2c, pos, d - pos);
	if (ppl_Pointset_Powerset_C_Polyhedron_intersection_assign(ps1c, ps2c) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_is_empty(ps1c)) {
		cmp = 0;
	} else {
		ppl_Pointset_Powerset_C_Polyhedron_space_dimension(ps1c, &sdim);
		if (ppl_new_Linear_Expression_with_dimension(&le, sdim) < 0)
			cloog_die("PPL error.\n");

		cloog_int_set_si(v, 1);
		ppl_assign_Coefficient_from_mpz_t(coeff_n, v);
		if (ppl_Linear_Expression_add_to_coefficient(le, pos, coeff_n) < 0)
			cloog_die("PPL error.\n");

		cloog_int_set_si(v, -1);
		ppl_assign_Coefficient_from_mpz_t(coeff_n, v);
		if (ppl_Linear_Expression_add_to_coefficient(le, d, coeff_n) < 0)
			cloog_die("PPL error.\n");

		bounded = ppl_Pointset_Powerset_C_Polyhedron_maximize(ps1c, le, coeff_n, coeff_d, &dummy);
		if (!bounded)
			cmp = 1;
		else {
			if (ppl_Coefficient_to_mpz_t(coeff_n, v) < 0)
				cloog_die("PPL error.\n");
			if (cloog_int_is_pos(v))
				cmp = 1;
			else
				cmp = -1;
		}

		ppl_delete_Linear_Expression(le);
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron(ps1c);
	ppl_delete_Pointset_Powerset_C_Polyhedron(ps2c);

	ppl_delete_Coefficient(coeff_n);
	ppl_delete_Coefficient(coeff_d);
	cloog_int_clear(v);

	return cmp;
}


/**
 * cloog_domain_sort function:
 * This function topologically sorts (nb_doms) domains. Here (doms) is an
 * array of pointers to CloogDomains, (nb_doms) is the number of domains,
 * (level) is the level to consider for partial ordering, (permut) if not
 * NULL, is an array of (nb_doms) integers that contains a permutation
 * specification after call in order to apply the topological sorting. 
 */
void cloog_domain_sort(CloogDomain **doms, unsigned nb_doms, unsigned level,
			int *permut)
{
	int i, j, k, cmp;
	unsigned char **follows;
	int d;

	if (!nb_doms)
		return;
	d = cloog_domain_dimension(doms[0]);

	follows = ALLOCN(unsigned char *, nb_doms);
	if (!follows)
		cloog_die("memory overflow.\n");
	for (i = 0; i < nb_doms; ++i) {
		follows[i] = ALLOCN(unsigned char, nb_doms);
		if (!follows[i])
			cloog_die("memory overflow.\n");
		for (j = 0; j < nb_doms; ++j)
			follows[i][j] = 0;
	}

	for (i = 1; i < nb_doms; ++i) {
		for (j = 0; j < i; ++j) {
			if (follows[i][j] || follows[j][i])
				continue;
			cmp = powerset_compare_at(doms[i]->ps, doms[j]->ps,
						  level - 1, d);
			if (!cmp)
				continue;
			if (cmp > 0) {
				follows[i][j] = 1;
				for (k = 0; k < i; ++k)
					follows[i][k] |= follows[j][k];
			} else {
				follows[j][i] = 1;
				for (k = 0; k < i; ++k)
					follows[k][i] |= follows[k][j];
			}
		}
	}

	for (i = 0, j = 0; i < nb_doms; j = (j + 1) % nb_doms) {
		for (k = 0; k < nb_doms; ++k)
			if (follows[j][k])
				break;
		if (k < nb_doms)
			continue;
		for (k = 0; k < nb_doms; ++k)
			follows[k][j] = 0;
		follows[j][j] = 1;
		permut[i] = 1 + j;
		++i;
	}

	for (i = 0; i < nb_doms; ++i)
		free(follows[i]);
	free(follows);
}


/**
 * Check whether there is or may be any value of dom1 at the given level
 * that is greater than or equal to a value of dom2 at the same level.
 *
 * Return
 *	 1 is there is or may be a greater-than pair.
 *	 0 if there is no greater-than pair, but there may be an equal-to pair
 *	-1 if there is definitely no such pair
 */
int cloog_domain_follows(CloogDomain *dom1, CloogDomain *dom2, unsigned level)
{
	return 1;
}


/**
 * cloog_domain_empty function:
 * Returns an empty domain of the same dimensions as template.
 */
CloogDomain *cloog_domain_empty(CloogDomain *template)
{
	int dim = cloog_domain_dimension(template) +
		  cloog_domain_parameter_dimension(template);
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_space_dimension(&ps, dim, 1) < 0)
		cloog_die("PPL error.\n");

	return cloog_domain_from_powerset(template->state, ps, template->nb_par);
}


static int polyhedron_is_bounded(ppl_const_Polyhedron_t pol, unsigned level)
{
	int lower = 0, upper = 0;
	ppl_Coefficient_t coeff;
	ppl_const_Constraint_System_t cs;
	ppl_Constraint_System_const_iterator_t cit, cit_end;
	mpz_t v;

	mpz_init(v);

	if (ppl_Polyhedron_get_minimized_constraints(pol, &cs) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint_System_const_iterator(&cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Constraint_System_const_iterator(&cit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_begin(cs, cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_end(cs, cit_end) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Constraint_System_const_iterator_equal_test(cit, cit_end)) {
		ppl_const_Constraint_t con;

		if (ppl_Constraint_System_const_iterator_dereference(cit, &con) < 0)
			cloog_die("PPL error.\n");

		if (ppl_Constraint_coefficient(con, level - 1, coeff) < 0)
			cloog_die("PPL error.\n");
		if (ppl_Coefficient_to_mpz_t(coeff, v) < 0)
			cloog_die("PPL error.\n");

		if (mpz_sgn(v) != 0) {
			if (ppl_Constraint_type(con) == PPL_CONSTRAINT_TYPE_EQUAL) {
				lower = upper = 1;
				break;
			}
			if (mpz_sgn(v) > 0)
				lower = 1;
			else
				upper = 1;
			if (lower && upper)
				break;
		}

		if (ppl_Constraint_System_const_iterator_increment(cit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Coefficient(coeff);
	ppl_delete_Constraint_System_const_iterator(cit_end);
	ppl_delete_Constraint_System_const_iterator(cit);

	mpz_clear(v);

	return lower && upper;
}


/**
 * Return 1 if the specified dimension has both an upper and a lower bound.
 */
int cloog_domain_is_bounded(CloogDomain *dom, unsigned level)
{
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;
	ppl_Pointset_Powerset_C_Polyhedron_t ps = dom->ps;
	int bounded = 1;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(ps, pit_end) < 0)
		cloog_die("PPL error.\n");
	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		ppl_const_Polyhedron_t pol;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");

		if (!polyhedron_is_bounded(pol, level)) {
			bounded = 0;
			break;
		}

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);

	return bounded;
}


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_domain_print_structure :
 * this function is a more human-friendly way to display the CloogDomain data
 * structure, it only shows the constraint system and includes an indentation
 * level (level) in order to work with others print_structure functions.
 */
void cloog_domain_print_structure(FILE *file, CloogDomain *domain, int level,
				  const char *name)
{
	int i;
	char *suffix = " ]";
	char *prefix;
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;
	ppl_Pointset_Powerset_C_Polyhedron_t ps = domain->ps;

	for (i = 0; i < level; i++)
		fprintf(file, "|\t");
  
	if (!domain->ps) {
		fprintf(file, "+-- Null CloogDomain\n");
		return;
	}
	fprintf(file, "+-- %s\n", name);
	prefix = ALLOCN(char, 2 * (level + 1) + 3);
	if (!prefix)
		cloog_die("memory overflow.\n");
	for (i = 0; i < level + 1; ++i)
		memcpy(prefix + 2 * i, "|\t", 2);
	strcpy(prefix + 2 * (level + 1), "[ ");

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(ps, pit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		ppl_const_Polyhedron_t pol;
		CloogMatrix *M;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");

		M = cloog_matrix_from_polyhedron(pol);
		cloog_matrix_print_structure(file, M, prefix, suffix);
		cloog_matrix_free(M);

		prefix[2*(level+1)] = '\0';
		fprintf(file, "%s\n", prefix);
		prefix[2*(level+1)] = '[';

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	free(prefix);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


void cloog_domain_list_free(CloogDomainList *list)
{
	CloogDomainList *next;

	for ( ; list; list = next) {
		next = list->next;
		cloog_domain_free(list->domain);
		free(list);
	}
}


/**
 * cloog_scattering_list_free function:
 * This function frees the allocated memory for a CloogScatteringList structure.
 */
void cloog_scattering_list_free(CloogScatteringList *list)
{
	while (list != NULL) {
		CloogScatteringList *temp = list->next;
		cloog_scattering_free(list->scatt);
		free(list);
		list = temp;
	}
}


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/


/*
 * Construct a polyhedron from the constraint in M (PolyLib format).
 */
static ppl_Polyhedron_t cloog_matrix_to_polyhedron(CloogMatrix *M)
{
	int row;
	unsigned dim;
	ppl_Polyhedron_t pol;
	ppl_Coefficient_t coeff;
	cloog_int_t v;

	if (M->NbColumns < 2)
		cloog_die("Input error.\n");
	dim = M->NbColumns - 2;
	if (ppl_new_C_Polyhedron_from_space_dimension(&pol, dim, 0) < 0)
		cloog_die("PPL error.\n");

	cloog_int_init(v);

	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	for (row = 0; row < M->NbRows; ++row) {
		int i;
		ppl_Linear_Expression_t le;
		ppl_Constraint_t con;
		int eq = cloog_int_is_zero(M->p[row][0]);

		if (ppl_new_Linear_Expression_with_dimension(&le, dim) < 0)
			cloog_die("PPL error.\n");

		for (i = 0; i < dim; ++i) {
			ppl_assign_Coefficient_from_mpz_t(coeff, M->p[row][1 + i]);

			if (ppl_Linear_Expression_add_to_coefficient(le, i, coeff) < 0)
				cloog_die("PPL error.\n");
		}

		ppl_assign_Coefficient_from_mpz_t(coeff, M->p[row][1 + dim]);
		if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
			cloog_die("PPL error.\n");

		if (ppl_new_Constraint(&con, le,
				eq ? PPL_CONSTRAINT_TYPE_EQUAL :
				     PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Linear_Expression(le);

		if (ppl_Polyhedron_add_constraint(pol, con) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Constraint(con);
	}

	ppl_delete_Coefficient(coeff);
	cloog_int_clear(v);

	return pol;
}


static char *next_line(FILE *input, char *line, unsigned len)
{
	char *p;

	do {
		if (!(p = fgets(line, len, input)))
			return NULL;
		while (isspace(*p) && *p != '\n')
			++p;
	} while (*p == '#' || *p == '\n');

	return p;
}


static ppl_Pointset_Powerset_C_Polyhedron_t
read_more_disjuncts(FILE *input, ppl_Pointset_Powerset_C_Polyhedron_t ps, int n)
{
	CloogMatrix *M;
	ppl_const_Polyhedron_t pol;

	while (n-- > 0) {
		M = cloog_matrix_read(input);
		pol = cloog_matrix_to_polyhedron(M);
		cloog_matrix_free(M);

		if (ppl_Pointset_Powerset_C_Polyhedron_add_disjunct(ps, pol) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Polyhedron(pol);
	}

	return ps;
}


/**
 * cloog_domain_read_context function:
 * Read parameter domain.
 */
CloogDomain *cloog_domain_read_context(CloogState *state, FILE *input)
{
	char line[1024];
	CloogMatrix *M;
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	int nb_par;
	unsigned n_row, n_col;
	int n = 1;

	if (!next_line(input, line, sizeof(line)))
		cloog_die("Input error.\n");
	if (sscanf(line, "%u %u", &n_row, &n_col) == 2)
		M = cloog_matrix_read_of_size(input, n_row, n_col);
	else {
		if (sscanf(line, "%d", &n) != 1)
			cloog_die("Input error.\n");
		if (n < 1)
			cloog_die("Input error.\n");
		M = cloog_matrix_read(input);
	}

	pol = cloog_matrix_to_polyhedron(M);
	nb_par = M->NbColumns - 2;
	cloog_matrix_free(M);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Polyhedron(pol);

	ps = read_more_disjuncts(input, ps, n - 1);

	return cloog_domain_from_powerset(state, ps, nb_par);
}


/**
 * cloog_domain_from_context
 * Reinterpret context by turning parameters into variables.
 */
CloogDomain *cloog_domain_from_context(CloogDomain *context)
{
	CloogDomain *domain;

	domain = cloog_domain_duplicate(context);
	cloog_domain_free(context);

	domain->nb_par = 0;

	return domain;
}


/**
 * cloog_domain_union_read function:
 * This function reads a union of polyhedra into a file (input) and
 * returns a pointer to a CloogDomain containing the read information. 
 */
CloogDomain *cloog_domain_union_read(CloogState *state,
					FILE *input, int nb_parameters)
{
	char line[1024];
	int n;
	CloogMatrix *M;
	ppl_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (!next_line(input, line, sizeof(line)))
		cloog_die("Input error.\n");
	if (sscanf(line, "%d", &n) != 1)
		cloog_die("Input error.\n");
	if (n < 1)
		cloog_die("Input error.\n");

	M = cloog_matrix_read(input);
	pol = cloog_matrix_to_polyhedron(M);
	cloog_matrix_free(M);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");
	ppl_delete_Polyhedron(pol);

	ps = read_more_disjuncts(input, ps, n - 1);

	return cloog_domain_from_powerset(state, ps, nb_parameters);
}


/**
 * cloog_domain_read_scattering function:
 * This function reads in a scattering function from the file input.
 */
CloogScattering *cloog_domain_read_scattering(CloogDomain *domain, FILE *input)
{
	char line[1024];
	CloogMatrix *M;
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	unsigned n_row, n_col;
	int n = 1;

	if (!next_line(input, line, sizeof(line)))
		cloog_die("Input error.\n");
	if (sscanf(line, "%u %u", &n_row, &n_col) == 2)
		M = cloog_matrix_read_of_size(input, n_row, n_col);
	else {
		if (sscanf(line, "%d", &n) != 1)
			cloog_die("Input error.\n");
		if (n < 1)
			cloog_die("Input error.\n");
		M = cloog_matrix_read(input);
	}

	pol = cloog_matrix_to_polyhedron(M);
	cloog_matrix_free(M);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Polyhedron(pol);

	ps = read_more_disjuncts(input, ps, n - 1);

	return cloog_scattering_from_powerset(domain->state, ps, domain->nb_par);
}

/******************************************************************************
 *                      CloogMatrix Reading function                          *
 ******************************************************************************/

/**
 * cloog_domain_from_cloog_matrix:
 * Create a CloogDomain containing the constraints described in matrix.
 * nparam is the number of parameters contained in the domain.
 * Returns a pointer to the CloogDomain if successful; NULL otherwise.
 */
CloogDomain *cloog_domain_from_cloog_matrix(CloogState *state,
	CloogMatrix *matrix, int nparam)
{
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	pol = cloog_matrix_to_polyhedron(matrix);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Polyhedron(pol);

	return cloog_domain_from_powerset(state, ps, nparam);
}

/**
 * cloog_scattering_from_cloog_matrix:
 * Create a CloogScattering containing the constraints described in matrix.
 * nparam is the number of parameters contained in the domain.
 * Returns a pointer to the CloogScattering if successful; NULL otherwise.
 */
CloogScattering *cloog_scattering_from_cloog_matrix(CloogState *state,
	CloogMatrix *matrix, int nb_scat, int nb_par)
{
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	pol = cloog_matrix_to_polyhedron(matrix);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Polyhedron(pol);

	return cloog_scattering_from_powerset(state, ps, nb_par);
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/



/**
 * cloog_domain_isempty function:
 */ 
int cloog_domain_isempty(CloogDomain *domain)
{
	return ppl_Pointset_Powerset_C_Polyhedron_is_empty(domain->ps);
}


/**
 * cloog_domain_universe function:
 * This function returns the complete dim-dimensional space.
 */
CloogDomain *cloog_domain_universe(CloogState *state, unsigned dim)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_space_dimension(&ps, dim, 0) < 0)
		cloog_die("PPL error.\n");

	return cloog_domain_from_powerset(state, ps, 0);
}


/**
 * cloog_domain_project function:
 * This function returns the projection of
 * (domain) on the (level) first dimensions (i.e. outer loops).
 */ 
CloogDomain *cloog_domain_project(CloogDomain *domain, int level)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	int n = cloog_domain_dimension(domain) - level;
	ppl_dimension_type *ds;
	int i;

	ds = ALLOCN(ppl_dimension_type, n);
	if (!ds)
		cloog_die("memory overflow.\n");

	for (i = 0; i < n; ++i)
		ds[i] = level + i;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_remove_space_dimensions(ps, ds, n) < 0)
		cloog_die("PPL error.\n");

	free(ds);

	return cloog_domain_from_powerset(domain->state, ps, domain->nb_par);
}


/**
 * cloog_domain_extend function:
 * This function returns the (domain) given as input with (dim)
 * dimensions and (nb_par) parameters.
 * This function does not free (domain), and returns a new CloogDomain.
 */ 
CloogDomain *cloog_domain_extend(CloogDomain *domain, int dim)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	int n = dim - cloog_domain_dimension(domain);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");
	powerset_insert_dims(ps, cloog_domain_dimension(domain), n);
	return cloog_domain_from_powerset(domain->state, ps, domain->nb_par);
}


/*
 * Check if the given polyhedron lies on any hyperplane that doesn't
 * contain any integer point.
 */
static int polyhedron_never_integral(ppl_const_Polyhedron_t pol)
{
	int empty = 0;
	CloogMatrix *M;
	int i;
	cloog_int_t gcd;

	cloog_int_init(gcd);
	M = cloog_matrix_from_polyhedron(pol);
	for (i = 0; i < M->NbRows; ++i) {
		if (!cloog_int_is_zero(M->p[i][0]))
			continue;
		if (cloog_int_is_zero(M->p[i][M->NbColumns - 1]))
			continue;
		cloog_seq_gcd(M->p[i] + 1, M->NbColumns - 2, &gcd);
		cloog_int_fdiv_r(gcd, M->p[i][M->NbColumns - 1], gcd);
		if (cloog_int_is_zero(gcd))
			continue;
		empty = 1;
		break;
	}
	cloog_matrix_free(M);
	cloog_int_clear(gcd);

	return empty;
}


/**
 * cloog_domain_never_integral function:
 * For us, an equality like 3*i -4 = 0 is always false since 4%3 != 0.
 */
int cloog_domain_never_integral(CloogDomain *domain)
{
	int empty = 1;
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;
	ppl_Pointset_Powerset_C_Polyhedron_t ps = domain->ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(ps, pit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		ppl_const_Polyhedron_t pol;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");

		if (!polyhedron_never_integral(pol)) {
			empty = 0;
			break;
		}

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);

	return empty;
}


/*
 * Check if the given equality constraint imposes a stride on
 * the iterator i identified by strided_level.  If so, set *stride
 * to the stride and *offset to a value c such that (i + c) % stride = 0.
 * Otherwise, leave *stride alone.
 */
static void equality_stride(ppl_const_Constraint_t eq, int dim,
	int strided_level, cloog_int_t *stride, cloog_int_t *offset)
{
	int i;
	ppl_Coefficient_t coeff;
	cloog_int_t v;
	cloog_int_t cst;
	cloog_int_t gcd;
	cloog_int_t tmp;

	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");
	cloog_int_init(v);
	cloog_int_init(cst);
	cloog_int_init(gcd);
	cloog_int_init(tmp);

	if (ppl_Constraint_coefficient(eq, strided_level - 1, coeff) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Coefficient_to_mpz_t(coeff, v) < 0)
		cloog_die("PPL error.\n");

	if (cloog_int_is_zero(v))
		goto done;

	if (ppl_Constraint_inhomogeneous_term(eq, coeff) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Coefficient_to_mpz_t(coeff, cst) < 0)
		cloog_die("PPL error.\n");

	/* Too complicated for our quick heuristic. */
	if (!cloog_int_is_divisible_by(cst, v))
		goto done;

	cloog_int_divexact(cst, cst, v);

	cloog_int_set_si(gcd, 0);
	for (i = strided_level; i < dim; ++i) {
		if (ppl_Constraint_coefficient(eq, i, coeff) < 0)
			cloog_die("PPL error.\n");
		if (ppl_Coefficient_to_mpz_t(coeff, tmp) < 0)
			cloog_die("PPL error.\n");
		if (!cloog_int_is_zero(tmp))
			cloog_int_gcd(gcd, gcd, tmp);
	}

	if (cloog_int_is_zero(gcd) || cloog_int_is_one(gcd))
		goto done;

	cloog_int_gcd(tmp, gcd, v);
	cloog_int_divexact(gcd, gcd, tmp);

	cloog_int_set(*stride, gcd);
	cloog_int_fdiv_r(*offset, cst, gcd);
done:
	cloog_int_clear(v);
	cloog_int_clear(cst);
	cloog_int_clear(gcd);
	cloog_int_clear(tmp);
	ppl_delete_Coefficient(coeff);
}


/*
 * Check if any of the equalities of pol imposes a stride on
 * the iterator i identified by strided_level.  If so, set *stride
 * to the stride and *offset to a value c such that (i + c) % stride = 0.
 * Otherwise, set *stride to 1 and *offset to 0.
 */
static void polyhedron_stride(ppl_const_Polyhedron_t pol, int dim,
	int strided_level, cloog_int_t *stride, cloog_int_t *offset)
{
	ppl_const_Constraint_System_t cs;
	ppl_Constraint_System_const_iterator_t cit, cit_end;

	cloog_int_set_si(*stride, 0);
	cloog_int_set_si(*offset, 0);

	if (ppl_Polyhedron_get_minimized_constraints(pol, &cs) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Constraint_System_const_iterator(&cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Constraint_System_const_iterator(&cit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_begin(cs, cit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Constraint_System_end(cs, cit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Constraint_System_const_iterator_equal_test(cit, cit_end)) {
		ppl_const_Constraint_t con;

		if (ppl_Constraint_System_const_iterator_dereference(cit, &con) < 0)
			cloog_die("PPL error.\n");

		if (ppl_Constraint_type(con) == PPL_CONSTRAINT_TYPE_EQUAL)
			equality_stride(con, dim, strided_level, stride, offset);
		if (!cloog_int_is_zero(*stride))
			break;

		if (ppl_Constraint_System_const_iterator_increment(cit) < 0)
			cloog_die("PPL error.\n");
	}

	if (cloog_int_is_zero(*stride))
		cloog_int_set_si(*stride, 1);

	ppl_delete_Constraint_System_const_iterator(cit_end);
	ppl_delete_Constraint_System_const_iterator(cit);
}


/**
 * Check whether the loop at "level" is executed at most once.
 * We conservatively assume that the loop may be executed more than once.
 */
int cloog_domain_is_otl(CloogDomain *domain, int level)
{
	return 0;
}


/**
 * cloog_domain_stride function:
 * This function finds the stride imposed to unknown with the column number
 * 'strided_level' in order to be integral. For instance, if we have a
 * constraint like -i - 2j + 2k = 0, and we consider k, then k can be integral
 * only if (i + 2j)%2 = 0. Then only if i%2 = 0. Then k imposes a stride 2 to
 * the unknown i. The function returns the imposed stride in a parameter field.
 * - domain is the set of constraint we have to consider,
 * - strided_level is the column number of the unknown for which a stride have
 *   to be found,
 * - stride is the stride that is returned back as a function parameter. 
 * - offset is the value of the constant c if the condition is of the shape
 *   (i + c)%s = 0, s being the stride.
 *
 * We only look for simple cases.  In particular, we don't even bother
 * if the domain is a union and we only look for strides imposed by
 * a single equality and not any combination of the equalities.
 */
void cloog_domain_stride(CloogDomain *domain, int strided_level,
	cloog_int_t *stride, cloog_int_t *offset)
{
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit;
	int dim = cloog_domain_dimension(domain);
	size_t n;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(domain->ps, pit) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_size(domain->ps, &n) < 0)
		cloog_die("PPL error.\n");

	if (n != 1)
		goto no_stride;

	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
		cloog_die("PPL error.\n");

	polyhedron_stride(pol, dim, strided_level, stride, offset);

	if (0) {
no_stride:
		cloog_int_set_si(*stride, 1);
		cloog_int_set_si(*offset, 0);
	}
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
}


/**
 * Return 1 if CLooG is allowed to perform stride detection on level "level"
 * and 0 otherwise.
 * In particular, stride detection should only be performed when the lower
 * bound at the given level is an integral constant.
 */
int cloog_domain_can_stride(CloogDomain *domain, int level)
{
	int found = 0;
	int row;
	CloogMatrix *M;
	size_t n;
	int dim = cloog_domain_dimension(domain);
	int nb_par = cloog_domain_parameter_dimension(domain);

	if (ppl_Pointset_Powerset_C_Polyhedron_size(domain->ps, &n) < 0)
		cloog_die("PPL error.\n");

	if (n != 1)
		return 0;

	M = cloog_matrix_from_powerset(domain->ps);

	for (row = 0; row < M->NbRows; ++row) {
		if (cloog_int_is_zero(M->p[row][level]))
			continue;
		if (cloog_int_is_zero(M->p[row][0]))
			break;
		if (cloog_int_is_neg(M->p[row][level]))
			continue;
		if (found)
			break;
		if (cloog_seq_first_non_zero(M->p[row] + 1, level - 1) != -1)
			break;
		if (cloog_seq_first_non_zero(M->p[row] + level + 1,
					      dim - level + nb_par) != -1)
			break;
		found = 1;
	}
	if (row < M->NbRows)
		found = 0;

	cloog_matrix_free(M);

	return found;
}


/**
 * Update the lower bounds at level "level" to the given stride information.
 * That is, make sure that the remainder on division by "stride"
 * is equal to "offset".
 * Since we only allow stride detection on loops with a fixed integral
 * lower bound, we don't actually need to change domain here as the
 * fixed lower bound will be updated by update_lower_bound in clast.c.
 */
CloogDomain *cloog_domain_stride_lower_bound(CloogDomain *domain, int level,
	CloogStride *stride)
{
	return domain;
}


/**
 * cloog_domain_lazy_equal function:
 * This function returns 1 if the domains given as input are the same, 0 if it
 * is unable to decide.
 */
int cloog_domain_lazy_equal(CloogDomain *d1, CloogDomain *d2)
{
	return ppl_Pointset_Powerset_C_Polyhedron_equals_Pointset_Powerset_C_Polyhedron(d1->ps, d2->ps);
}


/**
 * Return a union of sets S_i such that the convex hull of "dom",
 * when intersected with one the sets S_i, will have an upper and
 * lower bound for the dimension at "level" (provided "dom" itself
 * has such bounds for the dimensions).
 *
 * We currently take a very simple approach and split all parameters
 * into a negative and a positive part.
 */
CloogDomain *cloog_domain_bound_splitter(CloogDomain *dom, int level)
{
	int i;
	unsigned dim = cloog_domain_dimension(dom) + dom->nb_par;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	ppl_Pointset_Powerset_C_Polyhedron_t pos, neg;
	ppl_Coefficient_t coeff;
	ppl_Linear_Expression_t le;
	ppl_Constraint_t con;
	mpz_t v;

	mpz_init(v);
	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_space_dimension(&ps, dim, 0) < 0)
		cloog_die("PPL error.\n");

	for (i = 0; i < dom->nb_par; ++i) {
		if (ppl_new_Pointset_Powerset_C_Polyhedron_from_space_dimension(&pos, dim, 0) < 0)
			cloog_die("PPL error.\n");
		if (ppl_new_Linear_Expression(&le) < 0)
			cloog_die("PPL error.\n");
		mpz_set_si(v, 1);
		ppl_assign_Coefficient_from_mpz_t(coeff, v);
		if (ppl_Linear_Expression_add_to_coefficient(le, dim - dom->nb_par + i, coeff) < 0)
			cloog_die("PPL error.\n");
		if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
			cloog_die("PPL error.\n");
		ppl_delete_Linear_Expression(le);
		if (ppl_Pointset_Powerset_C_Polyhedron_add_constraint(pos, con) < 0)
			cloog_die("PPL error.\n");
		ppl_delete_Constraint(con);

		if (ppl_new_Pointset_Powerset_C_Polyhedron_from_space_dimension(&neg, dim, 0) < 0)
			cloog_die("PPL error.\n");
		if (ppl_new_Linear_Expression(&le) < 0)
			cloog_die("PPL error.\n");
		mpz_set_si(v, -1);
		ppl_assign_Coefficient_from_mpz_t(coeff, v);
		if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
			cloog_die("PPL error.\n");
		if (ppl_Linear_Expression_add_to_coefficient(le, dim - dom->nb_par + i, coeff) < 0)
			cloog_die("PPL error.\n");
		if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
			cloog_die("PPL error.\n");
		ppl_delete_Linear_Expression(le);
		if (ppl_Pointset_Powerset_C_Polyhedron_add_constraint(neg, con) < 0)
			cloog_die("PPL error.\n");
		ppl_delete_Constraint(con);

		if (ppl_Pointset_Powerset_C_Polyhedron_upper_bound_assign(pos, neg) < 0)
			cloog_die("PPL error.\n");
		ppl_delete_Pointset_Powerset_C_Polyhedron(neg);

		if (ppl_Pointset_Powerset_C_Polyhedron_intersection_assign(ps, pos) < 0)
			cloog_die("PPL error.\n");
		ppl_delete_Pointset_Powerset_C_Polyhedron(pos);
	}

	ppl_delete_Coefficient(coeff);
	mpz_clear(v);

	return cloog_domain_from_powerset(dom->state, ps, dom->nb_par);
}


/**
 * cloog_scattering_lazy_block function:
 * This function returns 1 if the two scattering functions s1 and s2 given
 * as input are the same (except possibly for the final dimension, where we
 * allow a difference of 1), assuming that the domains on which this
 * scatterings are applied are the same.
 * In fact this function answers the question "can I
 * safely consider the two domains as only one with two statements (a block) ?".
 * - s1 and s2 are the two domains to check for blocking,
 * - scattering is the linked list of all domains,
 * - scattdims is the total number of scattering dimentions.
 *
 * It's actually not that easy to get this right, so we'll just assume
 * that we can't block the two statements.
 */
int cloog_scattering_lazy_block(CloogScattering *s1, CloogScattering *s2,
			    CloogScatteringList *scattering, int scattdims)
{
	return 0;
}


/**
 * cloog_domain_lazy_disjoint function:
 * This function returns 1 if the domains given as input are disjoint, 0 if it
 * is unable to decide.
 */
int cloog_domain_lazy_disjoint(CloogDomain *d1, CloogDomain *d2)
{
	return ppl_Pointset_Powerset_C_Polyhedron_is_disjoint_from_Pointset_Powerset_C_Polyhedron(d1->ps, d2->ps);
} 
 
 
/**
 * cloog_scattering_list_lazy_same function:
 * This function returns 1 if two domains in the list are the same, 0 if it
 * is unable to decide.
 */
int cloog_scattering_list_lazy_same(CloogScatteringList *list)
{
	CloogScatteringList *one, *other;

	for (one = list; one; one = one->next)
		for (other = one->next; other; other = other->next)
			if (ppl_Pointset_Powerset_C_Polyhedron_equals_Pointset_Powerset_C_Polyhedron(one->scatt->dom.ps, other->scatt->dom.ps))
				return 1;
	return 0;
}

int cloog_domain_dimension(CloogDomain *domain)
{
	ppl_dimension_type dim;

	ppl_Pointset_Powerset_C_Polyhedron_space_dimension(domain->ps, &dim);

	return dim - domain->nb_par;
}

int cloog_domain_parameter_dimension(CloogDomain *domain)
{
	return domain->nb_par;
}

int cloog_scattering_dimension(CloogScattering *scatt, CloogDomain *domain)
{
	ppl_dimension_type dim;

	ppl_Pointset_Powerset_C_Polyhedron_space_dimension(scatt->dom.ps, &dim);

	return dim - cloog_domain_dimension(domain) - domain->nb_par;
}

int cloog_domain_isconvex(CloogDomain *domain)
{
	size_t n;

	if (ppl_Pointset_Powerset_C_Polyhedron_omega_reduce(domain->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_size(domain->ps, &n) < 0)
		cloog_die("PPL error.\n");

	return n <= 1;
}


/**
 * cloog_domain_cut_first function:
 * This function splits off and returns the first convex set in the
 * union "domain".  The remainder of the union is returned in rest.
 * The original "domain" itself is destroyed and may not be used
 * after a call to this function.
 */
CloogDomain *cloog_domain_cut_first(CloogDomain *domain, CloogDomain **rest)
{
	ppl_const_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps, ps_first;
	ppl_Pointset_Powerset_C_Polyhedron_iterator_t pit;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_iterator(&pit) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_iterator_dereference(pit, &pol) < 0)
		cloog_die("PPL error.\n");

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps_first, pol) < 0)
		cloog_die("PPL error.\n");

	if (ppl_Pointset_Powerset_C_Polyhedron_drop_disjunct(ps, pit, pit) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Pointset_Powerset_C_Polyhedron_iterator(pit);

	*rest = cloog_domain_from_powerset(domain->state, ps, domain->nb_par);
	cloog_domain_free(domain);
	return cloog_domain_from_powerset((*rest)->state, ps_first, (*rest)->nb_par);
}


/**
 * Given a union domain, try to find a simpler representation
 * using fewer sets in the union.
 * The original "domain" itself is destroyed and may not be used
 * after a call to this function.
 */
CloogDomain *cloog_domain_simplify_union(CloogDomain *domain)
{
	CloogDomain *res;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_pairwise_reduce(ps) < 0)
		cloog_die("PPL error.\n");

	res = cloog_domain_from_powerset(domain->state, ps, domain->nb_par);

	cloog_domain_free(domain);

	return res;
}


/*
 * Check whether the iterator at position pos in polyhedron pol
 * has a fixed value.  If so, return 1 and assign the this value to *value.
 * If update is set, then the fixed value is required to be identical
 * to the original value in *value.
 */
static int polyhedron_dim_is_fixed(ppl_const_Polyhedron_t pol, int pos,
	cloog_int_t *value, int update)
{
	int fixed = 0;
	int row;
	CloogMatrix *M;

	M = cloog_matrix_from_polyhedron(pol);

	for (row = 0; row < M->NbRows; ++row) {
		if (!cloog_int_is_one(M->p[row][1 + pos]) &&
		    !cloog_int_is_neg_one(M->p[row][1 + pos]))
			continue;
		if (!cloog_int_is_zero(M->p[row][0]))
			continue;
		if (cloog_seq_first_non_zero(M->p[row] + 1, pos) != -1)
			continue;
		if (cloog_seq_first_non_zero(M->p[row] + 1 + pos + 1,
					      M->NbColumns - 1 - pos - 2) != -1)
			continue;
		if (cloog_int_is_one(M->p[row][1 + pos]))
			cloog_int_neg(M->p[row][M->NbColumns - 1],
					M->p[row][M->NbColumns - 1]);
		if (update && cloog_int_ne(M->p[row][M->NbColumns - 1], *value))
			break;
		cloog_int_set(*value, M->p[row][M->NbColumns - 1]);
		fixed = 1;
		break;
	}

	cloog_matrix_free(M);

	return fixed;
}


/*
 * Check whether the iterator at position pos in powerset ps
 * has a fixed value.  If so, return 1 and assign the this value to *value.
 */
static int powerset_dim_is_fixed(ppl_Pointset_Powerset_C_Polyhedron_t ps,
	int pos, cloog_int_t *value)
{
	int fixed = 1;
	int first = 1;
	ppl_Pointset_Powerset_C_Polyhedron_const_iterator_t pit, pit_end;
	cloog_int_t v;

	cloog_int_init(v);
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_begin(ps, pit) < 0)
		cloog_die("PPL error.\n");
	if (ppl_new_Pointset_Powerset_C_Polyhedron_const_iterator(&pit_end) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_end(ps, pit_end) < 0)
		cloog_die("PPL error.\n");

	while (!ppl_Pointset_Powerset_C_Polyhedron_const_iterator_equal_test(pit, pit_end)) {
		ppl_const_Polyhedron_t pol;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_dereference(pit, &pol) < 0)
			cloog_die("PPL error.\n");

		if (!polyhedron_dim_is_fixed(pol, pos, &v, !first)) {
			fixed = 0;
			break;
		}

		first = 0;

		if (ppl_Pointset_Powerset_C_Polyhedron_const_iterator_increment(pit) < 0)
			cloog_die("PPL error.\n");
	}

	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit);
	ppl_delete_Pointset_Powerset_C_Polyhedron_const_iterator(pit_end);

	if (fixed && value)
		cloog_int_set(*value, v);

	cloog_int_clear(v);

	return fixed;
}


/**
 * cloog_scattering_lazy_isscalar function:
 * this function returns 1 if the scattering dimension 'dimension' in the
 * scattering 'scatt' is constant.
 * If value is not NULL, then it is set to the constant value of dimension.
 */
int cloog_scattering_lazy_isscalar(CloogScattering *scatt, int dimension,
					cloog_int_t *value)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps = scatt->dom.ps;

	return powerset_dim_is_fixed(ps, dimension, value);
}


/**
 * cloog_domain_lazy_isconstant function:
 * this function returns 1 if the dimension 'dimension' in the
 * domain 'domain' is constant.
 * If value is not NULL, then it is set to the constant value of dimension.
 */
int cloog_domain_lazy_isconstant(CloogDomain *domain, int dimension,
			         cloog_int_t *value)
{
	return powerset_dim_is_fixed(domain->ps, dimension, value);
}


/**
 * cloog_scattering_erase_dimension function:
 * this function returns a CloogDomain structure builds from 'domain' where
 * we removed the dimension 'dimension' and every constraint involving this
 * dimension.
 */
CloogScattering *cloog_scattering_erase_dimension(CloogScattering *scatt,
						int dimension)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	ppl_dimension_type d = dimension;

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, scatt->dom.ps) < 0)
		cloog_die("PPL error.\n");
	if (ppl_Pointset_Powerset_C_Polyhedron_remove_space_dimensions(ps, &d, 1) < 0)
		cloog_die("PPL error.\n");

	return cloog_scattering_from_powerset(scatt->dom.state, ps, scatt->dom.nb_par);
}

/**
 * cloog_domain_cube:
 * Construct and return a dim-dimensional cube, with values ranging
 * between min and max in each dimension.
 */
CloogDomain *cloog_domain_cube(CloogState *state,
				int dim, cloog_int_t min, cloog_int_t max)
{
	int i;
	ppl_Polyhedron_t pol;
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	ppl_Coefficient_t coeff;
	cloog_int_t v;

	if (dim == 0)
		return cloog_domain_universe(state, dim);

	if (ppl_new_C_Polyhedron_from_space_dimension(&pol, dim, 0) < 0)
		cloog_die("PPL error.\n");

	cloog_int_init(v);

	if (ppl_new_Coefficient(&coeff) < 0)
		cloog_die("PPL error.\n");

	for (i = 0; i < dim; ++i) {
		ppl_Linear_Expression_t le;
		ppl_Constraint_t con;

		if (ppl_new_Linear_Expression_with_dimension(&le, dim) < 0)
			cloog_die("PPL error.\n");
		cloog_int_set_si(v, 1);
		ppl_assign_Coefficient_from_mpz_t(coeff, v);
		if (ppl_Linear_Expression_add_to_coefficient(le, i, coeff) < 0)
			cloog_die("PPL error.\n");
		cloog_int_neg(v, min);
		ppl_assign_Coefficient_from_mpz_t(coeff, v);
		if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
			cloog_die("PPL error.\n");

		if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Linear_Expression(le);

		if (ppl_Polyhedron_add_constraint(pol, con) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Constraint(con);
	
		if (ppl_new_Linear_Expression_with_dimension(&le, dim) < 0)
			cloog_die("PPL error.\n");
		cloog_int_set_si(v, -1);
		ppl_assign_Coefficient_from_mpz_t(coeff, v);
		if (ppl_Linear_Expression_add_to_coefficient(le, i, coeff) < 0)
			cloog_die("PPL error.\n");
		ppl_assign_Coefficient_from_mpz_t(coeff, max);
		if (ppl_Linear_Expression_add_to_inhomogeneous(le, coeff) < 0)
			cloog_die("PPL error.\n");

		if (ppl_new_Constraint(&con, le, PPL_CONSTRAINT_TYPE_GREATER_OR_EQUAL) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Linear_Expression(le);

		if (ppl_Polyhedron_add_constraint(pol, con) < 0)
			cloog_die("PPL error.\n");

		ppl_delete_Constraint(con);
	}

	ppl_delete_Coefficient(coeff);
	cloog_int_clear(v);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_C_Polyhedron(&ps, pol) < 0)
		cloog_die("PPL error.\n");

	ppl_delete_Polyhedron(pol);

	return cloog_domain_from_powerset(state, ps, 0);
}


/**
 * cloog_domain_scatter function:
 * This function add the scattering (scheduling) informations to a domain.
 */
CloogDomain *cloog_domain_scatter(CloogDomain *domain, CloogScattering *scatt)
{
	ppl_Pointset_Powerset_C_Polyhedron_t ps;
	CloogDomain *res;
	int nb_scat = cloog_scattering_dimension(scatt, domain);

	if (ppl_new_Pointset_Powerset_C_Polyhedron_from_Pointset_Powerset_C_Polyhedron(&ps, domain->ps) < 0)
		cloog_die("PPL error.\n");
	powerset_insert_dims(ps, 0, nb_scat);
	if (ppl_Pointset_Powerset_C_Polyhedron_intersection_assign(ps, scatt->dom.ps) < 0)
		cloog_die("PPL error.\n");

	res = cloog_domain_from_powerset(domain->state, ps, domain->nb_par);

	cloog_domain_free(domain);

	return res;
}

/* Check if the given list of domains has a common stride on the given level.
 * If so, return a pointer to a CloogStride object.  If not, return NULL.
 *
 * We conservatively return NULL in this backend.
 */
CloogStride *cloog_domain_list_stride(CloogDomainList *list, int level)
{
	return NULL;
}
