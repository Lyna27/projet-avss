#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cstddef>
#include <NTL/ZZ_pX.h>
#include "field.hpp"

/*
 * Alias pour l'anneau des polynomes univaries F_p[X].
 */
using Uni_Poly = NTL::ZZ_pX;

/*
 * Structure de polynomes bivaries P(X,Y) de degre total <= t.
 * Representation asymetrique par vecteur de polynomes univaries :
 * P(X,Y) = somme_{i=0}^{t} X^i * c_i(Y), avec deg(c_i) <= t - i.
 * Vu comme un Poly univarie en X dont les t+1 coeffs c_i(Y) sont
 * des Polys univaries en Y.
 */
struct BivPoly {
    std::vector<Uni_Poly> coeffs_X;

    FieldElement evaluate(const FieldElement& x, const FieldElement& y) const;
};

/*
 * Tire au hasard le polynome secret phi(X,Y). 
 * Determine le terme constant a_00 comme etant le secret secret s.
 */
BivPoly rand_biv_poly(size_t t, const FieldElement& a_00);

/*
 * Tire le polynome de masquage r(X,Y) au hasard qui permettra 
 * de cacher le polynome secret phi(X,Y).
 */
BivPoly rand_biv_poly(size_t t); 

// Addition et multiplication par un scalaire 
BivPoly add(const BivPoly& a, const BivPoly& b);
BivPoly scalar_mul(const BivPoly& p, const FieldElement& scalar);

/*
 * Eevaluations partielles:P(x_i, Y) et P(X, y_i)
 */
Uni_Poly eval_X(const BivPoly& P, const FieldElement& x_i);
Uni_Poly eval_Y(const BivPoly& P, const FieldElement& y_i);

#endif // POLYNOMIAL_H