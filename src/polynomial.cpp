#include "polynomial.hpp"
#include <algorithm>

/*
 * Evaluation point a point d'un polynome bivarie P(X,Y).
 * Applique la methode de Horner sur la variable X.
 * Chaque coefficient c_i(Y) est prealablement evalue en y.
 */
FieldElement BivPoly::evaluate(const FieldElement& x, const FieldElement& y) const {
    FieldElement result; 
    for (long i = static_cast<long>(coeffs_X.size()) - 1; i >= 0; --i) {
        FieldElement c_i_at_y;
        NTL::eval(c_i_at_y, coeffs_X[static_cast<size_t>(i)], y);
        result = result * x + c_i_at_y;
    }
    return result;
}

/*
 * Genere aleatoirement le polynome de masquage r(X,Y) de 
 * degre total <= t
 * Chaque coefficient de X^i est un polynome univarie de degre 
 * <= t-i en Y, dont les termes sont tires uniformement 
 * au hasard.
 */
BivPoly rand_biv_poly(size_t t) {
    BivPoly P;
    P.coeffs_X.resize(t + 1);
    for (size_t i = 0; i <= t; ++i) {
        long degree_bound = static_cast<long>(t - i);
        for (long j = 0; j <= degree_bound; ++j) {
            NTL::SetCoeff(P.coeffs_X[i], j, random_element());
        }
    }
    return P;
}

/*
 * Generation du polynome secret phi(X,Y).
 * Tire d'abord un polynome aleatoire, puis determine 
 * son terme constant a_00 comme secret.
 */
BivPoly rand_biv_poly(size_t t, const FieldElement& a_00) {
    BivPoly P = rand_biv_poly(t);
    NTL::SetCoeff(P.coeffs_X[0], 0, a_00);
    return P;
}

/*
 * Additionne deux polynomes bivaries en effectuant
 * la somme terme a terme des sous-polynomes. 
 * Le vecteur resultat est automatiquement ajuste a la dimension 
 * maximale, avec un remplissage par le polynome nul pour combler 
 * les ecarts de degre.
 */
BivPoly add(const BivPoly& a, const BivPoly& b) {
    BivPoly result;
    size_t max_size = std::max(a.coeffs_X.size(), b.coeffs_X.size());
    result.coeffs_X.resize(max_size);
    Uni_Poly zero; 
    for (size_t i = 0; i < max_size; ++i) {
        const Uni_Poly& ai = (i < a.coeffs_X.size()) ? a.coeffs_X[i] : zero;
        const Uni_Poly& bi = (i < b.coeffs_X.size()) ? b.coeffs_X[i] : zero;
        result.coeffs_X[i] = ai + bi;
    }
    return result;
}

/*
 * Multiplie un polynome bivarie par un scalaire en
 * multipliant chaque coefficient de X^i par le scalaire.
 */
BivPoly scalar_mul(const BivPoly& p, const FieldElement& scalar) {
    BivPoly result;
    result.coeffs_X.resize(p.coeffs_X.size());
    for (size_t i = 0; i < p.coeffs_X.size(); ++i) {
        NTL::mul(result.coeffs_X[i], p.coeffs_X[i], scalar);
    }
    return result;
}


/*
 * Calcul de P(x_i, Y) par la methode de Horner.
 */
Uni_Poly eval_X(const BivPoly& P, const FieldElement& x_i) {
    Uni_Poly result; 
    for (long k = static_cast<long>(P.coeffs_X.size()) - 1; k >= 0; --k) {
        Uni_Poly scaled;
        NTL::mul(scaled, result, x_i);
        result = scaled + P.coeffs_X[static_cast<size_t>(k)];
    }
    return result;
}

/*
 * Calcul de P(X, y_i) par evaluation des coefficients.
 */
Uni_Poly eval_Y(const BivPoly& P, const FieldElement& y_i) {
    Uni_Poly result;
    for (size_t k = 0; k < P.coeffs_X.size(); ++k) {
        FieldElement scalar_coeff;
        NTL::eval(scalar_coeff, P.coeffs_X[k], y_i);
        NTL::SetCoeff(result, static_cast<long>(k), scalar_coeff);
    }
    return result;
}