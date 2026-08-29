/*
 * Tests unitaires pour le fichier polynomial.cpp.
 * Validation de l'evaluation et generation aleatoire 
 * des polynomes bivaries.
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "polynomial.hpp"
#include "field.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>

// Fonctions Utilitaires pour les tests

/* Fonction qui convertit un entier v de type long en un
 * element de corps de type FieldElement avec fe(v) en
 * utilisant NTL::conv. 
 */
static FieldElement fe(long v) {
    FieldElement x;
    NTL::conv(x, v);
    return x;
}

/* Construit P(X,Y) = 1 + 2X + 3Y en fixant son degré et 
 * ses coefficients en X.
 */
static BivPoly create_poly() {
    BivPoly P;
    P.coeffs_X.resize(2); // P a 2 coeffs en X, donc de degré 1
    
    // Construction du terme de degre 0 en X : c_0(Y) = 1 + 3Y
    NTL::SetCoeff(P.coeffs_X[0], 0, fe(1)); // c_0,0(Y)= 1
    NTL::SetCoeff(P.coeffs_X[0], 1, fe(3)); // c_0,1(Y)= 3

    // Construction du terme de degre 1 en X : c_1(Y) = 2
    NTL::SetCoeff(P.coeffs_X[1], 0, fe(2)); // c_1,0(Y)= 2
    return P;
}

// Tests Unitaires sur la fonction evaluate()

TEST_CASE("Le polynome vide retourne toujours 0") {
    init_field(NTL::ZZ(101)); // Z_101
    BivPoly zero_poly; // P(X,Y) = 0
    CHECK(zero_poly.evaluate(fe(0), fe(0)) == fe(0));
    CHECK(zero_poly.evaluate(fe(5), fe(7)) == fe(0));
}

TEST_CASE("Le polynome constant retourne toujours sa valeur") {
    init_field(NTL::ZZ(101));
    BivPoly P;
    P.coeffs_X.resize(1); // deg(P) = 0 en X
    NTL::SetCoeff(P.coeffs_X[0], 0, fe(42)); // P(X,Y) = 42
    CHECK(P.evaluate(fe(0), fe(0)) == fe(42));
    CHECK(P.evaluate(fe(5), fe(7)) == fe(42));
    CHECK(P.evaluate(fe(99), fe(99)) == fe(42));
}

// Test d'evaluation sans reduction modulaire 
TEST_CASE("L'evaluation du polynome P(X,Y) = 1 + 2X + 3Y retourne les bonnes valeurs") {
    init_field(NTL::ZZ(101));
    BivPoly P = create_poly(); // P(X,Y) = 1 + 2X + 3Y (mod 101)
    // Tests d'evaluation
    CHECK(P.evaluate(fe(0), fe(0)) == fe(1));   // 0 + 1 = 1
    CHECK(P.evaluate(fe(1), fe(1)) == fe(6));   // 1 + 2 + 3 = 6
    CHECK(P.evaluate(fe(2), fe(3)) == fe(14));  // 1 + 4 + 9 = 14
    CHECK(P.evaluate(fe(0), fe(5)) == fe(16));  // 1 + 15 = 16
}

// Test d'evaluation avec reduction modulaire
TEST_CASE("La reduction modulaire est correcte") {
    init_field(NTL::ZZ(7)); // Z_7
    BivPoly P = create_poly(); // P = 1 + 2X + 3Y (mod 7)
    CHECK(P.evaluate(fe(3), fe(4)) == fe(5)); // 19 ≡ 5 (mod 7)
    CHECK(P.evaluate(fe(6), fe(6)) == fe(3)); // 31 ≡ 3 (mod 7)
}

// Tests Unitaires sur la fonction rand_biv_poly() 

TEST_CASE(" Le coefficient constant est bien c_0(Y) = phi(0,0) = secret ") {
    init_field(NTL::ZZ(101)); // Z_101
    for (size_t t = 0; t <= 4; ++t) { // t = deg(P)
        BivPoly phi = rand_biv_poly(t, fe(42)); // phi(0,0) = 42 = secret
        CHECK(phi.evaluate(fe(0), fe(0)) == fe(42));
    }
}

TEST_CASE("Le nombre de coefficients est bien t+1 pour tout degré t") {
    init_field(NTL::ZZ(101)); 
    for (size_t t = 0; t <= 5; ++t) {
        BivPoly P = rand_biv_poly(t, fe(1)); // phi(0,0) = 1
        CHECK(P.coeffs_X.size() == t + 1);
    }
}

TEST_CASE("Le degre de chaque terme X^iY^j est respecte : i+j <= t donc j <= t-i") {
    init_field(NTL::ZZ(101));
    for (size_t t = 1; t <= 4; ++t) {
        BivPoly P = rand_biv_poly(t, fe(1));
        for (size_t i = 0; i <= t; ++i) {
            long d = NTL::deg(P.coeffs_X[i]); // deg(c_i(Y))
            // NTL::deg retourne -1 pour le polynome nul 
            CHECK(d <= static_cast<long>(t - i));
        }
    }
}

TEST_CASE(" Le polynome constant est toujours obtenu pour t=0 ") {
    init_field(NTL::ZZ(101));
    FieldElement s = fe(77); // secret
    BivPoly P = rand_biv_poly(size_t(0), s); // P(X,Y) = s
    CHECK(P.coeffs_X.size() == size_t(1));
    // phi(0,0) = s
    CHECK(P.evaluate(fe(0), fe(0)) == s);
    // Renvoie s en tout point
    CHECK(P.evaluate(fe(5), fe(5)) == s);
    CHECK(P.evaluate(fe(99), fe(0)) == s);
}

// Tests unitaires sur la fonction add() 

TEST_CASE("L'evaluation de la somme est egale a la somme des evaluations") {
    init_field(NTL::ZZ(101));
    BivPoly a = rand_biv_poly(size_t(3), fe(10));
    BivPoly b = rand_biv_poly(size_t(3), fe(20));
    BivPoly s = add(a, b); // s = a + b
    // Verification que (a+b)(x,y) = a(x,y) + b(x,y) = s(x,y)
    for (long xv = 0; xv < 6; ++xv)
        for (long yv = 0; yv < 6; ++yv)
            CHECK(s.evaluate(fe(xv), fe(yv)) ==
                  a.evaluate(fe(xv), fe(yv)) + b.evaluate(fe(xv), fe(yv)));
}

TEST_CASE("L'addition est bien commutative : a + b = b + a en tout point") {
    init_field(NTL::ZZ(101));
    BivPoly a = rand_biv_poly(size_t(2), fe(5));
    BivPoly b = rand_biv_poly(size_t(2), fe(9));
    BivPoly ab = add(a, b); // a + b
    BivPoly ba = add(b, a); // b + a
    for (long xv = 0; xv < 5; ++xv)
        for (long yv = 0; yv < 5; ++yv)
            CHECK(ab.evaluate(fe(xv), fe(yv)) == ba.evaluate(fe(xv), fe(yv)));
}

TEST_CASE("Le polynome nul est bien l'element neutre de l'addition : a + 0 = a") {
    init_field(NTL::ZZ(101));
    BivPoly a = rand_biv_poly(size_t(2), fe(7));
    BivPoly zero_poly; 
    BivPoly r = add(a, zero_poly); // r = a + 0
    for (long xv = 0; xv < 5; ++xv)
        for (long yv = 0; yv < 5; ++yv)
            CHECK(r.evaluate(fe(xv), fe(yv)) == a.evaluate(fe(xv), fe(yv)));
}

TEST_CASE("La taille de l'addition est correctement ajustee malgre des degres differents") {
    init_field(NTL::ZZ(101));
    BivPoly a = rand_biv_poly(size_t(1), fe(3)); // deg(a) = 1 dont 2 coeffs en X
    BivPoly b = rand_biv_poly(size_t(3), fe(8)); // deg(b) = 3 dont 4 coeffs en X
    BivPoly s = add(a, b);
    // Verification que la taille de l'addition est max(2, 4) = 4
    CHECK(s.coeffs_X.size() == size_t(4)); 
    for (long xv = 0; xv < 5; ++xv)
        for (long yv = 0; yv < 5; ++yv)
            CHECK(s.evaluate(fe(xv), fe(yv)) ==
                  a.evaluate(fe(xv), fe(yv)) + b.evaluate(fe(xv), fe(yv)));
}

// Tests unitaires sur la fonction scalar_mul() 

TEST_CASE("L'equation (alpha * P)(x,y) = alpha * P(x,y) est bien correcte") {
    init_field(NTL::ZZ(101));
    BivPoly P = rand_biv_poly(size_t(3), fe(1));
    FieldElement alpha = fe(13);
    BivPoly sc = scalar_mul(P, alpha); // sc = 13 * P
    for (long xv = 0; xv < 6; ++xv)
        for (long yv = 0; yv < 6; ++yv)
            CHECK(sc.evaluate(fe(xv), fe(yv)) ==
                  alpha * P.evaluate(fe(xv), fe(yv)));
}

TEST_CASE("L'equation 0 * P(x,y) = 0 pour tout (x,y) est bien correcte") {
    init_field(NTL::ZZ(101));
    BivPoly P = rand_biv_poly(size_t(3), fe(42));
    BivPoly r = scalar_mul(P, fe(0));
    // Verification de l'equation 
    for (long xv = 0; xv < 6; ++xv)
        for (long yv = 0; yv < 6; ++yv)
            CHECK(r.evaluate(fe(xv), fe(yv)) == fe(0)); 
}

TEST_CASE("1 est bien l'element neutre de la multiplication : 1 * P = P") {
    init_field(NTL::ZZ(101));
    BivPoly P = rand_biv_poly(size_t(2), fe(33));
    BivPoly r = scalar_mul(P, fe(1)); // r = 1 * P
    for (long xv = 0; xv < 5; ++xv)
        for (long yv = 0; yv < 5; ++yv)
            CHECK(r.evaluate(fe(xv), fe(yv)) == P.evaluate(fe(xv), fe(yv)));
}

// Tests unitaires sur la fonction eval_X() 

TEST_CASE("eval_X(x,Y) pour Y = y est identique a evaluate(x,y)") {
    init_field(NTL::ZZ(101));
    BivPoly P = rand_biv_poly(size_t(3), fe(7));
    for (long xv = 0; xv < 6; ++xv) {
        FieldElement x = fe(xv);
        // X est fixe a x et le polynome univarie rx = P(x, Y) en Y est cree 
        Uni_Poly rx = eval_X(P, x);
        for (long yv = 0; yv < 6; ++yv) {
            FieldElement y = fe(yv), via_r, direct; 
            NTL::eval(via_r, rx, y); // Evaluation de rx en y
            direct = P.evaluate(x, y); // Evaluation bivariee de P en (x,y)
            CHECK(via_r == direct);
        }
    }
}

TEST_CASE("P_0(Y) = c_0(Y) = coeffs_X[0] est correcte") {
    init_field(NTL::ZZ(101));
    BivPoly P = rand_biv_poly(size_t(3), fe(5));
    Uni_Poly r0 = eval_X(P, fe(0)); // r0 = P(0,Y)
    for (long yv = 0; yv < 10; ++yv) {
        FieldElement y = fe(yv), via_r, via_c;
        NTL::eval(via_r, r0, y); // Evaluation de P(0,Y) en Y = y
        NTL::eval(via_c, P.coeffs_X[0], y); // Evaluation du coefficient c_0(Y) en Y = y
        CHECK(via_r == via_c);
    }
}

// Tests unitaires sur la fonction eval_Y() 

TEST_CASE("eval_Y(X,y) pour X = x est identique a evaluate(x,y)") {
    init_field(NTL::ZZ(101));
    BivPoly P = rand_biv_poly(size_t(3), fe(9));
    for (long yv = 0; yv < 6; ++yv) {
        FieldElement y = fe(yv);
        // Y est fixe a y et le polynome univarie ry = P(X, y) en X est cree 
        Uni_Poly ry = eval_Y(P, y);
        for (long xv = 0; xv < 6; ++xv) {
            FieldElement x = fe(xv), via_r, direct;
            NTL::eval(via_r, ry, x); // Evaluation de ry en x
            direct = P.evaluate(x, y); // Evaluation bivariee de P en (x,y)
            CHECK(via_r == direct);
        }
    }
}

TEST_CASE("L'equation psi = phi + alpha * r est correcte") {
    init_field(NTL::ZZ(101));
    FieldElement secret = fe(42);
    size_t t = 3;
    // Polynome secret phi(X,Y) avec phi(0,0) = secret
    BivPoly phi = rand_biv_poly(t, secret);

    // Polynome de masquage r(X,Y)
    BivPoly r = rand_biv_poly(t);
    // alpha est un element aleatoire du corps
    FieldElement alpha = random_field_element();
    // Polynome psi(X,Y) = phi(X,Y) + alpha * r(X,Y) utilisé par le Dealer
    BivPoly psi = add(phi, scalar_mul(r, alpha));

    // Verification que phi(0,0) est toujours egale a la valeur du secret
    CHECK(phi.evaluate(fe(0), fe(0)) == secret);

    // Verification que psi(x,y) = phi(x,y) + alpha * r(x,y) en tout point
    for (long xv = 0; xv < 6; ++xv)
        for (long yv = 0; yv < 6; ++yv) {
            FieldElement x = fe(xv), y = fe(yv);
            CHECK(psi.evaluate(x, y) ==
                  phi.evaluate(x, y) + alpha * r.evaluate(x, y));
        }

    // Verification que psi(0,0) = secret + alpha*r(0,0) pas forcement = secret
    CHECK(psi.evaluate(fe(0), fe(0)) ==
          secret + alpha * r.evaluate(fe(0), fe(0)));
}