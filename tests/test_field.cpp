/*
 * Tests unitaires pour le fichier field.cpp.
 * Validation de la configuration du corps fini F_p et 
 * de la robustesse des tirages aleatoires.
 * L'infrastructure de test utilisee est doctest.
 */
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "field.hpp"
#include <NTL/ZZ.h>


/*
 * Premier TEST_CASE a lancer qui verifie si les operations  
 * sont bloquees si F_p n'est pas encore configure.
 * Appelle random_element() avant init_field() et s'attend a une erreur.
 */
TEST_CASE("random_element() ne peut fonctionner avant l'appel de init_field()") {
    CHECK_THROWS_AS(random_element(), std::runtime_error);
}


/*
 * Verifie que init_field() fixe correctement le module p.
 * Choisit un nombre premier 101 en parametre de la fonction init_field(), 
 * puis verifie directement avec NTL::ZZ_p::modulus() pour s'assurer que
 * p est bien egal a 101.
 */
TEST_CASE("le module p est bien fixé par init_field()") {
    NTL::ZZ p(101);
    init_field(p);
    CHECK(NTL::ZZ_p::modulus() == p);
}

/*
 * Verifie que init_field() peut etre rappele avec un module
 * different apres une premiere utilisation. 
 * Pour changer la taille de F_p, Il doit etre possible d'ecraser
 * l'ancien module p = 101 par un nouveau p' = 7919.
 * */
TEST_CASE("init_field() peut etre rappelee avec un module different") {
    init_field(NTL::ZZ(101));
    CHECK(NTL::ZZ_p::modulus() == NTL::ZZ(101));

    init_field(NTL::ZZ(7919)); // autre premier
    CHECK(NTL::ZZ_p::modulus() == NTL::ZZ(7919));
}

/*
 * S'assure que x = random_element() est toujours compris entre 0 et p.
 * Effectue 200 tirages aleatoires et verifie que chaque element genere
 * est bien dans l'intervalle [0, p[ en utilisant NTL::rep(x) pour obtenir
 * le representant entier canonique de x.
 */
TEST_CASE("random_element() est toujours une valeur de F_p") {
    NTL::ZZ p(101);
    init_field(p);

    for (int i = 0; i < 200; ++i) {
        FieldElement x = random_element();
        NTL::ZZ representative = NTL::rep(x); 
        CHECK(representative >= NTL::ZZ(0));
        CHECK(representative < p);
    }
}

/* 
 * Verifie que la sequence d'elements produite par random_element() 
 * ne retourne pas toujours la meme valeur.
 * Applique 50 tirages aleatoires et s'assure qu'au moins un element est  
 * different du premier tirage.
 * Avec p = 101, la probabilite que 50 tirages uniformes tombent tous
 * sur la meme valeur que le premier est negligeable : (1/101)^50
*/
TEST_CASE("random_element() ne converge pas vers une valeur constante") {
    init_field(NTL::ZZ(101));

    FieldElement first = random_element();
    bool found_different = false;
    for (int i = 0; i < 50 && !found_different; ++i) {
        if (random_element() != first) {
            found_different = true;
        }
    }
    CHECK(found_different);
}

/*
 * Verifie que random_element() fonctionne pour le corps 
 * fini minimal F_2 = {0, 1}, et qu'il n'y a pas de 
 * concentration sur une seule des deux valeurs 
 * pour 50 tirages aleatoires. 
 */
TEST_CASE("random_element() fonctionne correctement avec p = 2") {
    init_field(NTL::ZZ(2));

    bool found_zero = false;
    bool found_one = false;

    for (int i = 0; i < 50; ++i) {
        FieldElement x = random_element();
        NTL::ZZ rep = NTL::rep(x);
        
        CHECK(rep >= NTL::ZZ(0));
        CHECK(rep < NTL::ZZ(2));

        if (rep == NTL::ZZ(0)) found_zero = true;
        if (rep == NTL::ZZ(1)) found_one = true;
    }

    CHECK(found_zero);
    CHECK(found_one);
}