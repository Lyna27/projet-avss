#ifndef FIELD_H
#define FIELD_H

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

/* 
* Alias pour le corps fini F_p, base sur NTL::ZZ_p qui 
* fournit l'arithmetique modulaire, principalement les
* operateurs +, -, *, ==, != et l'inversion 1/a.
 */
using FieldElement = NTL::ZZ_p;

/*
 * Definit le module p du corps fini F_p.
 * Ce module doit etre choisi suffisamment grand selon le 
 * niveau de securite exige et le besoin d'une FFT pour 
 * les prochaines preuves de proximite.
 * Doit etre appelee une seule fois avant toute operation 
 * algebrique sur FieldElement.
 */
void init_field(const NTL::ZZ& prime);

/* 
 * Fonction qui genere un element aleatoire 
 * uniformement dans F_p, qui doit etre configure au
 * prealable par init_field.
 */
FieldElement random_element();

#endif // FIELD_H