#ifndef PROXIMITY_PROOF_H
#define PROXIMITY_PROOF_H
#include "polynomial.hpp"
#include "merkle.hpp"
#include <vector>
#include <cstddef>

// Structure le la preuve de proximite
struct ProximityProof {
    FieldElement rho; 
    Uni_Poly h_X; // = psi(X,rho)
    struct  RowOpen {
        // au moins t+1 indices reveles
        std::vector<size_t> indices; 
        // valeurs de psi correspondantes
        std::vector<FieldElement> evals; 
        // leurs chemins
        std::vector<MerklePath> paths;   
    };
    // Ouverture pour chaque P_i
    std::vector< RowOpen> openings; 
};


/*
 * Rend l'IOPP non-interactive avec Fiat-Shamir : 
 * D cree son propre defi rho en hachant ses propres 
 * engagements publics. 
 * S'il modifie son polynome secret pour tricher, le hache
 * change, rho aussi, et la triche est demasquee.
 */
FieldElement fiat_shamir_rho(
    // racine de psi dans MerkleTree
    const Commitment& c_psi,
    // racine de psi(x_i,Y)
    const std::vector<Commitment>& c_psi_Y
);

/* Prouve que deg(psi(X,Y)) <= t en calculant h_X = psi(X,rho)
 * Ouvre au moins t+1 cases du MerkleTree et fournit leurs
 * chemins jusqu'aux racines c_psi_Y pour chaque P_i.
 */
ProximityProof prove_prox(
    const BivPoly& psi,
    size_t t,
    const std::vector<FieldElement>& eval_pts,
    const Commitment& c_psi,
    const std::vector<Commitment>& c_psi_Y,
    const std::vector<MerkleTree::Tree>& row_trees
);

/*
 * Fonction pour la verification de la validite des shares en 
 * verifiant la preuve de proximite reçue.
 * Les P_i : Valident le calcul de Fiat-Shamir, contrôlent le degre de 
 * h_X, et verifie toutes les preuves d'inclusion Merkle.
 */
bool verify_prox(
    const ProximityProof& proof,
    size_t t,
    const std::vector<FieldElement>& eval_pts,
    const Commitment& c_psi,
    const std::vector<Commitment>& c_psi_Y
);

#endif // PROXIMITY_PROOF_H