#ifndef DEAL_H
#define DEAL_H

#include <vector>
#include <cstddef>
#include "field.hpp"
#include "polynomial.hpp"
#include "merkle.hpp"
#include "proximity_proof.hpp"  

/* 
 * Structure regroupant les donnees publiques diffusees 
 * dans l'ensemble du reseau.
 */
struct BroadcastMessage {
    Commitment c; // racine de l'arbre global lie a psi
    ProximityProof pi; // preuve NI-IOPP d'evaluation
    Commitment c_phi; // racine de l'arbre global lie a phi
    std::vector<Commitment> c_r_X; // engagements sur les colonnes de r
    std::vector<Commitment> c_r_Y; // engagements sur les lignes de r
    std::vector<Commitment> c_psi_X; // engagements sur les colonnes : psi(X, x_i)
    std::vector<Commitment> c_psi_Y; // engagements sur les shares : psi(x_i,Y)
};

// Structure contenant les parts privees destinees aux P_i.
struct PrivateMessage {
    std::vector<FieldElement> a_i; // a_{i,j} = phi(x_i, x_j)_j  
    std::vector<FieldElement> b_i; // b_{i,j} = r(x_i, x_j)_j 
    MerklePath pi_i;  // chemin MT de la feuille c_a(x_i,Y) vers c_phi 
};

// Etat interne du Dealer
struct DealerInternalState {
    std::vector<Commitment> c_a_X_Y; // feuilles : c_a(x_i,Y)_i
    MerkleTree::Tree tree_c_phi; // arbre complet 
};

// structure des donnees generees par le Dealer
struct DealResult {
    BroadcastMessage broadcast_data; // message public
    std::vector<PrivateMessage> private_data; // tableau des shares prives de taille n
    DealerInternalState internal_state; // tests locaux du Dealer 
};

// Fonction principale du Dealer
DealResult run_deal(
    const FieldElement& secret,
    size_t n,
    size_t t,
    const std::vector<FieldElement>& x_coords
);

/* 
 * Fonction qui genere le defi non-interactif via Fiat-Shamir.
 * Combine le masque r et l'engagement sur phi pour obtenir :
 * alpha = H(c_r(x1,Y) || ... || c_r(xn,Y) || c_phi)
 */
FieldElement compute_alpha(
    const std::vector<Commitment>& c_r_Y,
    const Commitment& c_phi
);

#endif // DEAL_H