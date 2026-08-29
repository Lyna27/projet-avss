#ifndef GETSHARE_H
#define GETSHARE_H

#include <vector>
#include <optional>
#include <cstddef>

#include "polynomial.hpp"
#include "deal.hpp"
#include "proximity_proof.hpp" 

/* 
 * Structure du message diffusé par P_i pour P_j
 * Contient les subshares et leurs preuves.
 */
struct EchoMessage {
    size_t from_index; // index de P_i
    FieldElement a_ij; // a_{i,j} = \phi(x_i, x_j)
    FieldElement b_ij; // b_{i,j} = r(x_i, x_j)
    MerklePath path_r; // Preuve reliant les feuilles b_{i,j} à c_r_Y[i]
    MerklePath path_psi; // Preuve reliant les feuilles psi(x_i, x_j) à c_psi_Y[i]
};

 /*
 * Structure d'une paire (a,b)_{j,i} vérifiée et validée 
 * par P_i après réception du subshare de P_j.
 * P_i récolte ainsi 2 colonnes dans V : 
 * V_phi= {a_{1,i}, a_{2,i} ... a_{n,i}} 
 * V_r= {b_{1,i}, b_{2,i} ...b_{n,i}} 
 */
struct ValidPair {
    size_t from_index; // index de P_j
    FieldElement a_ji; // phi(x_j, x_i)
    FieldElement b_ji; // r(x_j, x_i)
};

/*
 * Structure qui contient le résultat final du protocole GetShare :
 * Pour chauqe P_i ses polynômes univariés locaux de degré t en X  
 * phi(X,x_i) et r(X,x_i) obtenus par interpolation de Lagrange 
 * sur V_phi et V_r.
 * phi(0, x_i), r(0, x_i) sont donc les secrets locaux des P_i 
 */
struct GetShareResult {
    Uni_Poly phi_share; // phi(X,x_i)
    Uni_Poly r_share;  // r(X,x_i)
};

/*
 * Structure de l'etat local complet d'un P_i durant l'AVSS 
*/
struct ParticipantState {
    // Parametres d'identite
    size_t my_index;
    size_t n_Pi;
    size_t threshold;
    std::vector<FieldElement> x_coords; 

    // Donnees issues de Deal Algo
    bool hon_d; 
    FieldElement alpha; 

    // Phase Echo
    std::vector<ValidPair> V; 
    // S'assure d'avoir une partcipation par P_i pour avoir
    // une bonne interpolation sans doublons
    std::vector<bool> seen_echo;

    // Phase Ready
    size_t ready_count;
    bool ready_sent;
    std::vector<bool> rdy_counted;

    // Engagements publics
    std::vector<Commitment> c_r_Y; // Racines des colonnes r(X, x_j)
    std::vector<Commitment> c_psi_Y; // Racines des colonnes \psi(X, x_j)
};

// Fonction d'initialisation de l'etat d'un P_i
ParticipantState init_participant(
    size_t my_index,
    size_t n_Pi,
    size_t threshold,
    const std::vector<FieldElement>&  x_coords
);

/* 
 * Fonction qui fait le traitement du message envoye par le Dealer :
 * Verifie la preuve de proximite et les chemins Merkle personnels.
 * Si valide, met hon_d a true et genere les messages Echo pour les autres.
 */
std::vector<EchoMessage> receive_deal(
    ParticipantState& state,
    const BroadcastMessage& broadcast,
    const PrivateMessage& my_share
);

/*
 * Fonction executee par P_j à la reception d'un Echo de P_i.
 * Verifie les chemins Merkle envoyes par P_j. 
 * Si valide, l'ajoute à l'ensemble V.
 * Retourne true s'il faut déclencher l'envoi 
 * du message Ready des que |V| > 2t.
 */
bool receive_echo(
    ParticipantState& state,
    const EchoMessage& echo
);

/*
 * Fonction qui traite un message Ready.
 * Incrémente le compteur et retourne true si
 * >= t+1 reçus et envoi son propre message Ready.
 */
bool receive_rdy(
    ParticipantState& state,
    size_t from_index
);

/*
 * Fonction de reconstruction.
 * Vérifie si >= 2t+1 Ready et |V| > t.
 * Si oui, effectue l'interpolation de Lagrange et 
 * retourne les polynômes univariés phi(X,x_i) et r(X,x_i).
 */
std::optional<GetShareResult> reconstruct(
    const ParticipantState& state
);

/*
 * Structure qui contient le temps d'exécution et nombre de
 * hachages de l'algorithme GetShare.
 */
struct BenchResult {
    double temps;
    long long nbr_hash;
};

/*
 * Fonction pour exécuter le protocole GetShare complet.
 * avec le temps d'exécution et le nombre de hachages effectues.
 */
BenchResult run_getshare(
    size_t n,
    size_t t,
    const std::vector<FieldElement>& pts,
    const DealResult& res
);

#endif // GETSHARE_H