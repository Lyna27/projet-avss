#include "proximity_proof.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>
#include <stdexcept>
#include <cstring>

/*
 * Fonction qui geenere le defi non interactif Fiat-Shamir rho.
 * Hache la concatenation d'un tag de domaine et des engagements publics.
 */
FieldElement fiat_shamir_rho(
    // Racine de l'arbre global psi(X,Y)
    const Commitment& c_psi,
    // Racines des arbres de psi(X, x_i) = share de chaque P_i
    const std::vector<Commitment>& c_psi_Y
) {

    /*
     * Tags de separation de domaine pour la protection contre le rejeu.
     * Evite que Hash(a) donne le meme resultat lors d'un 2eme hachage.
     * On utilise les 13 octets en code ASCII "AVSS_PROX_RHO" 
     */ 
    static const uint8_t TAG[] = {
        'P','R','O','X','_','R','H','O'
    };

 // Preallocation et construction du buf : TAG || c_psi || c_psi_Y
    std::vector<uint8_t> buf;
    buf.reserve(sizeof(TAG) + HASH_SIZE * (1 + c_psi_Y.size()));
    buf.insert(buf.end(), TAG, TAG + sizeof(TAG));
    buf.insert(buf.end(), c_psi.begin(), c_psi.end());
    for (const auto& c : c_psi_Y)
        buf.insert(buf.end(), c.begin(), c.end());

    // Hachage du tableau buf
    Hash_type h = hash_bytes(buf.data(), buf.size());
    NTL::ZZ z;
    NTL::ZZFromBytes(z, h.data(), static_cast<long>(HASH_SIZE));
    FieldElement rho;
    // Reduction modulaire dans F_p
    NTL::conv(rho, z % NTL::ZZ_p::modulus());
    return rho;
}

ProximityProof prove_prox(
    const BivPoly& psi,
    size_t t,
    const std::vector<FieldElement>& eval_pts,
    const Commitment& c_psi,
    const std::vector<Commitment>& c_psi_Y,
    const std::vector<MerkleTree::Tree>& row_trees
) {
    const size_t n = eval_pts.size(); // nbr de P_i

    if (n < t + 1) {
        throw std::invalid_argument(
            "erreur : n < t+1 ; pas assez de points pour interpoler.");
    }
    if (row_trees.size() != n || c_psi_Y.size() != n) {
        throw std::invalid_argument(
            "erreur : tailles incoherentes entre n, row_trees et c_psi_Y).");
    }

    ProximityProof proof;
    // Generation du defi non-interactif via Fiat-Shamir
    proof.rho = fiat_shamir_rho(c_psi, c_psi_Y);
    // h(X) = psi(X, rho)
    proof.h_X = eval_Y(psi, proof.rho);

    // Generation des t+1 ouvertures pour chaque P_i.
    proof.openings.resize(n);
    for (size_t i = 0; i < n; ++i) {
        ProximityProof:: RowOpen& op = proof.openings[i];

        for (size_t j = 0; j <= t; ++j) {
            op.indices.push_back(j);
            // calcul de psi(x_i, x_j) 
            op.evals.push_back(psi.evaluate(eval_pts[i], eval_pts[j]));
            // Extraction du chemin de l'arbre pour chaque feuille j
            op.paths.push_back(row_trees[i].path(j));
        }
    }

    return proof;
}

bool verify_prox(
    const ProximityProof& proof,
    size_t t,
    const std::vector<FieldElement>& eval_pts,
    const Commitment& c_psi,
    const std::vector<Commitment>& c_psi_Y
) {
    const size_t n = eval_pts.size();

    // Verification des dimensions
    if (proof.openings.size() != n || c_psi_Y.size() != n) {
        return false;
    }
 
    // Verification du defi non-interactif de Fiat-Shamir
    FieldElement expected_rho = fiat_shamir_rho(c_psi, c_psi_Y);
    if (proof.rho != expected_rho) {
        return false;
    }

    // Verification que deg(h(X)) <= t
    if (NTL::deg(proof.h_X) > static_cast<long>(t)) {
        return false;
    }

    for (size_t i = 0; i < n; ++i) {
        const ProximityProof::RowOpen& op = proof.openings[i];

        // Verification des t+1 points
        if (op.indices.size() != t + 1 ||
            op.evals.size() != t + 1 ||
            op.paths.size() != t + 1) {
            return false;
        }

        /* 
         * Verification que les feuilles correspondent a 
         * racine avec le chemin dans l'arbre de Merkle 
         */
        for (size_t k = 0; k <= t; ++k) {
            Hash_type leaf = hash_element(op.evals[k]);
            if (!verify_path(c_psi_Y[i], leaf, op.paths[k])) {
                return false;
            }
        }

        // Reconstruction du share univarie p_i par interpolation de Lagrange
        NTL::vec_ZZ_p xa, ya;
        xa.SetLength(static_cast<long>(t + 1));
        ya.SetLength(static_cast<long>(t + 1));
        for (size_t k = 0; k <= t; ++k) {
            xa[static_cast<long>(k)] = eval_pts[op.indices[k]];
            ya[static_cast<long>(k)] = op.evals[k];
        }
        Uni_Poly p_i;
        NTL::interpolate(p_i, xa, ya);

        // Test d'intersection : p_i(rho) == h(x_i)
        FieldElement p_i_rho, h_xi;
        NTL::eval(p_i_rho, p_i, proof.rho);
        NTL::eval(h_xi, proof.h_X, eval_pts[i]);

        if (p_i_rho != h_xi) {
            return false;
        }
    }

    return true;
}