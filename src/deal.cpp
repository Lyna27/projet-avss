#include "deal.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <stdexcept>

FieldElement compute_alpha(
    const std::vector<Commitment>& c_r_Y,
    const Commitment& c_phi
) {
    // Tag de separation de domaine de fiat_shamir_rho pour eviter les collisions.
    static const uint8_t TAG[] = {
        'A','L','P','H','A'
    };
    std::vector<uint8_t> buf;
    buf.insert(buf.end(), TAG, TAG + sizeof(TAG));
    for (const auto& c : c_r_Y)
        buf.insert(buf.end(), c.begin(), c.end());
    buf.insert(buf.end(), c_phi.begin(), c_phi.end());

    Hash_type h = hash_bytes(buf.data(), buf.size());
    NTL::ZZ z;
    NTL::ZZFromBytes(z, h.data(), static_cast<long>(HASH_SIZE));
    FieldElement alpha;
    NTL::conv(alpha, z % NTL::ZZ_p::modulus());
    return alpha;
}

DealResult run_deal(
    const FieldElement& secret,
    size_t n,
    size_t t,
    const std::vector<FieldElement>& x_coords
) {
    if (x_coords.size() != n)
        throw std::invalid_argument("erreur : x_coords.size() != n");
    if (n < t + 1)
        throw std::invalid_argument("erreur : n < t+1");

    // phi(X,Y) genere aleatoirement avec phi(0,0) = s 
    BivPoly phi = rand_biv_poly(t, secret);

    // r(X,Y) genere aleatoirement 
    BivPoly r = rand_biv_poly(t);
   
    std::vector<std::vector<FieldElement>> matrix_a(n, std::vector<FieldElement>(n));
    std::vector<std::vector<FieldElement>> matrix_b(n, std::vector<FieldElement>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            // matrix_a[i][j] = phi(x_i,x_j) 
            matrix_a[i][j] = phi.evaluate(x_coords[i], x_coords[j]);
            // matrix_b[i][j] = r(x_i,x_j)  
            matrix_b[i][j] = r.evaluate(x_coords[i], x_coords[j]);
        }

    // c_{r(X,x_i)} = MT.Commit({b_{j,i}}_{j}) 
    std::vector<Commitment> c_r_X(n);
    for (size_t i = 0; i < n; ++i) {
        std::vector<FieldElement> col;
        for (size_t j = 0; j < n; ++j) col.push_back(matrix_b[j][i]);
        c_r_X[i] = MerkleTree::commit(col);
    }

    // c_{r(x_i,Y)} = MT.Commit({b_{i,j}}_{j}) 
    std::vector<Commitment> c_r_Y(n);
    for (size_t i = 0; i < n; ++i)
        c_r_Y[i] = MerkleTree::commit(matrix_b[i]);

    // c_{a(x_i,Y)} = MT.Commit({a_{i,j}}_{j})
    std::vector<MerkleTree::Tree> tree_a(n);
    std::vector<Commitment> c_a_X_Y(n);
    for (size_t i = 0; i < n; ++i) {
        tree_a[i]  = MerkleTree::tree_field(matrix_a[i]);
        c_a_X_Y[i] = tree_a[i].root();
    }

    // c_phi = MT.Commit({c_{a(x_i,Y)}}_{i}) 
    MerkleTree::Tree tree_c_phi = MerkleTree::tree_commits(c_a_X_Y);
    Commitment c_phi = tree_c_phi.root();

    // pi_i = chemin de la feuille c_{a(x_i,Y)} vers c_phi 
    std::vector<MerklePath> pi_paths(n);
    for (size_t i = 0; i < n; ++i)
        pi_paths[i] = tree_c_phi.path(i);

    // alpha = H(c_r_Y[0] || ... || c_r_Y[n-1] || c_phi) 
    FieldElement alpha = compute_alpha(c_r_Y, c_phi);
    // psi = phi + alpha * r 
    BivPoly psi = add(phi, scalar_mul(r, alpha));

    //  Matrice d'evaluation de psi 
    std::vector<std::vector<FieldElement>> psi_mat(n, std::vector<FieldElement>(n));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            psi_mat[i][j] = psi.evaluate(x_coords[i], x_coords[j]);

    // c_{psi(X,x_i)} = MT.Commit({psi(x_j,x_i)}_{j}) 
    // Feuille j = psi_mat[j][i], colonne i de la grille psi
    std::vector<Commitment> c_psi_X(n);
    for (size_t i = 0; i < n; ++i) {
        std::vector<FieldElement> col;
        for (size_t j = 0; j < n; ++j) col.push_back(psi_mat[j][i]);
        c_psi_X[i] = MerkleTree::commit(col);
    }

    // c_{psi(x_i,Y)} = MT.Commit({psi(x_i,x_j)}_{j}) 
    // Arbres complets conserves  
    std::vector<MerkleTree::Tree> row_trees_psi(n);
    std::vector<Commitment> c_psi_Y(n);
    for (size_t i = 0; i < n; ++i) {
        row_trees_psi[i] = MerkleTree::tree_field(psi_mat[i]);
        c_psi_Y[i] = row_trees_psi[i].root();
    }

    // c = MT.Commit({c_{psi(X,x_i)}}_{i}) 
    Commitment c = MerkleTree::commit(c_psi_X);

    // pi = ProximityProof(psi, t) 
    ProximityProof pi = prove_prox(psi, t, x_coords, c, c_psi_Y, row_trees_psi);

    // Broadcast c, pi, c_phi, {c_r_X, c_r_Y, c_psi_X, c_psi_Y}_i
    BroadcastMessage bcast;
    bcast.c = c;
    bcast.pi = std::move(pi);
    bcast.c_phi = c_phi;
    bcast.c_r_X = std::move(c_r_X);
    bcast.c_r_Y = std::move(c_r_Y);
    bcast.c_psi_X = std::move(c_psi_X);
    bcast.c_psi_Y = std::move(c_psi_Y);

    // Envoie {a_{i,j}, b_{i,j}}_j, pi_i a chaque P_i 
    std::vector<PrivateMessage> privates(n);
    for (size_t i = 0; i < n; ++i) {
        privates[i].a_i = matrix_a[i];
        privates[i].b_i = matrix_b[i];
        privates[i].pi_i = pi_paths[i];
    }

    DealerInternalState internal;
    internal.c_a_X_Y = std::move(c_a_X_Y);
    internal.tree_c_phi = std::move(tree_c_phi);

    DealResult res;
    res.broadcast_data = std::move(bcast);
    res.private_data = std::move(privates);
    res.internal_state = std::move(internal);
    return res;
}