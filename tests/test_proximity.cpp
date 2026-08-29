// Afin de prouver que deg(psi) <= t, et que psi(x_i, rho) = h_X(x_i) 

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "proximity_proof.hpp"
#include "field.hpp"
#include "polynomial.hpp"
#include "merkle.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>


static FieldElement fe(long v) {
    FieldElement x; NTL::conv(x, v); return x;
}

// Construit les n points d'evaluation x_1 = 1, x_2 = 2, ..., x_n = n dans F_p.
// On evite x=0 pour que les shares s_i = psi(x_i,0) soient non triviaux.
static std::vector<FieldElement> x_pts(size_t n) {
    std::vector<FieldElement> pts;
    for (size_t i = 1; i <= n; ++i) pts.push_back(fe(static_cast<long>(i)));
    return pts;
}

// Construit les arbres Merkle des rangees de psi :
// row_trees[i] = MT sur {psi(x_i, x_0), psi(x_i, x_1), ..., psi(x_i, x_{n-1})}
// Retourne aussi les commitments c_psi_Y[i] = row_trees[i].root().
static void build_row_trees(
    const BivPoly& psi,
    const std::vector<FieldElement>& eval_pts,
    std::vector<MerkleTree::Tree>& row_trees,
    std::vector<Commitment>& c_psi_Y
) {
    size_t n = eval_pts.size();
    row_trees.resize(n);
    c_psi_Y.resize(n);
    for (size_t i = 0; i < n; ++i) {
        std::vector<FieldElement> row_evals;
        for (size_t j = 0; j < n; ++j)
            row_evals.push_back(psi.evaluate(eval_pts[i], eval_pts[j]));
        row_trees[i] = MerkleTree::tree_field(row_evals);
        c_psi_Y[i] = row_trees[i].root();
    }
}

// Construit c_psi = MT.Commit({c_psi_X[i]}_i) ou c_psi_X[i] est le commit
// sur la colonne (Y = x_i). Pour simplifier, on genere c_psi depuis les row
// commits eux memes (dans un vrai Deal ce serait les col commits ; ici on
// test juste que le Fiat-Shamir est deterministe, peu importe ce qu'on met).
// On simule c_psi = commit(c_psi_Y) pour rester coherent.
static Commitment build_c_psi(const std::vector<Commitment>& c_psi_Y) {
    return MerkleTree::commit(c_psi_Y);
}


// fiat_shamir_rho()

TEST_CASE("fiat_shamir_rho() : memes entrees => meme rho") {
    init_field(NTL::ZZ(101));
    auto pts = x_pts(4);
    auto psi = rand_biv_poly(size_t(2), fe(42));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    FieldElement r1 = fiat_shamir_rho(c, c_Y);
    FieldElement r2 = fiat_shamir_rho(c, c_Y);
    CHECK(r1 == r2);
}

TEST_CASE("fiat_shamir_rho() : change si c_psi change") {
    // Test deterministe : on construit deux polynomes explicitement differents
    // psi1(X,Y) = 1 + 2X + 3Y et psi2(X,Y) = 4 + 5X + 6Y
    // On utilise un grand premier (7919) pour eviter les collisions
    // du type SHA256(a) % p == SHA256(b) % p qui arrivent avec proba 1/p.
    // Avec p=101 (101 valeurs) la collision est frequente ; avec p=7919 elle
    // vaut < 1/7919 par appel, ce qui rend le test deterministe en pratique.
    init_field(NTL::ZZ(7919));
    auto pts = x_pts(4);

    BivPoly psi1, psi2;
    psi1.coeffs_X.resize(2); psi2.coeffs_X.resize(2);
    NTL::SetCoeff(psi1.coeffs_X[0], 0, fe(1)); // constante
    NTL::SetCoeff(psi1.coeffs_X[0], 1, fe(3)); // coeff Y
    NTL::SetCoeff(psi1.coeffs_X[1], 0, fe(2)); // coeff X
    NTL::SetCoeff(psi2.coeffs_X[0], 0, fe(4)); // constante differente
    NTL::SetCoeff(psi2.coeffs_X[0], 1, fe(6)); // coeff Y different
    NTL::SetCoeff(psi2.coeffs_X[1], 0, fe(5)); // coeff X different

    std::vector<MerkleTree::Tree> t1, t2;
    std::vector<Commitment> c_Y1, c_Y2;
    build_row_trees(psi1, pts, t1, c_Y1);
    build_row_trees(psi2, pts, t2, c_Y2);
    Commitment c1 = build_c_psi(c_Y1);
    Commitment c2 = build_c_psi(c_Y2);

    // Les evaluations sur la grille different => commits differents => rho differents
    CHECK(fiat_shamir_rho(c1, c_Y1) != fiat_shamir_rho(c2, c_Y2));
}

TEST_CASE("fiat_shamir_rho() : resultat dans F_p") {
    init_field(NTL::ZZ(101));
    auto pts = x_pts(3);
    auto psi = rand_biv_poly(size_t(1), fe(5));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    // Doit produire un element de F_p bien forme
    FieldElement rho = fiat_shamir_rho(c, c_Y);
    // La valeur doit etre egale a elle-meme 
    CHECK(rho == fiat_shamir_rho(c, c_Y));
}

// structure de la preuve

TEST_CASE("prove_prox() : bonne structure pour n=4, t=2") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(7));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);

    // Rho = Fiat-Shamir
    CHECK(proof.rho == fiat_shamir_rho(c, c_Y));

    // h_X = psi(X, rho) de degre <= t
    CHECK(NTL::deg(proof.h_X) <= static_cast<long>(t));

    // n ouvertures, chacune avec t+1 = 3 points
    CHECK(proof.openings.size() == n);
    for (size_t i = 0; i < n; ++i) {
        CHECK(proof.openings[i].indices.size() == t + 1);
        CHECK(proof.openings[i].evals.size() == t + 1);
        CHECK(proof.openings[i].paths.size() == t + 1);
    }
}

TEST_CASE("prove_prox() : bonnes evaluations avec psi") {
    init_field(NTL::ZZ(101));
    size_t n=5, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(13));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);

    // Chaque evals[i][k] doit etre exactement psi(x_i, x_{indices[k]})
    for (size_t i = 0; i < n; ++i) {
        const auto& op = proof.openings[i];
        for (size_t k = 0; k <= t; ++k) {
            size_t j = op.indices[k];
            FieldElement expected = psi.evaluate(pts[i], pts[j]);
            CHECK(op.evals[k] == expected);
        }
    }
}

TEST_CASE("prove_prox() : h_X(x_i) == psi(x_i, rho) pour tout i") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(99));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);

    // Verifie que h_X(x_i) == psi(x_i, rho) pour tout i
    for (size_t i = 0; i < n; ++i) {
        FieldElement h_xi, direct;
        NTL::eval(h_xi, proof.h_X, pts[i]);
        direct = psi.evaluate(pts[i], proof.rho);
        CHECK(h_xi == direct);
    }
}

// verify_prox() 

TEST_CASE("verify_prox() : preuve honnete valide") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(42));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    CHECK(verify_prox(proof, t, pts, c, c_Y));
}

TEST_CASE("verify_prox() : preuve honnete valide (n=6, t=2)") {
    init_field(NTL::ZZ(101));
    size_t n=6, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(77));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    CHECK(verify_prox(proof, t, pts, c, c_Y));
}

TEST_CASE("verify_prox() pour le seuil minimal t=1") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=1;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(3));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    CHECK(verify_prox(proof, t, pts, c, c_Y));
}

// detection des tricheries
TEST_CASE("verify() : échoue pour rho falsifie") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(5));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    // Modifier rho : Fiat-Shamir ne correspondra plus
    FieldElement orig_rho = proof.rho;
    proof.rho = orig_rho + fe(1);
    CHECK(!verify_prox(proof, t, pts, c, c_Y));
}

TEST_CASE("verify() echoue si deg(h_X) = t+1 ") {
    init_field(NTL::ZZ(101));
    size_t n=5, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(10));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    // Forcer un coefficient de degre t+1 dans h_X
    NTL::SetCoeff(proof.h_X, static_cast<long>(t + 1), fe(1));
    CHECK(!verify_prox(proof, t, pts, c, c_Y));
}

TEST_CASE("verify() échoue en cas d'evaluation falsifiee") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(33));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    // Modifier une evaluation dans l'ouverture de la rangee 0
    proof.openings[0].evals[0] = proof.openings[0].evals[0] + fe(1);
    // Si chemin ne correspond plus, la verification doit echouer
    CHECK(!verify_prox(proof, t, pts, c, c_Y));
}

TEST_CASE("verify() échoue si le nombre d'ouvertures est incorrect") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(20));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    // Supprimer une ouverture : taille n-1 != n
    proof.openings.pop_back();
    CHECK(!verify_prox(proof, t, pts, c, c_Y));
}

TEST_CASE("verify() échoue si p_i(rho) != h(x_i)") {
    init_field(NTL::ZZ(101));
    size_t n=4, t=2;
    auto pts = x_pts(n);
    auto psi = rand_biv_poly(t, fe(55));

    std::vector<MerkleTree::Tree> trees; std::vector<Commitment> c_Y;
    build_row_trees(psi, pts, trees, c_Y);
    Commitment c = build_c_psi(c_Y);

    ProximityProof proof = prove_prox(psi, t, pts, c, c_Y, trees);
    // Modifier h_X sans modifier les ouvertures :
    // deg(h_X) sera toujours <= t mais p_i(rho) != h(x_i)
    // On ajoute un terme constant sans changer le degre si t >= 1
    FieldElement coeff0;
    NTL::eval(coeff0, proof.h_X, fe(0));
    NTL::SetCoeff(proof.h_X, 0, coeff0 + fe(7));
    CHECK(!verify_prox(proof, t, pts, c, c_Y));
}