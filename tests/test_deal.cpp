#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "deal.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>

static FieldElement fe(long v) { FieldElement x; NTL::conv(x,v); return x; }

static std::vector<FieldElement> make_pts(size_t n) {
    std::vector<FieldElement> pts;
    for (size_t i = 1; i <= n; ++i) pts.push_back(fe(static_cast<long>(i)));
    return pts;
}

static FieldElement lagrange_eval(
    const std::vector<FieldElement>& xs,
    const std::vector<FieldElement>& ys,
    const FieldElement& x_star
) {
    NTL::vec_ZZ_p xa, ya;
    xa.SetLength(static_cast<long>(xs.size()));
    ya.SetLength(static_cast<long>(ys.size()));
    for (size_t k = 0; k < xs.size(); ++k) {
        xa[static_cast<long>(k)] = xs[k];
        ya[static_cast<long>(k)] = ys[k];
    }
    Uni_Poly f;
    NTL::interpolate(f, xa, ya);
    FieldElement r; NTL::eval(r, f, x_star); return r;
}

// Structure
TEST_CASE("run_deal() - execution sans exception (n=4, t=1)") {
    init_field(NTL::ZZ(7919));
    CHECK_NOTHROW(run_deal(fe(42), 4, 1, make_pts(4)));
}

TEST_CASE("run_deal() - tailles BroadcastMessage coherentes") {
    init_field(NTL::ZZ(7919));
    size_t n=5, t=2;
    auto res = run_deal(fe(7), n, t, make_pts(n));
    const auto& b = res.broadcast_data;
    CHECK(b.c_r_X  .size() == n);
    CHECK(b.c_r_Y  .size() == n);
    CHECK(b.c_psi_X.size() == n);
    CHECK(b.c_psi_Y.size() == n);
    CHECK(b.c    .size() == HASH_SIZE);
    CHECK(b.c_phi.size() == HASH_SIZE);
}

TEST_CASE("run_deal() - tailles PrivateMessages coherentes") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto res = run_deal(fe(99), n, t, make_pts(n));
    CHECK(res.private_data.size() == n);
    for (size_t i=0;i<n;++i) {
        CHECK(res.private_data[i].a_i.size() == n);
        CHECK(res.private_data[i].b_i.size() == n);
    }
}

// ProximityProof 
TEST_CASE("run_deal() - ProximityProof verifie (n=4, t=1)") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto pts = make_pts(n);
    auto res = run_deal(fe(42), n, t, pts);
    const auto& b = res.broadcast_data;
    CHECK(verify_prox(b.pi, t, pts, b.c, b.c_psi_Y));
}

TEST_CASE("run_deal() - ProximityProof verifie (n=6, t=2)") {
    init_field(NTL::ZZ(7919));
    size_t n=6, t=2;
    auto pts = make_pts(n);
    auto res = run_deal(fe(13), n, t, pts);
    const auto& b = res.broadcast_data;
    CHECK(verify_prox(b.pi, t, pts, b.c, b.c_psi_Y));
}

// Chemins Merkle pi_i 
TEST_CASE("run_deal() - chemins pi_i valides contre c_phi") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto pts = make_pts(n);
    auto res = run_deal(fe(55), n, t, pts);
    const auto& b = res.broadcast_data;
    for (size_t i=0;i<n;++i) {
        Commitment c_a = MerkleTree::commit(res.private_data[i].a_i);
        CHECK(verify_path(b.c_phi, c_a, res.private_data[i].pi_i));
    }
}

// Coherence a + alpha*b = psi 
TEST_CASE("run_deal() - a+alpha*b coherent avec c_psi_Y") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto pts = make_pts(n);
    auto res = run_deal(fe(33), n, t, pts);
    const auto& b = res.broadcast_data;
    const auto& p = res.private_data;
    FieldElement alpha = compute_alpha(b.c_r_Y, b.c_phi);
    for (size_t i=0;i<n;++i) {
        std::vector<FieldElement> psi_row(n);
        for (size_t j=0;j<n;++j)
            psi_row[j] = p[i].a_i[j] + alpha * p[i].b_i[j];
        CHECK(MerkleTree::commit(psi_row) == b.c_psi_Y[i]);
    }
}

// Coherence colonnes psi avec c_psi_X et c 
TEST_CASE("run_deal() - colonnes psi coherentes avec c_psi_X et c") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto pts = make_pts(n);
    auto res = run_deal(fe(11), n, t, pts);
    const auto& b = res.broadcast_data;
    const auto& p = res.private_data;
    FieldElement alpha = compute_alpha(b.c_r_Y, b.c_phi);
    std::vector<Commitment> c_psi_X_recomp(n);
    for (size_t i=0;i<n;++i) {
        std::vector<FieldElement> col(n);
        for (size_t j=0;j<n;++j)
            col[j] = p[j].a_i[i] + alpha * p[j].b_i[i];
        c_psi_X_recomp[i] = MerkleTree::commit(col);
        CHECK(c_psi_X_recomp[i] == b.c_psi_X[i]);
    }
    CHECK(MerkleTree::commit(c_psi_X_recomp) == b.c);
}

// Coherence c_phi depuis les a_{i,j} 
TEST_CASE("run_deal() - c_phi reconstruit depuis les a_{i,j}") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto res = run_deal(fe(77), n, t, make_pts(n));
    const auto& b = res.broadcast_data;
    const auto& p = res.private_data;
    std::vector<Commitment> c_a(n);
    for (size_t i=0;i<n;++i) c_a[i] = MerkleTree::commit(p[i].a_i);
    CHECK(MerkleTree::commit(c_a) == b.c_phi);
}

// Reconstruction du secret par Lagrange (deux niveaux)
// phi(0,0) = secret.  Avec t+1 messages prives :
// phi(x_i,0) = Lagrange({x_j},{a_{i,j}})(0) pour i=0..t
// phi(0,0)   = Lagrange({x_i},{phi(x_i,0)})(0)
TEST_CASE("run_deal() - secret recuperable par double Lagrange (n=4, t=1)") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto pts = make_pts(n);
    FieldElement secret = fe(42);
    auto res = run_deal(secret, n, t, pts);
    const auto& p = res.private_data;

    std::vector<FieldElement> phi_xi_0(t+1);
    for (size_t i=0;i<=t;++i)
        phi_xi_0[i] = lagrange_eval(pts, p[i].a_i, fe(0));

    std::vector<FieldElement> pts_t1(pts.begin(), pts.begin() + static_cast<long>(t+1));
    CHECK(lagrange_eval(pts_t1, phi_xi_0, fe(0)) == secret);
}

TEST_CASE("run_deal() - secret recuperable (n=6, t=2)") {
    init_field(NTL::ZZ(7919));
    size_t n=6, t=2;
    auto pts = make_pts(n);
    FieldElement secret = fe(123);
    auto res = run_deal(secret, n, t, pts);
    const auto& p = res.private_data;

    std::vector<FieldElement> phi_xi_0(t+1);
    for (size_t i=0;i<=t;++i)
        phi_xi_0[i] = lagrange_eval(pts, p[i].a_i, fe(0));

    std::vector<FieldElement> pts_t1(pts.begin(), pts.begin() + static_cast<long>(t+1));
    CHECK(lagrange_eval(pts_t1, phi_xi_0, fe(0)) == secret);
}

TEST_CASE("run_deal() - secret=0 recuperable") {
    init_field(NTL::ZZ(7919));
    size_t n=4, t=1;
    auto pts = make_pts(n);
    auto res = run_deal(fe(0), n, t, pts);
    const auto& p = res.private_data;
    std::vector<FieldElement> phi_xi_0(t+1);
    for (size_t i=0;i<=t;++i)
        phi_xi_0[i] = lagrange_eval(pts, p[i].a_i, fe(0));
    std::vector<FieldElement> pts_t1(pts.begin(), pts.begin() + static_cast<long>(t+1));
    CHECK(lagrange_eval(pts_t1, phi_xi_0, fe(0)) == fe(0));
}

// Parametres invalides
TEST_CASE("run_deal() - x_coords.size() != n : exception") {
    init_field(NTL::ZZ(7919));
    CHECK_THROWS_AS(run_deal(fe(1), 4, 1, make_pts(3)), std::invalid_argument);
}

TEST_CASE("run_deal() - n < t+1 : exception") {
    init_field(NTL::ZZ(7919));
    CHECK_THROWS_AS(run_deal(fe(1), 2, 2, make_pts(2)), std::invalid_argument);
}