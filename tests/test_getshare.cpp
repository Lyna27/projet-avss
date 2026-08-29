#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "getshare.hpp"
#include "deal.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_p.h>
#include <optional>

static FieldElement fe(long valeur) {
    FieldElement x; NTL::conv(x, valeur); return x;
}

static std::vector<FieldElement> creer_points(size_t n) {
    std::vector<FieldElement> pts;
    for (size_t i = 1; i <= n; ++i) pts.push_back(fe(static_cast<long>(i)));
    return pts;
}

static FieldElement evaluer(const Uni_Poly& poly, const FieldElement& point) {
    FieldElement res; NTL::eval(res, poly, point); return res;
}

static FieldElement interpoler_en(
    const std::vector<FieldElement>& xs,
    const std::vector<FieldElement>& ys,
    const FieldElement& cible
) {
    NTL::vec_ZZ_p xa, ya;
    xa.SetLength(static_cast<long>(xs.size()));
    ya.SetLength(static_cast<long>(ys.size()));
    for (size_t k = 0; k < xs.size(); ++k) {
        xa[static_cast<long>(k)] = xs[k];
        ya[static_cast<long>(k)] = ys[k];
    }
    Uni_Poly f; NTL::interpolate(f, xa, ya);
    FieldElement res; NTL::eval(res, f, cible); return res;
}

struct ResultatSimulation {
    bool tout_reussi;
    std::vector<GetShareResult> parts;
};

static ResultatSimulation simuler_reseau_honnete(
    size_t n_Pi,
    size_t seuil,
    const FieldElement& valeur_secrete
) {
    auto points = creer_points(n_Pi);
    DealResult deal = run_deal(valeur_secrete, n_Pi, seuil, points);

    std::vector<ParticipantState> reseau;
    reseau.reserve(n_Pi);
    for (size_t qui = 0; qui < n_Pi; ++qui)
        reseau.push_back(init_participant(qui, n_Pi, seuil, points));

    std::vector<std::vector<EchoMessage>> echo_2snd(n_Pi);
    for (size_t expediteur = 0; expediteur < n_Pi; ++expediteur) {
        echo_2snd[expediteur] = receive_deal(
            reseau[expediteur], deal.broadcast_data, deal.private_data[expediteur]
        );
        if (echo_2snd[expediteur].size() != n_Pi - 1) return {false, {}};
    }

    std::vector<bool> rdy_2send(n_Pi, false);
    for (size_t expediteur = 0; expediteur < n_Pi; ++expediteur) {
        size_t indice_echo = 0;
        for (size_t destinataire = 0; destinataire < n_Pi; ++destinataire) {
            if (destinataire == expediteur) continue;
            bool signal_ready = receive_echo(reseau[destinataire], echo_2snd[expediteur][indice_echo]);
            if (signal_ready) rdy_2send[destinataire] = true;
            indice_echo++;
        }
    }

    for (int tour = 0; tour < 3; ++tour) {
        std::vector<bool> new_ready(n_Pi, false);
        for (size_t expediteur = 0; expediteur < n_Pi; ++expediteur) {
            if (!rdy_2send[expediteur]) continue;
            for (size_t destinataire = 0; destinataire < n_Pi; ++destinataire) {
                if (destinataire == expediteur) continue;
                bool signal_ready = receive_rdy(reseau[destinataire], expediteur);
                if (signal_ready) new_ready[destinataire] = true;
            }
        }
        rdy_2send = new_ready;
    }

    ResultatSimulation resultat;
    resultat.tout_reussi = true;
    for (size_t qui = 0; qui < n_Pi; ++qui) {
        auto share = reconstruct(reseau[qui]);
        if (!share) {
            resultat.tout_reussi = false;
        } else {
            resultat.parts.push_back(*share);
        }
    }
    return resultat;
}

static std::pair<ParticipantState, EchoMessage> preparer_etat_et_echo_de_p1(
    const DealResult& deal, size_t n, size_t seuil, const std::vector<FieldElement>& points
) {
    auto etat_p0 = init_participant(0, n, seuil, points);
    receive_deal(etat_p0, deal.broadcast_data, deal.private_data[0]);

    auto etat_p1 = init_participant(1, n, seuil, points);
    auto echos_de_p1 = receive_deal(etat_p1, deal.broadcast_data, deal.private_data[1]);
    return {etat_p0, echos_de_p1[0]};
}

TEST_CASE("init_participant - état initial propre") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    ParticipantState etat = init_participant(2, n, seuil, points);
    CHECK(etat.my_index == 2);
    CHECK(etat.hon_d == false);
}

TEST_CASE("receive_deal - dealer honnête") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    DealResult deal = run_deal(fe(42), n, seuil, points);
    auto etat_p0 = init_participant(0, n, seuil, points);
    auto echos = receive_deal(etat_p0, deal.broadcast_data, deal.private_data[0]);
    CHECK(echos.size() == n - 1);
    CHECK(etat_p0.hon_d);
}

TEST_CASE("receive_deal - ProximityProof falsifiée -> abort") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    DealResult deal = run_deal(fe(42), n, seuil, points);

    BroadcastMessage broadcast_falsifie = deal.broadcast_data;
    broadcast_falsifie.pi.rho = broadcast_falsifie.pi.rho + fe(1);

    auto etat = init_participant(0, n, seuil, points);
    auto echos = receive_deal(etat, broadcast_falsifie, deal.private_data[0]);
    CHECK(echos.empty());
}

TEST_CASE("receive_deal - valeurs a_i falsifiées -> abort") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    DealResult deal = run_deal(fe(99), n, seuil, points);

    PrivateMessage share_trafique = deal.private_data[0];
    share_trafique.a_i[0] = share_trafique.a_i[0] + fe(1);

    auto etat = init_participant(0, n, seuil, points);
    auto echos = receive_deal(etat, deal.broadcast_data, share_trafique);
    CHECK(echos.empty());
}

TEST_CASE("receive_deal - valeurs b_i falsifiées -> abort") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    DealResult deal = run_deal(fe(55), n, seuil, points);

    PrivateMessage share_trafique = deal.private_data[0];
    share_trafique.b_i[0] = share_trafique.b_i[0] + fe(1);

    auto etat = init_participant(0, n, seuil, points);
    auto echos = receive_deal(etat, deal.broadcast_data, share_trafique);
    CHECK(echos.empty());
}

TEST_CASE("receive_echo - Echo valide ajouté à V") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    DealResult deal = run_deal(fe(42), n, seuil, points);
    auto [etat_p0, echo_de_p1] = preparer_etat_et_echo_de_p1(deal, n, seuil, points);
    bool signal_ready = receive_echo(etat_p0, echo_de_p1);
    CHECK(etat_p0.V.size() == 2); 
    CHECK(signal_ready == false);
}

TEST_CASE("receive_echo - valeur b_ij falsifiée -> rejeté") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    auto points = creer_points(n);
    DealResult deal = run_deal(fe(42), n, seuil, points);
    auto [etat_p0, echo_trafique] = preparer_etat_et_echo_de_p1(deal, n, seuil, points);
    
    echo_trafique.b_ij = echo_trafique.b_ij + fe(1);

    receive_echo(etat_p0, echo_trafique);
    CHECK(etat_p0.V.size() == 1);
}

TEST_CASE("Simulation honnête - secret=42 récupérable") {
    init_field(NTL::ZZ(7919));
    size_t n = 4, seuil = 1;
    FieldElement valeur_secrete = fe(42);
    auto sim = simuler_reseau_honnete(n, seuil, valeur_secrete);
    REQUIRE(sim.tout_reussi);
}