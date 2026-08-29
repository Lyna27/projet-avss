// Compilation : make avss_bench puis ./avss_bench n_max t_max

#include "matplotlibcpp.h"
#include "deal.hpp"
#include "getshare.hpp"
#include <iostream>
#include <vector>
#include <ctime>

namespace plt = matplotlibcpp;
extern long long hash_count; // defini dans merkle.cpp

// Pour generer les points d'evaluation
static std::vector<FieldElement> x_pts(size_t n) {
    std::vector<FieldElement> pts; // pts = []
    for (size_t i = 1; i <= n; ++i) {
        FieldElement x; NTL::conv(x, static_cast<long>(i));
        pts.push_back(x);
    }
    return pts;
}

// Separateur visuel
static void sep(char c = '-', int w = 56) {
    std::cout << std::string(w, c) << "\n";
}

int main(int argc, char* argv[]) {

    // Parametres par defaut
    size_t n_max = 100;
    size_t t_max = 30; 
    long p = 7919;

    if (argc >= 2) n_max = std::atoi(argv[1]);
    if (argc >= 3) t_max = std::atoi(argv[2]);
    if (argc >= 4) p = std::atol(argv[3]);

    if (t_max >= n_max) {
        std::cerr << "Erreur : t_max doit etre inferieur a n_max.\n";
        return 1;
    }

    if (t_max * 3 >= n_max) { // si t >= 1/3 n, ça ne marche pas
    std::cerr << "Erreur : la condition t < n/3 n'est pas respectée.\n";
    return 1;
}

    // Ratio t/n fixe a 30% car t ne doit pas depasser un tiers de n
    double ratio = static_cast<double>(t_max) / static_cast<double>(n_max);
    size_t pas = (n_max <= 50) ? 5 : 10; // if n_max <= 50, pas = 5, else pas = 10

    init_field(NTL::ZZ(p));
    FieldElement secret; NTL::conv(secret, 62L); // s = 62

    std::vector<double> n_vals, temps_deal, temps_getshare, hashes_deal, hashes_getshare;

    std::cout << "Paramètre p=" << p << " et ratio t/n=" << ratio << "\n";
    sep();

    for (size_t n = pas; n <= n_max; n += pas) { // nbre de pas adapte a n_max

        size_t t = std::max(size_t(1), static_cast<size_t>(ratio * n));
        if (t >= n) t = n - 1; // pour ne pas depasser sur la courbe

        std::vector<FieldElement> pts = x_pts(n);

        // Deal 
        hash_count = 0;
        clock_t d_start = clock();
        DealResult res = run_deal(secret, n, t, pts);
        clock_t d_end = clock();
        double t_deal = (double)(d_end - d_start) / CLOCKS_PER_SEC; // d_end - d_start nbr de tics / tics par s = s
        long long h_deal = hash_count; // nbr de hash dans Deal

        // GetShare 
        hash_count = 0;
        BenchResult bench = run_getshare(n, t, pts, res);

        // Affichage des deux resultats
        std::cout << " n=" << n << ", t=" << t
                  << " | Deal : " << t_deal << "s, " << h_deal << " hash"
                  << " | GetShare : " << bench.temps << "s, " << bench.nbr_hash << " hash\n";

        // Allignement des resultats pour le graphe
        n_vals.push_back((double)n);
        temps_deal.push_back(t_deal);
        temps_getshare.push_back(bench.temps);
        hashes_deal.push_back((double)h_deal);
        hashes_getshare.push_back((double)bench.nbr_hash);
    }

    // Figure 1 : temps d'exécution de Deal et GetShare en fonction de n
    plt::figure_size(800, 500);
    plt::named_plot("Deal", n_vals, temps_deal,"bo-"); // (nom courbe, x, y, style)
    plt::named_plot("GetShare", n_vals, temps_getshare, "rs-");
    plt::title("Temps d'execution en fonction de n");
    plt::xlabel("n : participants");
    plt::ylabel("Temps : secondes");
    plt::legend();
    plt::grid(true);
    plt::save("avss_temps.png");
    plt::show();

    // Figure 2 : nombre de fonctions de hachages de Deal et GetShare en fonction de n P_i
    plt::figure_size(800, 500);
    plt::named_plot("Deal", n_vals, hashes_deal, "bo-");
    plt::named_plot("GetShare", n_vals, hashes_getshare, "rs-");
    plt::title("Nombre de hachages en fonction de n");
    plt::xlabel("n : participants");
    plt::ylabel("Nombre de hachages");
    plt::legend();
    plt::grid(true);
    plt::save("avss_hachages.png");
    plt::show();

    // Figure 3 : Deal vs GetShare pour participant
    std::vector<double> hashes_par_Pi;
    for (size_t k = 0; k < n_vals.size(); ++k) {
        hashes_par_Pi.push_back(hashes_getshare[k] / n_vals[k]);
    }

    plt::figure_size(800, 500);
    plt::named_plot("Deal (Dealer seul)", n_vals, hashes_deal, "bo-");
    plt::named_plot("GetShare (par P_i)", n_vals, hashes_par_Pi, "rs-");
    plt::title("Complexite par acteur : Deal vs GetShare par P_i");
    plt::xlabel("n : participants");
    plt::ylabel("Nombre de hachages");
    plt::legend();
    plt::grid(true);
    plt::save("avss_complexite_par_acteur.png");
    plt::show();

    return 0;
}