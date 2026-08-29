// Compilation : make avss_demo puis ./avss_demo
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <time.h>
#include "deal.hpp"
#include "getshare.hpp" 

// Compteur de hash defini dans merkle.cpp
extern long long hash_count; // 64 bits

//Paramètres par defaut pouvant etre modifies lors de l'execution
long p = 7919;
size_t n = 5; 
size_t t = 2;

//Fonctions utilitaires d'affichage et de conversion pour la demonstration

/*
 * Affiche un Hash_type de 32 octets = 64 caracteres en hex et les ecrit
 * sous forme de (1er caractere hex ... dernier caractere hex) : 35eb5246...4fb3ff2c
*/ 
static std::string short_hex(const Hash_type& d) {
    std::ostringstream s;
    for (size_t i = 0; i < 4; ++i)
        s << std::hex << std::setw(2) << std::setfill('0') << (int)d[i];
    s << "...";
    for (size_t i = 28; i < 32; ++i)
        s << std::hex << std::setw(2) << std::setfill('0') << (int)d[i];
    return s.str();
}

// Affiche un FieldElement comme entier 
static long fe_to_long(const FieldElement& x) {
    return NTL::to_long(NTL::rep(x));
}

// Separateur visuel
static void sep(char c = '-', int w = 56) {
    std::cout << std::string(w, c) << "\n";
}

// Interpolation de Lagrange 
static FieldElement lagrange_eval(
    const std::vector<FieldElement>& xs,
    const std::vector<FieldElement>& ys,
    const FieldElement& x
) {
    FieldElement result;
    NTL::conv(result, 0); // result = 0
    
    for (size_t i = 0; i < xs.size(); ++i) {
        FieldElement num, den;
        NTL::conv(num, 1); // numerateur = 1
        NTL::conv(den, 1); // denominateur = 1
        
        for (size_t j = 0; j < xs.size(); ++j) {
            if (i != j) {
                num *= (x - xs[j]); // num * (x - x_j)
                den *= (xs[i] - xs[j]); // den * (x_i - x_j)
            }
        }
        // result = Sum (y_i * (num / den))
        // Il suffira d'évaluer x en 0 dans le main
        result += ys[i] * (num / den);
    }
    return result;
}

static std::vector<FieldElement> x_pts(size_t n) {
    std::vector<FieldElement> pts; // pts = [0,...,0]
    for (size_t i = 1; i <= n; ++i) { // starts from i = 1 
        // Conversion de i dans N en FieldElement
        FieldElement x; NTL::conv(x, static_cast<long>(i));
        pts.push_back(x); // pts = [1, 2,..., n]
    }
    return pts;
}


int main(int argc, char* argv[]) {
    if (argc >= 2) n = std::atoi(argv[1]); // ASCII to int
    if (argc >= 3) t = std::atoi(argv[2]); 
    if (argc >= 4) p = std::atol(argv[3]); // ASCII to long pour un grand p

    // Vérification de sécurité
    if (t >= n) {
        std::cerr << "Erreur : le seuil t doit etre inferieur a n.\n";
        return 1;
    }

    if (t * 3 >= n) {
    std::cerr << "Erreur : la condition t < n/3 n'est pas respectée.\n";
    return 1;
}
    init_field(NTL::ZZ(p));
    

    FieldElement secret; NTL::conv(secret, 62L); // secret = 62 de type long

    // Affichage des paramètres
    std::cout << " p = " << p << "\n";
    std::cout << " n = " << n << "\n";
    std::cout << " t = " << t << "\n";
    std::cout << " secret s = " << fe_to_long(secret) << "\n";

    std::vector<FieldElement> pts = x_pts(n);
    std::cout << " points d'evaluation : x_i = 1, 2, ..., " << n << "\n";

    clock_t start_time = clock();
    // Execution de l'algorithme Deal

    sep();
    std::cout << " Exécution de l'algorithme Deal : \n";
    DealResult res = run_deal(secret, n, t, pts);
    std::cout << " phi(X,Y) et r(X,Y) de degré <= " << t << " sont générés, tel que phi(0,0) = s\n";
    std::cout << " Les évaluation de phi et r sont représentés par les matrices a_{i,j} et b_{i,j}, calculées en " << n << "x" << n << " = " << n*n << " évaluations\n";
    std::cout << " Les engagements Merkle construits\n";
    std::cout << " // alpha = Hash(c_phi || racine(c_r_X[i]) || racine(c_r_Y[i])) est calculé\n";
    std::cout << " psi = phi + alpha*r est évalué\n";
    std::cout << " La preuve de proximite est generée\n";

    // BroadcastMessage : donnees publiques que tous les P_i recoivent de la part du Dealer
    sep();
    // Extraction directe des résultats du BroadcastMessage avec le format de short_hex()
    const auto& bcast = res.broadcast_data;
    std::cout << " Donnees publiques envoyées par le Dealer dans BroadcastMessage :\n";
    std::cout << " c_phi est la racine Merkle de phi : " << short_hex(bcast.c_phi) << "\n";
    std::cout << " c est la racine Merkle de psi : " << short_hex(bcast.c) << "\n";
    std::cout << " rho est le défi Fiat-Shamir : " << fe_to_long(bcast.pi.rho) << "\n";
    std::cout << " Il y a " << n << " engagements c_r_X[i], c_r_Y[i], c_psi_X[i], c_psi_Y[i]\n";
    std::cout << " Exemple : c_r_X[0] : " << short_hex(bcast.c_r_X[0]) << "\n";
    std::cout << " Exemple : c_psi_Y[0]: " << short_hex(bcast.c_psi_Y[0]) << "\n";

    // PrivateMessages : ce que chaque P_i recoit en prive
    sep();
    std::cout << " Données privees phi(x_i, Y) et r(x_i, Y) que D envoie à chaque P_idans PrivateMessages :\n";
    for (size_t i = 0; i < n; ++i) {
        const auto& private_m = res.private_data[i];
        // Affichage de la ligne b_i 
        std::cout << " P_" << i+1
                  << " : a_{(" << i+1 << ",1)..n} = [";
        for (size_t j = 0; j < n; ++j) {
            std::cout << fe_to_long(private_m.a_i[j]); // private_m.a_i[j] = phi(x_i, Y)
            if (j < n-1) std::cout << ", ";
        }
        std::cout << "]\n";

        // Affichage de la ligne b_i 
        std::cout << "       b_{(" << i+1 << ",1)..n} = [";
        for (size_t j = 0; j < n; ++j) {
            std::cout << fe_to_long(private_m.b_i[j]); // private_m.b_i[j] = r(x_i, Y)
            if (j < n-1) std::cout << ", ";
        }
        std::cout << "]\n";
    }

    // Verification de la preuve de proximite
    sep();
    std::cout << " Vérification de la Preuve de Proximite pour le degré des polynomes et leurs engagements : \n";
    // pi contient l'ouverture ou le chemin juste à cote de chaque feuille c_psi_Y[i] pour chaque P_i
    bool prox_verif = verify_prox(
        bcast.pi, t, pts, bcast.c, bcast.c_psi_Y);
    std::cout << " verify_prox(bcast.pi,t, pts, bcast.c, bcast.c_psi_Y) = "
            // if prox_verif = true : Valide, sinon : Invalide
              << (prox_verif ? "Valide " : "Invalide ") << "\n"; 
              
    // Verification des chemins Merkle pi_i (P_i peut verifier sa part)
    sep();
    std::cout << "Verification des chemins Merkle pi_i\n";
    bool paths_verified = true;
    for (size_t i = 0; i < n; ++i) {
        Commitment c_a_i = MerkleTree::commit(res.private_data[i].a_i); // c_a_i = Hash(a_i)
        bool check_path = verify_path(bcast.c_phi, c_a_i, res.private_data[i].pi_i);
        std::cout << " P_" << i+1 << " : pi_" << i+1 << " dans c_phi = "
                // if check_path = true : Valide, sinon : Invalide
                  << (check_path ? "Valide" : "Invalide") << "\n";
        // Si un seul chemin (check_path) est invalide, paths_verified = false
        if (!check_path) paths_verified = false;
    }

    // Reconstruction theorique du secret depuis le Deal Algo
    sep();
    std::cout << "Vérification de la bonne reconstruction du secret avec une simulation de t+1 participants honnetes :\n";

    // Lagrange sur les t+1 premiers points pour reconstruire les secrets locaux phi(x_i,0) pour chauqe P_i
    FieldElement zero; NTL::conv(zero, 0L);
    std::vector<FieldElement> phi_x_0(t + 1); // phi_xi_Y = [phi(0,x_1), phi(0,x_2), ..., phi(0,x_{t+1})]
    for (size_t i = 0; i <= t; ++i) {
        phi_x_0[i] = lagrange_eval(pts, res.private_data[i].a_i, zero); // phi_xi_Y[0] = phi_xi_Y(x_i, 0)
    }
    // Lagrange sur les t+1 premiers secrets locaux phi_xi_Y[0]pour retrouver phi(0,0)
    std::vector<FieldElement> pts_t1(pts.begin(), pts.begin() + static_cast<long>(t+1));
    FieldElement recovered = lagrange_eval(pts_t1, phi_x_0, zero);
    bool good_secret = (recovered == secret);
    std::cout << "struit est : " << fe_to_long(recovered) << " " 
              << (good_secret ? "Valide" : "Invalide") << "\n\n";

    // Execution de l'algorithme GetShare par le reseau de P_i
    sep();
    std::cout << "Exécution de l'algorithme GetShare par les P_i :\n";
    std::cout << "Ils vérifient les shares reçus et preparent leurs Echos (subshares + preuves Merkle)\n";

    // Initialisation de l'état interne de chaque P_i après reception des shares 
    // avec hond_d = 0, V = {} et ready_count = 0
    std::vector<ParticipantState> reseau;
    for (size_t i = 0; i < n; ++i) {
        reseau.push_back(init_participant(i, n, t, pts));
    }

    // Traitement du message du Dealer
    std::vector<std::vector<EchoMessage>> echos_sub(n);
    bool D_hon = true;
    for (size_t exp = 0; exp < n; ++exp) {
        echos_sub[exp] = receive_deal(reseau[exp], bcast, res.private_data[exp]);
        if (echos_sub[exp].empty()) D_hon = false;
    }
    std::cout <<"hon_d = " << (D_hon ? "1 : Valide, Dealer honnête" : "0 : Erreur, Dealer rejeté") << "\n";

    // Echange des Echos de P_i vers P_j
    sep();
    std::cout << "Echange des messages Echo contenant les subshares de P_i vers P_j : \n";
    std::vector<bool> rdy_2send(n, false); // au début personne n'est pret
    for (size_t exp = 0; exp < n; ++exp) {
        size_t indice_echo = 0;
        // P_i distribue ses Echos à tous les autres P_j
        for (size_t j = 0; j < n; ++j) {
            if (j == exp) continue;
            // Si echo valide, a_i,j et b_I,j ajouté dans V 
            //Retourne true si card(V) > 2t et hon_d=1
            bool signal_ready = receive_echo(reseau[j], echos_sub[exp][indice_echo]);
            if (signal_ready) rdy_2send[j] = true; // Si true, P_j est pret envoie Ready
            indice_echo++;
        }
    }
    std::cout << "Les subshares des Echos sont verifiés et ajoutés à l'ensemble V.\n";
    // a la fin rdy_2send = [F, T, F, T, T], P2, P4 et P5 sont prets, MAIS PAS P1 et P3

    // Accumulation des Ready
    sep();
    std::cout << "Rassemblement de 2t+1 Ready : \n";
    for (int tour = 1; tour <= 2; ++tour) { // 2 tours juste pour l'exmple avec d honnete et t+1 P_i honnetes 
        std::vector<bool> new_ready(n, false);
        for (size_t exp = 0; exp < n; ++exp) {
            if (!rdy_2send[exp]) continue; // On saute les P_i qui ne sont pas prets apres lecture de rdy_2send
            for (size_t j = 0; j < n; ++j) {
                // Chaque P_i pret envoie son Ready à tous les autres P_j 
                if (receive_rdy(reseau[j], exp)) new_ready[j] = true;
            }
        }
        rdy_2send = new_ready; // P1 et P3 recoivent >= t+1 Ready et deviennent prets, donc rdy_2send = [T, F, T, F, F]
    } // il faut un 2eme tout pour rdy_2send = [T, T, T, T, T]
    std::cout << " " << (2*t+1) << "/" << n 
          << " Il y a au moiuns 2t+1 Ready atteints, la reconstruction du secret est possible.\n";

    // Reconstruction des polynomes phi(X,xj) "pas a l'origine du secret" 
    sep();
    std::cout << "P_i reconstruisent phi_i(X) = phi(X, x_i) et r_i(X) = r(X, x_i)\n";
    std::cout << " a partir des subshares phi(x_j, x_i) recues de chaque P_j\n";
    bool good_sub = true;
    std::vector<GetShareResult> final_sub(n); // Contiendra phi_i(X) et r_i(X)
    for (size_t i = 0; i < n; ++i) {
        auto sub = reconstruct(reseau[i]); // P_i interpole les bons phi(xj, xi) in V
        if (!sub) { // Si reconstruct = false, P_i n'a pas assez de subshares valides dans V
            good_sub = false;
        } else {
            final_sub[i] = sub.value(); // Sinon, P_i a reconstruit phi(X, xi) et r(X, xi)
            std::cout << " P_" << i+1 << " a reconstruit ses polynomes phi_" << i+1 << "(X) et r_" << i+1 << "(X)\n";
        }
    }

    // Juste pour s'asurer que phi(0,0) peut etre aussi retrouve avec les phi(X,x_i) de GetShare
    sep();
    std::cout << "Vérification que le secret phi(0,0) est bien retrouvé depuis les phi(X,x_i) de GetShare :\n";
    std::vector<FieldElement> phi_xi_0(t + 1); // vecteur de phi(X,xi) pour i=1,...,t+1
    for (size_t i = 0; i <= t; ++i) {
        FieldElement eval_0;
        NTL::eval(eval_0, final_sub[i].phi_share, zero); // evaluation de phi(0,xi) 
        phi_xi_0[i] = eval_0;
    }
    // Lagrange sur phi(0,xi)
    FieldElement s_found = lagrange_eval(pts_t1, phi_xi_0, zero);
    bool good_s = (s_found == secret); // Comparaison entre le vrai secret et celui qu'on a trouve
    
    std::cout << "Le secret d'origine est : " << fe_to_long(secret) << "\n";
    std::cout << "struit par le réseau est : " << fe_to_long(s_found) << "\n";


    bool all_verif = prox_verif && paths_verified && good_secret && D_hon && good_sub && good_s;
    
    sep();
    if (all_verif) {
        std::cout << " Les protocoles Deal et GetShare fonctionnent\n";
    } else {
        std::cout << " erreur : Une étape du protocole a echoué.\n";
    }

    sep();
    clock_t end_time = clock();
    double time_dif = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    std::cout << "Le temps d'éxécution est de : " << time_dif << " secondes\n";
    std::cout << "Nombre de hachage realisés: " << hash_count << "\n";

    return all_verif ? 0 : 1;
}