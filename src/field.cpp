#include "field.hpp"
#include <fstream>
#include <stdexcept>
#include <vector>

namespace {

bool is_init = false;

/* 
 * Fournit l’entropie de depart du generateur pseudo-aleatoire
 * interne de NTL en utilisant /dev/urandom, et renvoie un message 
 * d'erreur si l'ouverture ou la lecture echoue.
 */
void seed_rng() {
    constexpr std::size_t SEED_BYTES = 32;
    std::ifstream urandom("/dev/urandom", std::ios::binary);
    if (!urandom) {
        throw std::runtime_error(
            "erreur : impossible d'ouvrir /dev/urandom.");
    }
    std::vector<unsigned char> seed(SEED_BYTES);
    urandom.read(reinterpret_cast<char*>(seed.data()),
                 static_cast<std::streamsize>(SEED_BYTES));
    if (!urandom) {
        throw std::runtime_error(
            "erreur : lecture de /dev/urandom incomplète.");
    }
    NTL::SetSeed(seed.data(), static_cast<long>(seed.size()));
}

} // fin du namespace anonyme

void init_field(const NTL::ZZ& prime) {
    NTL::ZZ_p::init(prime);
    seed_rng();
    is_init = true;
}

FieldElement random_element() {
    if (!is_init) {
        throw std::runtime_error(
            "erreur : le corps n'est pas correctement initialisé,");
    }
    FieldElement x;
    NTL::random(x);
    return x;
}