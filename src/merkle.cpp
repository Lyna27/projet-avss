#include "merkle.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <cstring>
#include <stdexcept>
#include <vector>

long long hash_count = 0;

namespace {

/* Implementation interne de SHA-256 suivant la norme FIPS 180-4.
 * Suivre le lien suivant pour la specification officielle du NIST :
 * lien : https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.180-4.pdf
 */
static const uint32_t SHA256_K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
    0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
    0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
    0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

/* 
 * Vecteur d'initialisation H^(0), compose de 8 valeurs 
 * de 32 bits, utilise au tout debut du hachage.
 */
static const uint32_t SHA256_IV[8] = {
    0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
    0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
};

/* 
 * Fonction qui applique une rotation circulaire de n bits
 * vers la droite sur un mot de 32 bits.
 */
inline uint32_t rotr32(uint32_t x, unsigned n) {
    return (x >> n) | (x << (32u - n));
}

//Fonction de compression de SHA-256 
void sha256_block(uint32_t state[8], const uint8_t block[64]) {
    uint32_t w[64]; // message schedule : W_0, W_1, ..., W_63

    // Initialisation des 16 premiers mots W_0, ..., W_15 
    for (int i = 0; i < 16; ++i) { 
        // Recomposition de 4 octets en un mot de 32 bits avec Big-Endian
        w[i] = (uint32_t(block[4*i+0]) << 24) // AA000000
              | (uint32_t(block[4*i+1]) << 16) // +00BB0000
              | (uint32_t(block[4*i+2]) <<  8) // +0000CC00
              |  uint32_t(block[4*i+3]);  // + 000000DD = AABBCCDD
    }

    /* 
     * Expansion du message pour avoir W_16, ..., W_63 grace a la
     * formule : W_i = g1(W_{i-2}) + W_{i-7} + g0(W_{i-15}) + W_{i-16}
     * Pour avoir en tout 64 mots pour les prochains 64 tours de compression. 
     */
    for (int i = 16; i < 64; ++i) {
        uint32_t g0 = rotr32(w[i-15], 7) ^ rotr32(w[i-15], 18) ^ (w[i-15] >>  3);
        uint32_t g1 = rotr32(w[i- 2], 17) ^ rotr32(w[i- 2], 19) ^ (w[i- 2] >> 10);
        w[i] = g1 + w[i-7] + g0 + w[i-16];
    }

    // Initialisation des variables de travail avec l'etat actuel
    uint32_t a=state[0], b=state[1], c=state[2], d=state[3];
    uint32_t e=state[4], f=state[5], g=state[6], h=state[7];

    /* 
     * Boucle de compression sur 64 tours qui melange, 
     * additionne et xor les variables de travail aux
     * messages w[i] et les constantes SHA256_K[i] mod 2^32.
     */
    for (int i = 0; i < 64; ++i) {
        uint32_t s1  = rotr32(e, 6)  ^ rotr32(e, 11) ^ rotr32(e, 25);

        /* 
         * Choix entre f et g selon le bit de e.
         * ch = (e AND f) XOR (NOT e AND g)
         * Si le bit de e est 1, on prend le bit de f correspondant 
         * a l'indice du bit de e, sinon on prend le bit de g.
         */
        uint32_t ch  = (e & f) ^ (~e & g);
        uint32_t t1  = h + s1 + ch + SHA256_K[i] + w[i];
        uint32_t s0  = rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22);

        // dans le vecteur (a,b,c) on prend le bit majoritaire 
        // 0 ou 1 entre a,b,c pour chaque position
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2  = s0 + maj;

        //Mise a jour des variables de travail pour le tour suivant
        h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }

    // Addition de l'ancien etat avec les variables de travail apres 64 tours
    state[0]+=a; state[1]+=b; state[2]+=c; state[3]+=d;
    state[4]+=e; state[5]+=f; state[6]+=g; state[7]+=h;
}

// Fonction SHA-256 complete avec padding 
Hash_type sha256_raw(const uint8_t* data, size_t len) {
    hash_count++;
    uint32_t state[8];
    // Initialisation du tableau des etats avec les valeurs de H^(0)
    for (int i = 0; i < 8; ++i) state[i] = SHA256_IV[i];

    /* 
     * Appel de sha256_block() sur tous les blocs d'au moins 64 octets
     * La boucle s'arrete s'il en reste moins a lire, laissant le reste  
     * dumessage pour l'etape du padding.
     */
    size_t offset = 0;
    for (; offset + 64 <= len; offset += 64) {
        sha256_block(state, data + offset);
    }

    /* 
     * Un tableau vide "buf" de 128 octets est cree.
     * Les octets restants "tail" sont copies depuis la 
     * position de lecture (data + offset) vers le debut de buf.
     * La fin du message est marquee par un bit 1 suivi de zeros 
     * pour completer le bloc.
     */
    uint8_t buf[128] = {}; 
    size_t tail = len - offset; 
    std::memcpy(buf, data + offset, tail); 
    buf[tail] = 0x80u; 

    /* 
     * Ajout de la taille totale du message original en bits.
     * La norme exige de placer cette taille sur les 8 derniers octets.
     * Si tail < 56, il y a assez de place dans le bloc actuel (64 octets).
     * Sinon, le padding deborde sur un deuxieme bloc, et la taille 
     * est ecrite a la fin de ce deuxieme bloc (128 octets au total).
     */
    uint64_t bit_len = static_cast<uint64_t>(len) * 8u;
    if (tail < 56) {
        // Tient dans un seul bloc final de 64 octets
        for (int j = 0; j < 8; ++j)
            buf[63 - j] = static_cast<uint8_t>(
                bit_len >> (8u * static_cast<unsigned>(j)));
        sha256_block(state, buf);
    } else {
        // Necessite deux blocs finaux de 128 octets
        for (int j = 0; j < 8; ++j)
            buf[127 - j] = static_cast<uint8_t>(
                bit_len >> (8u * static_cast<unsigned>(j)));
        sha256_block(state, buf);
        sha256_block(state, buf + 64);
    }

    /* 
     * Conversion de l'etat final en tableau d'octets "digest".
     * Le tableau state contient 8 mots de 32 bits. 
     * On decoupe chaque mot en 4 octets de 8 bits pour obtenir 
     * le resultat final standardise de 32 octets donc 256 bits.
     */
    Hash_type digest;
    for (int j = 0; j < 8; ++j) {
        digest[static_cast<size_t>(4*j+0)] = static_cast<uint8_t>(state[j] >> 24);
        digest[static_cast<size_t>(4*j+1)] = static_cast<uint8_t>(state[j] >> 16);
        digest[static_cast<size_t>(4*j+2)] = static_cast<uint8_t>(state[j] >>  8);
        digest[static_cast<size_t>(4*j+3)] = static_cast<uint8_t>(state[j]);
    }
    return digest;
}

// Fonction qui calcule la plus petite puissance de 2 superieure ou egale a n.
size_t next_pow2(size_t n) {
    size_t p = 1;
    // Left shift : p double a chaque iteration.
    while (p < n) p = p * 2;
    return p;
}

} // Fin du namespace anonyme

Hash_type hash_bytes(const uint8_t* data, size_t len) {
    return sha256_raw(data, len);
}

Hash_type hash_element(const FieldElement& x) {
    // Conversion de l'element de F_p en ZZ
    NTL::ZZ z = NTL::rep(x);  

    // Calcul de log_2(p)/8 = nbr d'octets necessaires pour representer p
    long mod_bytes = NTL::NumBytes(NTL::ZZ_p::modulus());
    if (mod_bytes <= 0) mod_bytes = 1; 

    // Creation d'un tableau d'octets de taille fixe pre-rempli de 0
    std::vector<uint8_t> buf(static_cast<size_t>(mod_bytes), 0u);
    NTL::BytesFromZZ(buf.data(), z, mod_bytes); // Ajout de z dans buf

    // Hachage final du tableau d'octets : 32 octets en sortie
    return hash_bytes(buf.data(), buf.size());
}

Hash_type hash_concat(const Hash_type& left, const Hash_type& right) {
    uint8_t buf[HASH_SIZE * 2]; 
    //Concatenation de left et right dans buf avant le hachage
    std::memcpy(buf, left.data(),  HASH_SIZE);
    std::memcpy(buf + HASH_SIZE, right.data(), HASH_SIZE);
    return hash_bytes(buf, sizeof(buf)); // hash de la concatenation
}


// MerkleTree

namespace MerkleTree {

Commitment Tree::root() const {
    // Verifie si l'arbre et la racine sont vide
    if (levels.empty() || levels.back().empty()) {
        throw std::logic_error(
            "erreur : L'arbre est vide ou mal construit.");
    }
    return levels.back()[0]; // recupere la racine
}

MerklePath Tree::path(size_t position) const {
    if (levels.empty()) {
        throw std::logic_error(
            "erreur : L'arbre est vide.");
    }
    if (position >= levels[0].size()) {
        throw std::out_of_range(
            "erreur : La position est hors borne.");
    }

    MerklePath p; // va contenir les hashs freres lors la montee
    size_t idx = position;

    for (size_t lvl = 0; lvl + 1 < levels.size(); ++lvl) {
        // noeud courant gauche : idx pair 
        if (idx % 2 == 0) { 
            size_t sib = idx + 1; // right sibling (true)

            // ajoute le hash de sib a p
            p.siblings.push_back(levels[lvl][sib]); 
            p.directions.push_back(true);

        // noeud courant droit : idx impair 
        } else {
            // idx-1 left sibling (false)
            p.siblings.push_back(levels[lvl][idx - 1]);
            p.directions.push_back(false);
        }
        idx >>= 1; // passe au niveau superieur
    }
    return p; // le chemin complet
}

Tree build(const std::vector<Hash_type>& leaves) {
    if (leaves.empty()) { // Test de verification
        throw std::invalid_argument(
            "erreur : liste de feuilles vide.");
    }

    Tree tree;

    // Padding des feuilles a la prochaine puissance de 2 avec des 0 
    size_t padded = next_pow2(leaves.size());

    // stocke les feuilles + padding (les 0) dans le tableau "current"
    std::vector<Hash_type> current(padded, Hash_type{});
    for (size_t i = 0; i < leaves.size(); ++i) current[i] = leaves[i];

    // Le prmeier niveau "current" = levels[0] est insere dans l'arbre
    tree.levels.push_back(current);

    // Construction niveau par niveau tant qu'il ya plus d'un noeud 
    // jusqu'au dernier qui est la racine.
    while (current.size() > 1) {
        std::vector<Hash_type> next;
        // Taille niveau(n+1) = (1/2) Taille (niveaun)
        next.reserve(current.size() / 2);

        // Concatenation avec hash_concat()
        for (size_t i = 0; i < current.size(); i += 2) { 
            next.push_back(hash_concat(current[i], current[i + 1]));
        }
        //Resultat stocke dans tree
        tree.levels.push_back(next);
        current = std::move(next);
    }

    return tree;
}

Tree tree_field(const std::vector<FieldElement>& values) {
    if (values.empty()) {
        throw std::invalid_argument(
            "erreur : liste de valeurs vide.");
    }
    std::vector<Hash_type> leaves;
    leaves.reserve(values.size());
    for (const auto& v : values) {
        leaves.push_back(hash_element(v));
    }
    return build(leaves);
}

Tree tree_commits(const std::vector<Commitment>& values) {
    if (values.empty()) {
        throw std::invalid_argument(
            "erreur : liste d'engagements polynomiaux vide.");
    }
    return build(values); 
}

/* 
 * Fonctions utilitaires commit() pour extraire la racine sans stocker
 * l'arbre en entier quel que soit le type de donnees manipulees.
 */
Commitment commit(const std::vector<FieldElement>& values) {
    return tree_field(values).root();
}

Commitment commit(const std::vector<Commitment>& values) {
    return tree_commits(values).root();
}

} // Fin du namespace MerkleTree

bool verify_path(const Commitment& root, const Hash_type& leaf, const MerklePath& path) {
    if (path.siblings.size() != path.directions.size()) {
        return false; // chemin malforme
    }
    Hash_type current = leaf;
    for (size_t i = 0; i < path.siblings.size(); ++i) {

        // Comme Hash(A || B) != Hash(B || A) :
        if (path.directions[i]) {
            current = hash_concat(current, path.siblings[i]);
        } else {
            current = hash_concat(path.siblings[i], current);
        }
    }
    return current == root;
}