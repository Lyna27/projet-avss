// test_merkle.cpp -- tests unitaires pour merkle.cpp

// Compilation (avec NTL installe) :
// g++ -std=c++17 field.cpp merkle.cpp test_merkle.cpp -o test_merkle -lntl -lgmp
// ./test_merkle

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include "merkle.hpp"
#include "field.hpp"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <string>


// Convertit un Hash_type en chaine hexadecimale (lisibilite des echecs).
static std::string to_hex(const Hash_type& d) {
    std::ostringstream s;
    for (auto b : d)
        s << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b);
    return s.str();
}

// Cree un Hash_type dont tous les octets valent `val`.
static Hash_type fill_digest(uint8_t val) {
    Hash_type d;
    d.fill(val);
    return d;
}

// Cree un FieldElement depuis un long via NTL::conv 
static FieldElement fe(long v) {
    FieldElement x;
    NTL::conv(x, v);
    return x;
}

// hash_bytes()

TEST_CASE("hash_bytes() : SHA256('') vecteur NIST") {
    // Longueur 0 : on passe un pointeur valide mais len=0 .
    uint8_t dummy = 0;
    Hash_type h = hash_bytes(&dummy, 0);
    CHECK(to_hex(h) == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}

TEST_CASE("hash_bytes() : SHA256('abc') confirme par Python/OpenSSL") {
    const uint8_t abc[] = {'a', 'b', 'c'};
    Hash_type h = hash_bytes(abc, 3);
    CHECK(to_hex(h) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}

TEST_CASE("hash_bytes() : SHA256(56 octets) vecteur NIST B.2 (chemin deux blocs)") {
    // Ce vecteur active le chemin de padding a deux blocs (tail >= 56).
    const char* msg = "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    Hash_type h = hash_bytes(reinterpret_cast<const uint8_t*>(msg), 56);
    CHECK(to_hex(h) == "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
}

TEST_CASE("hash_bytes() est deterministe : meme entree => meme sortie") {
    const uint8_t data[] = {0xDE, 0xAD, 0xBE, 0xEF};
    CHECK(hash_bytes(data, 4) == hash_bytes(data, 4));
}

TEST_CASE("hash_bytes() : entrees differentes => sorties differents") {
    const uint8_t a[] = {1, 2, 3};
    const uint8_t b[] = {1, 2, 4}; // seul le 3eme octet differe
    CHECK(hash_bytes(a, 3) != hash_bytes(b, 3));
}

TEST_CASE("hash_bytes() : longueurs differentes => sorties differents") {
    const uint8_t data[] = {0x61, 0x62, 0x63};
    CHECK(hash_bytes(data, 1) != hash_bytes(data, 2));
    CHECK(hash_bytes(data, 2) != hash_bytes(data, 3));
}

// hash_element()


TEST_CASE("hash_element() est deterministe") {
    init_field(NTL::ZZ(101));
    CHECK(hash_element(fe(42)) == hash_element(fe(42)));
    CHECK(hash_element(fe(0)) == hash_element(fe(0)));
}

TEST_CASE("hash_element() est non commutatif") {
    init_field(NTL::ZZ(101));
    // Verifier toutes les paires de {0,1,2,3,4,5}
    for (long v1 = 0; v1 < 5; ++v1)
        for (long v2 = v1 + 1; v2 < 6; ++v2)
            CHECK(hash_element(fe(v1)) != hash_element(fe(v2)));
}

TEST_CASE("hash_element() : reste deterministe pour 0") {
    init_field(NTL::ZZ(101));
    Hash_type h = hash_element(fe(0));
    CHECK(h == hash_element(fe(0)));
}

TEST_CASE("hash_element() : encodage depend du module") {
    // fe(5) mod 101 s'encode sur 1 octet ; fe(5) mod 7919 s'encode sur 2 octets.
    // Les digests peuvent donc differer -- l'important est la determinisme.
    init_field(NTL::ZZ(101));
    Hash_type h1 = hash_element(fe(5));
    CHECK(h1 == hash_element(fe(5))); // deterministe avec p=101

    init_field(NTL::ZZ(7919));
    Hash_type h2 = hash_element(fe(5));
    CHECK(h2 == hash_element(fe(5))); // deterministe avec p=7919
}

// hash_concat()

TEST_CASE("hash_concat() est egal a SHA256(left || right)") {
    const uint8_t raw_a[] = {0x11, 0x22};
    const uint8_t raw_b[] = {0x33, 0x44};
    Hash_type ha = hash_bytes(raw_a, 2);
    Hash_type hb = hash_bytes(raw_b, 2);

    uint8_t cat[HASH_SIZE * 2];
    std::memcpy(cat, ha.data(), HASH_SIZE);
    std::memcpy(cat + HASH_SIZE, hb.data(), HASH_SIZE);
    Hash_type expected = hash_bytes(cat, sizeof(cat));

    CHECK(hash_concat(ha, hb) == expected);
}

TEST_CASE("hash_concat() est non commutatif") {
    Hash_type ha = fill_digest(0xAA);
    Hash_type hb = fill_digest(0xBB);
    CHECK(hash_concat(ha, hb) != hash_concat(hb, ha));
}

TEST_CASE("hash_concat() est deterministe") {
    Hash_type ha = fill_digest(0x01);
    Hash_type hb = fill_digest(0x02);
    CHECK(hash_concat(ha, hb) == hash_concat(ha, hb));
}


// MerkleTree::build() + Tree::root()

TEST_CASE("build() : pour n=1, racine == feuille") {
    Hash_type leaf = fill_digest(0x42);
    auto tree = MerkleTree::build({leaf});
    CHECK(tree.levels.size() == 1);
    CHECK(tree.root() == leaf);
}

TEST_CASE("build() : pour n=2, 2 niveaux, racine == hash_concat(l0, l1)") {
    Hash_type l0 = fill_digest(0x01);
    Hash_type l1 = fill_digest(0x02);
    auto tree = MerkleTree::build({l0, l1});
    CHECK(tree.levels.size() == 2);
    CHECK(tree.levels[0].size() == 2);
    CHECK(tree.levels[1].size() == 1);
    CHECK(tree.root() == hash_concat(l0, l1));
}

TEST_CASE("build() : pour n=4, 3 niveaux, racine calculee a la main") {
    Hash_type l0=fill_digest(0xAA), l1=fill_digest(0xBB);
    Hash_type l2=fill_digest(0xCC), l3=fill_digest(0xDD);
    auto tree = MerkleTree::build({l0, l1, l2, l3});

    CHECK(tree.levels.size() == 3);
    CHECK(tree.levels[0].size() == 4);
    CHECK(tree.levels[1].size() == 2);
    CHECK(tree.levels[2].size() == 1);

    // Racine = hash(hash(l0,l1) || hash(l2,l3))
    Hash_type h01 = hash_concat(l0, l1);
    Hash_type h23 = hash_concat(l2, l3);
    Hash_type root = hash_concat(h01, h23);
    CHECK(tree.root() == root);
}

TEST_CASE("build() : padding a 4 feuilles") {
    std::vector<Hash_type> leaves;
    for (uint8_t i = 0; i < 3; ++i) leaves.push_back(fill_digest(i));
    auto tree = MerkleTree::build(leaves);
    CHECK(tree.levels[0].size() == 4); // padde a 4
    CHECK(tree.levels.size() == 3);
}

TEST_CASE("build() : padding a 8 feuilles, 4 niveaux") {
    std::vector<Hash_type> leaves;
    for (uint8_t i = 0; i < 5; ++i) leaves.push_back(fill_digest(i));
    auto tree = MerkleTree::build(leaves);
    CHECK(tree.levels[0].size() == 8);
    CHECK(tree.levels.size() == 4);
}

TEST_CASE("build() : liste vide") {
    CHECK_THROWS_AS(MerkleTree::build({}), std::invalid_argument);
}

// Tree::path() + verify_path()

TEST_CASE("path() + verify() : si chemin vide, feuille == racine") {
    Hash_type leaf = fill_digest(0x11);
    auto tree = MerkleTree::build({leaf});
    MerklePath p = tree.path(0);
    CHECK(p.siblings.empty());
    CHECK(p.directions.empty());
    CHECK(verify_path(tree.root(), leaf, p));
}

TEST_CASE("path() + verify() : directions de l0 et l1 sont correctes") {
    Hash_type l0 = fill_digest(0x01);
    Hash_type l1 = fill_digest(0x02);
    auto tree = MerkleTree::build({l0, l1});
    Commitment root = tree.root();

    // Feuille 0 (index pair) : sibling l1 est a droite alors direction = true
    MerklePath p0 = tree.path(0);
    CHECK(p0.siblings.size() == 1);
    CHECK(p0.siblings[0] == l1);
    CHECK(p0.directions[0] == true);
    CHECK(verify_path(root, l0, p0));

    // Feuille 1 (index impair) : sibling l0 est a gauche alors direction = false
    MerklePath p1 = tree.path(1);
    CHECK(p1.siblings.size() == 1);
    CHECK(p1.siblings[0] == l0);
    CHECK(p1.directions[0] == false);
    CHECK(verify_path(root, l1, p1));
}

TEST_CASE("path() + verify() : les 4 feuilles se verifient") {
    std::vector<Hash_type> leaves;
    for (uint8_t i = 0; i < 4; ++i) leaves.push_back(fill_digest(i));
    auto tree = MerkleTree::build(leaves);
    Commitment root = tree.root();
    for (size_t i = 0; i < 4; ++i) {
        MerklePath p = tree.path(i);
        CHECK(p.siblings.size() == 2); // log2(4) = 2 niveaux de sibling
        CHECK(verify_path(root, leaves[i], p));
    }
}

TEST_CASE("path() + verify() : les 5 feuilles se verifient malgré un padding)") {
    std::vector<Hash_type> leaves;
    for (uint8_t i = 0; i < 5; ++i) leaves.push_back(fill_digest(static_cast<uint8_t>(i * 10)));
    auto tree = MerkleTree::build(leaves);
    Commitment root = tree.root();
    for (size_t i = 0; i < 5; ++i) {
        CHECK(verify_path(root, leaves[i], tree.path(i)));
    }
}

TEST_CASE("verify_path() : feuille falsifiée detectée") {
    std::vector<Hash_type> leaves;
    for (uint8_t i = 0; i < 4; ++i) leaves.push_back(fill_digest(i));
    auto tree = MerkleTree::build(leaves);
    Commitment root = tree.root();

    Hash_type fake = leaves[0];
    fake[0] ^= 0xFF; // un seul bit flip suffit
    CHECK(!verify_path(root, fake, tree.path(0)));
}

TEST_CASE("verify_path() : chemin[i] ne correspond pas à feuille[j])") {
    std::vector<Hash_type> leaves;
    for (uint8_t i = 0; i < 4; ++i) leaves.push_back(fill_digest(i));
    auto tree = MerkleTree::build(leaves);
    Commitment root = tree.root();

    // Utiliser le chemin de la feuille 0 pour verifier la feuille 1 doit echouer
    CHECK(!verify_path(root, leaves[1], tree.path(0)));
    CHECK(!verify_path(root, leaves[0], tree.path(1)));
    CHECK(!verify_path(root, leaves[3], tree.path(0)));
}

TEST_CASE("verify_path() : chemin malformé rejete") {
    Hash_type leaf = fill_digest(0x42);
    auto tree = MerkleTree::build({leaf, fill_digest(0x43)});

    MerklePath bad;
    bad.siblings = {leaf};
    bad.directions = {}; // 1 sibling mais 0 directions : incoherent
    CHECK(!verify_path(tree.root(), leaf, bad));
}

TEST_CASE("path() : position hors borne") {
    auto tree = MerkleTree::build({fill_digest(0x01), fill_digest(0x02)});
    CHECK_THROWS_AS(tree.path(100), std::out_of_range);
}

// MerkleTree::tree_field()

TEST_CASE("tree_field() : coherent avec build(hash_element())") {
    init_field(NTL::ZZ(101));
    std::vector<FieldElement> vals = {fe(10), fe(20), fe(30), fe(40)};

    // Via tree_field
    Commitment root1 = MerkleTree::tree_field(vals).root();

    // Via build() avec les feuilles precalculees manuellement
    std::vector<Hash_type> leaves;
    for (const auto& v : vals) leaves.push_back(hash_element(v));
    Commitment root2 = MerkleTree::build(leaves).root();

    CHECK(root1 == root2);
}

TEST_CASE("tree_field() : toutes les feuilles sont verifiables") {
    init_field(NTL::ZZ(101));
    std::vector<FieldElement> vals = {fe(1), fe(2), fe(3), fe(4), fe(5)};
    auto tree = MerkleTree::tree_field(vals);
    Commitment root = tree.root();
    for (size_t i = 0; i < vals.size(); ++i) {
        Hash_type leaf = hash_element(vals[i]);
        CHECK(verify_path(root, leaf, tree.path(i)));
    }
}

TEST_CASE("tree_field() : liste vide") {
    init_field(NTL::ZZ(101));
    CHECK_THROWS_AS(
        MerkleTree::tree_field(std::vector<FieldElement>{}),
        std::invalid_argument);
}

// MerkleTree::tree_commits()

TEST_CASE("tree_commits() : (Commitment == Hash_type)") {
    std::vector<Commitment> comms = {
        fill_digest(0xAA), fill_digest(0xBB), fill_digest(0xCC)
    };
    CHECK(MerkleTree::tree_commits(comms).root() ==
          MerkleTree::build(comms).root());
}

TEST_CASE("tree_commits() : MT.Commit sur des c_a") {
    // Simule : c_phi = MT.Commit({c_a(x_i,Y)}_i) pour n=4 participants.
    // Verifie que le chemin vers le 2e commitment est valide dans c_phi.
    init_field(NTL::ZZ(101));
    std::vector<Commitment> c_a_commits;
    for (long i = 1; i <= 4; ++i)
        c_a_commits.push_back(
            MerkleTree::commit(std::vector<FieldElement>{fe(i), fe(i*2), fe(i*3), fe(i*4)}));

    auto tree_phi = MerkleTree::tree_commits(c_a_commits);
    Commitment c_phi = tree_phi.root();

    // La racine doit etre la meme que le raccourci commit()
    CHECK(c_phi == MerkleTree::commit(c_a_commits));

    // Chaque commitment c_a doit etre verifiable dans c_phi
    for (size_t i = 0; i < c_a_commits.size(); ++i) {
        CHECK(verify_path(c_phi, c_a_commits[i], tree_phi.path(i)));
    }
}

TEST_CASE("tree_commits() : liste vide") {
    CHECK_THROWS_AS(
        MerkleTree::tree_commits(std::vector<Commitment>{}),
        std::invalid_argument);
}

// MerkleTree::commit() raccourcis

TEST_CASE("commit(FieldElements) == tree_field().root()") {
    init_field(NTL::ZZ(101));
    std::vector<FieldElement> vals = {fe(3), fe(7), fe(11), fe(13)};
    CHECK(MerkleTree::commit(vals) ==
          MerkleTree::tree_field(vals).root());
}

TEST_CASE("commit(Commitments) == tree_commits().root()") {
    std::vector<Commitment> comms = {
        fill_digest(0x01), fill_digest(0x02), fill_digest(0x03)
    };
    CHECK(MerkleTree::commit(comms) ==
          MerkleTree::tree_commits(comms).root());
}