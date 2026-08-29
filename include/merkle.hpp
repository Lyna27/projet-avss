#ifndef MERKLE_H
#define MERKLE_H

#include <array>
#include <vector>
#include <cstdint>
#include <cstddef>
#include "field.hpp"

constexpr size_t HASH_SIZE = 32; // 32 octets = 256 bits pour SHA-256
using Hash_type = std::array<uint8_t, HASH_SIZE>;
using Commitment = Hash_type;

/*
 * Application de hachage standard H : {0,1}^* -> {0,1}^256.
 * Prend en entree un message de taille arbitraire data et sa taille len
 * et renvoie le hash SHA-256 sur 32 octets.
 */
Hash_type hash_bytes(const uint8_t* data, size_t len);

/* 
 * Hache un element de F_p en le convertissant en une chaine 
 * d'octets de taille fixe avant d'appliquer SHA-256.
 */
Hash_type hash_element(const FieldElement& x);

/*
 * Fonction de concatenation de deux hashs de 32 octets 
 * chacun pour obtenir un nouveau hash.
 */ 
Hash_type hash_concat(const Hash_type& left, const Hash_type& right);

// Structure d'un chemin dans un arbre de Merkle
struct MerklePath {
    std::vector<Hash_type> siblings;
    // Comme Hash(A,B) != Hash(B,A) : 
    std::vector<bool> directions; // true = sibling a droite, false sinon
};

namespace MerkleTree {

    // Structure d'un arbre de Merkle complet avec sauvegarde de tous les niveaux
    struct Tree {
        // Les différents niveaux : levels[0] = feuilles, ..., levels.back() = {root}
        std::vector<std::vector<Hash_type>> levels; 

        // Fonction qui renvoie la racine de l'arbre a partir des feuilles
        Commitment root() const; 

        /*
         * Fonction qui prend en entree la position d'une feuille et 
         * construit son chemin vers la racine.
         * C'est la preuve cryptographique qui permet de verifier 
         * l'appartenance de la feuille a l'arbre sans avoir a connaitre 
         * le reste des donnees.
         */
        MerklePath path(size_t position) const; 
    };

    /* Fonction qui prend en entrees des blocs de 32 octets (Hash_type)
     * et construit l'arbre complet a partir de feuilles deja hachees.
     */ 
    Tree build(const std::vector<Hash_type>& leaves);

    /*
     * Genere un arbre de Merkle a partir d'elements dans F_p
     * Appelle hash_element() pour le hachage puis build() 
     * pour la construction.
     */ 

    Tree tree_field(const std::vector<FieldElement>& values);

    /* 
     * Cree un arbre de Merkle a partir d'engagements polynomiaux.
     * Comme Commitment = Hash_type, il suffit d'appeler build().
     */
    Tree tree_commits(const std::vector<Commitment>& values);

    /* Calcule la racine de l'arbre Merkle a partir 
     * d'une liste d'elements de F_p sans enregistrer 
     les niveaux de l'arbre complet 
     */ 
    Commitment commit(const std::vector<FieldElement>& values);

    /* Calcule la racine de l'arbre Merkle a partir 
     * d'une liste de commitments sans enregistrer 
     les niveaux de l'arbre complet 
     */ 
    Commitment commit(const std::vector<Commitment>&   values);

} // fin de namespace MerkleTree

/* Fonction qui verifie qu'une feuille appartient bien a
 * l'arbre en recalculant la racine a partir du chemin 
 * MerklePath fourni.
 * Renvoie true si la racine recalculee correspond a la 
 * racine donnee, false si la preuve ou la feuille est invalide
 */ 
bool verify_path(const Commitment& root, const Hash_type& leaf, const MerklePath& path);

#endif // MERKLE_H