#include <iostream>
#include <NTL/ZZ.h>
using namespace std;
using namespace NTL;

int main() {
    ZZ a, b, c;
    a = 123456789;
    b = 987654321;
    c = a * b;
    cout << "Le produit est : " << c << endl;
    return 0;
}