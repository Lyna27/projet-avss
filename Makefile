CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Iinclude -Itests -isystem /opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lntl -lgmp

test_field: src/field.cpp tests/test_field.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp tests/test_field.cpp -o test_field $(LDFLAGS)

test_polynomial: src/field.cpp src/polynomial.cpp tests/test_polynomial.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/polynomial.cpp tests/test_polynomial.cpp -o test_polynomial $(LDFLAGS)

test_merkle: src/field.cpp src/merkle.cpp tests/test_merkle.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/merkle.cpp tests/test_merkle.cpp -o test_merkle $(LDFLAGS)

test_proximity: src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp tests/test_proximity.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp tests/test_proximity.cpp -o test_proximity $(LDFLAGS)

test_deal: src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp tests/test_deal.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp tests/test_deal.cpp -o test_deal $(LDFLAGS)

test_getshare: src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp src/getshare.cpp tests/test_getshare.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp src/getshare.cpp tests/test_getshare.cpp -o test_getshare $(LDFLAGS)

avss_demo: src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp src/getshare.cpp main.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp src/getshare.cpp main.cpp -o avss_demo $(LDFLAGS)

avss_bench: src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp src/getshare.cpp benchmarks/bench.cpp
	$(CXX) $(CXXFLAGS) src/field.cpp src/polynomial.cpp src/merkle.cpp src/proximity_proof.cpp src/deal.cpp src/getshare.cpp benchmarks/bench.cpp -o avss_bench $(LDFLAGS) -I/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/Headers -I/Users/apple/Library/Python/3.9/lib/python/site-packages/numpy/_core/include -F/Library/Developer/CommandLineTools/Library/Frameworks -framework Python3 -Wno-deprecated-declarations -DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION -Wl,-rpath,/Library/Developer/CommandLineTools/Library/Frameworks

clean:
	rm -f test_field test_ntl test_polynomial test_merkle test_proximity test_deal test_getshare avss_demo avss_bench 