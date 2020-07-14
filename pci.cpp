#include <iostream>
#include <fstream>
#include <map>
#include <complex>
#include <utility>
#include <Eigen>

using namespace std;

// for S_3
enum group { e, p12, p23, p13, p123, p132 };

// declaring multiplication in S_3 (very ugly rn)
inline group operator* (group& g1, group& g2) {
    if (g1 == group::e) {
        return g2;
    }
    else if (g2 == group::e) {
        return g1;
    }
    else if (g1 == g2 && g1 < group::p123) {
        return group::e;
    }
    else if (g1 == g2 && g1 == group::p123) {
        return group::p132;
    }
    else if (g1 == g2 && g1 == group::p132) {
        return group::p123;
    }
    else if (g1 == group::p12) {
        switch (g2) {
            case group::p23:
                return group::p123;
            case group::p13:
                return group::p132;
            case group::p123:
                return group::p23;
            case group::p132:
                return group::p13;
        }
    }
    else if (g1 == group::p23) {
        switch (g2) {
            case group::p12:
                return group::p132;
            case group::p13:
                return group::p123;
            case group::p123:
                return group::p13;
            case group::p132:
                return group::p12;
        }
    }
    else if (g1 == group::p13) {
        switch (g2) {
            case group::p12:
                return group::p123;
            case group::p23:
                return group::p132;
            case group::p123:
                return group::p12;
            case group::p132:
                return group::p23;
        }
    }
    else if (g1 == group::p123) {
        switch (g2) {
            case group::p12:
                return group::p13;
            case group::p23:
                return group::p12;
            case group::p13:
                return group::p23;
            case group::p132:
                return group::e;
        }
    }
    else if (g1 == group::p132) {
        switch (g2) {
            case group::p12:
                return group::p23;
            case group::p23:
                return group::p13;
            case group::p13:
                return group::p12;
            case group::p123:
                return group::e;
        }
    }
}

group invert(group g) {
    switch (g) {
        case group::e :
            return group::e;
        case group::p12 :
            return group::p12;
        case group::p23 :
            return group::p23;
        case group::p13 :
            return group::p13;
        case group::p123 :
            return group::p132;
        case group::p132 :
            return group::p123;

    }
}

int main() {
    // enumerate \{ V_{g,h} \} as V_j
    pair<group, group> V[36];
    pair<group, group> p;
    for (size_t i = 0; i < 36; ++i) {
        p.first = (group) (i/6);
        p.second = (group) (i % 6);
        V[i] = p;
    }

    // fill Z(A(G)) = Ker(Z)
    Eigen::Matrix<complex<double>, 1296, 36> Z;

    int x, y, z;
    
    group a, b, c, d;

    // d_1 = d_{xy}^z, d_2 = d_{yx}^z
    // d_{xy}^z = delta_{bab^-1, c} if V_z = V_{a, db} where V_x = V_a,b, V_y = V_c,d
    int d_1 = 0, d_2 = 0;

    // 1296 = 36^4 = max((i-1)36 + k), 1 \leq i,k \leq 36
    for (size_t i = 0; i < 1296; ++i) {
        x = i / 36 + 1;
        z = i % 36; 
        for (size_t j = 0; j < 36; ++j) {
            d_1 = 0;
            d_2 = 0;
            y = j;
            a = V[x].first;
            b = V[x].second;
            c = V[y].first;
            d = V[y].second;
            if (V[z].first == a && V[z].second == d*b) {
                d_1 = 1;
            }
            else if (V[z].first == c && V[z].second == b*d) {
                d_2 = 1;
            }
            Z(i,j) = d_1 - d_2;
        }
    }
    
    // calculate the basis of nullspace of Z
    Eigen::FullPivLU<Eigen::Matrix<complex<double>,1296, 36>> lu(Z);
    Eigen::Matrix<complex<double>, 1296, 36> Z_null_space = lu.kernel();

}