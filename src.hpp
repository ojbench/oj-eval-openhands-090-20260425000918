#ifndef SRC_HPP
#define SRC_HPP

#include <bits/stdc++.h>
#include "fraction.hpp"

class matrix {
private:
    int m, n;
    fraction **data;

    void allocate(int rows, int cols) {
        m = rows; n = cols;
        if (m <= 0 || n <= 0) { data = nullptr; return; }
        data = new fraction*[m];
        for (int i = 0; i < m; ++i) {
            data[i] = new fraction[n];
        }
    }

    void deallocate() {
        if (!data) return;
        for (int i = 0; i < m; ++i) delete[] data[i];
        delete[] data;
        data = nullptr; m = n = 0;
    }

public:
    matrix() { m = n = 0; data = nullptr; }

    matrix(int m_, int n_) { data = nullptr; allocate(m_, n_); }

    matrix(const matrix &obj) {
        data = nullptr; allocate(obj.m, obj.n);
        if (data) {
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    data[i][j] = obj.data[i][j];
        }
    }

    matrix(matrix &&obj) noexcept {
        m = obj.m; n = obj.n; data = obj.data;
        obj.m = obj.n = 0; obj.data = nullptr;
    }

    ~matrix() { deallocate(); }

    matrix &operator=(const matrix &obj) {
        if (this == &obj) return *this;
        if (m != obj.m || n != obj.n) {
            deallocate();
            allocate(obj.m, obj.n);
        }
        if (data) {
            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    data[i][j] = obj.data[i][j];
        }
        return *this;
    }

    fraction &operator()(int i, int j) {
        if (i < 1 || i > m || j < 0 || j >= n || data == nullptr) throw matrix_error();
        return data[i-1][j];
    }

    friend matrix operator*(const matrix &lhs, const matrix &rhs) {
        if (lhs.n == 0 || rhs.m == 0 || lhs.n != rhs.m || lhs.data == nullptr || rhs.data == nullptr)
            throw matrix_error();
        matrix res(lhs.m, rhs.n);
        for (int i = 0; i < lhs.m; ++i) {
            for (int k = 0; k < lhs.n; ++k) {
                const fraction &lik = lhs.data[i][k];
                if (lik == fraction(0)) continue;
                for (int j = 0; j < rhs.n; ++j) {
                    res.data[i][j] = res.data[i][j] + lik * rhs.data[k][j];
                }
            }
        }
        return res;
    }

    matrix transposition() {
        if (data == nullptr || m == 0 || n == 0) throw matrix_error();
        matrix t(n, m);
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                t.data[j][i] = data[i][j];
        return t;
    }

    fraction determination() {
        if (data == nullptr || m == 0 || n == 0 || m != n) throw matrix_error();
        std::vector<std::vector<fraction>> a(m, std::vector<fraction>(n));
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                a[i][j] = data[i][j];
        int sign = 1;
        int row = 0;
        for (int col = 0; col < n && row < n; ++col) {
            int piv = row;
            while (piv < n && a[piv][col] == fraction(0)) ++piv;
            if (piv == n) continue; // no pivot in this column
            if (piv != row) { std::swap(a[piv], a[row]); sign = -sign; }
            fraction pivval = a[row][col];
            for (int r = row + 1; r < n; ++r) if (!(a[r][col] == fraction(0))) {
                fraction factor = a[r][col] / pivval;
                for (int c = col; c < n; ++c) a[r][c] = a[r][c] - factor * a[row][c];
            }
            ++row;
        }
        if (row < n) return fraction(0);
        fraction det(1);
        for (int i = 0; i < n; ++i) det = det * a[i][i];
        if (sign == -1) det = fraction(0) - det;
        return det;
    }

    int rows() const { return m; }
    int cols() const { return n; }
};

class resistive_network {
private:
    int interface_size, connection_size;
    matrix adjacency, conduction; // A (n x m), C (m x m diag)

    std::vector<int> from_v, to_v;
    std::vector<fraction> r_v, g_v; // resistance and conductance

    matrix laplacian;     // L = A C A^T (n x n)
    matrix laplacian_red; // L' remove last row&col ( (n-1) x (n-1) )

    static std::vector<fraction> solve_linear(matrix A, const std::vector<fraction> &b) {
        int n = A.rows();
        if (n != A.cols() || (int)b.size() != n) throw matrix_error();
        std::vector<std::vector<fraction>> a(n, std::vector<fraction>(n + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) a[i][j] = A(i+1, j);
            a[i][n] = b[i];
        }
        for (int col = 0; col < n; ++col) {
            int piv = col;
            while (piv < n && a[piv][col] == fraction(0)) ++piv;
            if (piv == n) throw matrix_error();
            if (piv != col) std::swap(a[piv], a[col]);
            fraction pivval = a[col][col];
            for (int r = col + 1; r < n; ++r) if (!(a[r][col] == fraction(0))) {
                fraction factor = a[r][col] / pivval;
                for (int c = col; c <= n; ++c) a[r][c] = a[r][c] - factor * a[col][c];
            }
        }
        std::vector<fraction> x(n, fraction(0));
        for (int i = n - 1; i >= 0; --i) {
            fraction sum(0);
            for (int j = i + 1; j < n; ++j) sum = sum + a[i][j] * x[j];
            fraction rhs = a[i][n] - sum;
            x[i] = rhs / a[i][i];
        }
        return x;
    }

    static matrix build_laplacian(const matrix &A, const matrix &C) {
        // L = A C A^T
        matrix AT = matrix(A).transposition();
        return A * (C * AT);
    }

public:
    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) :
        interface_size(interface_size_), connection_size(connection_size_),
        adjacency(interface_size_, connection_size_), conduction(connection_size_, connection_size_) {
        from_v.resize(connection_size);
        to_v.resize(connection_size);
        r_v.resize(connection_size);
        g_v.resize(connection_size);
        for (int k = 0; k < connection_size; ++k) {
            from_v[k] = from[k];
            to_v[k] = to[k];
            r_v[k] = resistance[k];
            g_v[k] = fraction(1) / resistance[k];
        }
        // Build incidence matrix A (n x m): +1 at from-1, -1 at to-1
        for (int i = 1; i <= interface_size; ++i)
            for (int j = 0; j < connection_size; ++j)
                adjacency(i, j) = fraction(0);
        for (int k = 0; k < connection_size; ++k) {
            int u = from_v[k] - 1;
            int v = to_v[k] - 1;
            adjacency(u + 1, k) = adjacency(u + 1, k) + fraction(1);
            adjacency(v + 1, k) = adjacency(v + 1, k) - fraction(1);
        }
        // Build conduction diagonal matrix C (m x m)
        for (int i = 1; i <= connection_size; ++i)
            for (int j = 0; j < connection_size; ++j)
                conduction(i, j) = fraction(0);
        for (int k = 0; k < connection_size; ++k) conduction(k + 1, k) = g_v[k];

        laplacian = build_laplacian(adjacency, conduction); // n x n
        // Build reduced Laplacian removing last node
        laplacian_red = matrix(interface_size - 1, interface_size - 1);
        for (int i = 0; i < interface_size - 1; ++i)
            for (int j = 0; j < interface_size - 1; ++j)
                laplacian_red(i + 1, j) = laplacian(i + 1, j);
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        if (interface_id1 == interface_id2) return fraction(0);
        // Inject 1A at id1, -1A at id2, ground node n (voltage at node n is 0)
        std::vector<fraction> I(interface_size, fraction(0));
        I[interface_id1 - 1] = I[interface_id1 - 1] + fraction(1);
        I[interface_id2 - 1] = I[interface_id2 - 1] - fraction(1);
        // Reduced current vector excludes ground (node n)
        std::vector<fraction> Ired(interface_size - 1);
        for (int i = 0; i < interface_size - 1; ++i) Ired[i] = I[i];
        std::vector<fraction> Ured = solve_linear(laplacian_red, Ired);
        fraction ui = (interface_id1 == interface_size) ? fraction(0) : Ured[interface_id1 - 1];
        fraction uj = (interface_id2 == interface_size) ? fraction(0) : Ured[interface_id2 - 1];
        fraction req = ui - uj; // since current is 1A
        if (req == fraction(0)) return fraction(0);
        return req;
    }

    fraction get_voltage(int id, fraction current[]) {
        // Build I, with u_n = 0
        std::vector<fraction> Ired(interface_size - 1);
        for (int i = 0; i < interface_size - 1; ++i) Ired[i] = current[i];
        std::vector<fraction> Ured = solve_linear(laplacian_red, Ired);
        return Ured[id - 1];
    }

    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (int k = 0; k < connection_size; ++k) {
            int u = from_v[k] - 1;
            int v = to_v[k] - 1;
            fraction du = voltage[u] - voltage[v];
            total = total + g_v[k] * du * du;
        }
        return total;
    }
};


#endif //SRC_HPP
