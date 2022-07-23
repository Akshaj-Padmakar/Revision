/*
CONTENTS -  *> PRAGMAS STUFF
            *> MOD TEMPLATE/ MATHS TEMPLATE
            *> LONGEST SUBSEQUNECES STUFF
            *> COORDINATE COMPRESSION
            *> GEOMETRY POINT STRUCT
            *> SEGTREE STUFF (EX - LESS THAN ON LEFT )
            *> MATRIX EXPO
            *> LONG-DOUBLE-MARTIX EXPO  !
            *> FAST PRIME CHECK
*/
#pragma gcc optimise "trapv"
//for checking overflow -> gives rte rather than WA


fast pragma
#pragma GCC optimize("-funsafe-loop-optimizations")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-fwhole-program")
#pragma GCC optimize("Ofast,no-stack-protector")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fcse-skip-blocks")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-fstrict-overflow")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")





// -------------------------------<MOD TEMPLATE>-----------------------------

const int MOD = 998244353;

template<typename T> T gcd(T a, T b) {if (b == 0) return a; a %= b; return gcd(b, a);}
template<typename T> T lcm(T a, T b) {return (a * (b / gcd(a, b)));}
int add(int a, int b, int c = MOD) {int res = a + b; return (res >= c ? res % c : res);}
int sub(int a, int b, int c = MOD) {int res; if (abs(a - b) < c)res = a - b; else res = (a - b) % c; return (res < 0 ? res + c : res);}
int mul(int a, int b, int c = MOD) {int res = (long long)a * b; return (res >= c ? res % c : res);}

template<typename T>T binpow(T e, T n) {T x = 1, p = e; while (n) {if (n & 1)x = x * p; p = p * p; n >>= 1;} return x;}
template<typename T>T binpow2(T e, T n, T m = MOD) {T x = 1, p = e; while (n) {if (n & 1)x = mul(x, p, m); p = mul(p, p, m); n >>= 1;} return x;}
template<typename T>T extended_euclid(T a, T b, T &x, T &y) {
    T xx = 0, yy = 1; y = 0; x = 1; while (b) {
        T q = a / b, t = b; b = a % b; a = t; \
        t = xx; xx = x - q * xx; x = t; t = yy; yy = y - q * yy; y = t;
    } return a;
}
template<typename T>T mod_inverse(T a, T n = MOD) {T x, y, z = 0; T d = extended_euclid(a, n, x, y); return (d > 1 ? -1 : sub(x, z, n));}




//-------------<for dealing with -ve numbers but this is slower>---------

const int MOD = 998244353;

template<typename T> T gcd(T a, T b) {if (b == 0) return a; a %= b; return gcd(b, a);}
template<typename T> T lcm(T a, T b) {return (a * (b / gcd(a, b)));}

void checker(int &res) {res %= MOD; res += MOD; res %= MOD;}
int add(int a, int b, int c = MOD) {int res = a + b; checker(res); return res;}
int sub(int a, int b, int c = MOD) {int res = a - b; checker(res); return res;}
int mul(int a, int b, int c = MOD) {int res = (long long) a * b; checker(res); return res;}

template<typename T>T binpow(T e, T n) {T x = 1, p = e; while (n) {if (n & 1)x = x * p; p = p * p; n >>= 1;} return x;}
template<typename T>T binpow2(T e, T n, T m = MOD) {T x = 1, p = e; while (n) {if (n & 1)x = mul(x, p, m); p = mul(p, p, m); n >>= 1;} return x;}
template<typename T>T extended_euclid(T a, T b, T &x, T &y) {
    T xx = 0, yy = 1; y = 0; x = 1; while (b) {
        T q = a / b, t = b; b = a % b; a = t; \
        t = xx; xx = x - q * xx; x = t; t = yy; yy = y - q * yy; y = t;
    } return a;
}
template<typename T>T mod_inverse(T a, T n = MOD) {T x, y, z = 0; T d = extended_euclid(a, n, x, y); return (d > 1 ? -1 : sub(x, z, n));}





string to_binary(int n) {string r; while (n != 0) {r = (n % 2 == 0 ? "0" : "1") + r; n /= 2;} return r;}
int frm_binary(string s) {int val = 0, cnt = 0; while (s.size()) {int x = s.back() - '0'; s.pop_back(); val += x * (1ll << cnt); cnt++;} return val;}



// -------------------------------<MATH COMPLETE>-----------------------------

const int MOD = 998244353;
const int MAXN = 2e5 + 5;

template<typename T> T gcd(T a, T b) {if (b == 0) return a; a %= b; return gcd(b, a);}
template<typename T> T lcm(T a, T b) {return (a * (b / gcd(a, b)));}

int add(int a, int b, int c = MOD) {int res = a + b; return (res >= c ? res % c : res);}
int sub(int a, int b, int c = MOD) {int res; if (abs(a - b) < c)res = a - b; else res = (a - b) % c; return (res < 0 ? res + c : res);}
int mul(int a, int b, int c = MOD) {int res = (long long)a * b; return (res >= c ? res % c : res);}

template<typename T>T binpow(T e, T n) {T x = 1, p = e; while (n) {if (n & 1)x = x * p; p = p * p; n >>= 1;} return x;}
template<typename T>T binpow2(T e, T n, T m = MOD) {T x = 1, p = e; while (n) {if (n & 1)x = mul(x, p, m); p = mul(p, p, m); n >>= 1;} return x;}
template<typename T>T extended_euclid(T a, T b, T &x, T &y) {
    T xx = 0, yy = 1; y = 0; x = 1; while (b) {
        T q = a / b, t = b; b = a % b; a = t; \
        t = xx; xx = x - q * xx; x = t; t = yy; yy = y - q * yy; y = t;
    } return a;
}
template<typename T>T mod_inverse(T a, T n = MOD) {T x, y, z = 0; T d = extended_euclid(a, n, x, y); return (d > 1 ? -1 : sub(x, z, n));}


const int FACSZ = 5005; // Make sure to change this


vector<int> fact(FACSZ), ifact(FACSZ);
// factorial and inverse factorial
// till FACSZ

void precom(int c = MOD) {
    fact[0] = 1;
    for (int i = 1; i < FACSZ; i++) {
        fact[i] = mul(i, fact[i - 1], c);
    }

    ifact[FACSZ - 1] = mod_inverse(fact[FACSZ - 1], c);
    for (int i = FACSZ - 2; i >= 0; i--) {
        ifact[i] = mul(i + 1, ifact[i + 1], c);
    }
}

vector<int> primes;
void prime_precom() {
    /*
    size     last prime
    78499 -> ~1e6
    1e5   -> ~1299743 (1.2s)
    1e5   -> ~15485917 (8s)
    */
    primes.push_back(2);
    for (int x = 3; primes.size() <= MAXN; x += 2) {
        bool isPrime = true;
        for (auto p : primes) {
            if (x % p == 0) {
                isPrime = false; break;
            }
            if (p * p > x) {
                break;
            }
        }
        if (isPrime) {
            primes.push_back(x);
        }
    }
}
int ncr(int n, int k) {
    //O(k)
    if (n < k)
        return 0;
    if (k == 0)
        return 1;
    int res = 1; if (k > n - k)
        k = n - k;
    for (int i = 0; i < k; ++i) {
        res *= (n - i); res /= (i + 1);
    }
    return res;
}
int ncr_modp(int n, int k, int c = MOD) {
    //O(k)
    if (n < k) {
        return 0;
    }
    if (k == 0) {
        return 1;
    }

    int res = 1;

    if (k > n - k) {
        k = n - k;
    }
    for (int i = 0; i < k; ++i) {
        res = mul(res, n - i, c);
        res = mul(res, binpow2(i + 1, c - 2, c), c);
    }
    return res;
}


void factorize(int a, vector<int> &factors) {
    for (int i = 1; i * i <= a; i++) {
        if (a % i == 0) {
            factors.push_back(i);
            factors.push_back(a / i);
        }
    }
    sort(factors.begin(), factors.end());
}

int ncr_precom(int n, int r, int c = MOD) {
    //define fact and ifact first
    return mul(mul(ifact[r], ifact[n - r], c), fact[n], c);
}
int ceil(int a, int b) {
    return (a + b - 1) / b;
}
bool is_prime(int n) {
    for (int i = 2; i * i <= n; i++) {
        if ( n % i == 0 ) {
            return 0;
        }
    }
    return 1;
}

//could be MLE
vector<int> spf;
void sieve() {
    spf.resize(MAXN);
    spf[1] = 1;
    for (int i = 2; i < MAXN; i++) {
        spf[i] = i;
    }
    for (int i = 4; i < MAXN; i += 2) {
        spf[i] = 2;
    }
    for (int i = 3; i * i <= MAXN; i += 2) {
        if (spf[i] == i) {
            for (int j = i * i; j < MAXN; j += i) {
                if (spf[j] == j) {
                    spf[j] = i;
                }
            }
        }
    }
}
//Linear Sieve !
void sieve() {
    spf.resize(MAXN);
    spf[1] = 1;
    for (int i = 2; i < MAXN; i += 2) {
        spf[i] = 2;
    }

    for (int i = 3; i < MAXN; i += 2) {
        spf[i] = i;
    }
    vector<int> pp;
    for (int i = 3; i < MAXN; i += 2) {
        if (spf[i] == i) {
            pp.push_back(i);
        }
        for (int j = 0; j < pp.size() && spf[j] * i < MAXN; j++) {
            if (i % pp[j] == 0) {
                //not a spf now !
                break;
            }
            spf[i * pp[j]] = pp[j];
        }
    }
}


void pfactor(int x, vector<int> &pfac) {
    while (x != 1) {
        pfac.push_back(spf[x]);
        int z = spf[x];
        while (x % z == 0) {
            x /= z;
        }

    }
}

// -------------------------------</MATH COMPLETE>-----------------------------
// -------------------------------</MOD TEMPLATE>-----------------------------


// -------------------------------<LONGEST_SUBSEQUNECE>-----------------------------
int l_non_ds(vector<int> &a) {
    // longest non-decreasing subsequence
    int n = a.size();
    vector<int> d(n + 1, INF);
    // d[i] -> stores the min last element for i-length subsequence
    d[0] = -INF;
    for (int i = 0; i < n; i++) {
        int j = upper_bound(d.begin(), d.end(), a[i]) - d.begin();
        if (d[j - 1] <= a[i] && a[i] <= d[j]) {
            d[j] = a[i];
        }
    }
    int ans = 0;
    for (int i = 0; i <= n; i++) {
        if (d[i] < INF) {
            ans = i;
        }
    }
    return ans;
}
int lis(vector<int> &a) {
    // longest increasing subsequence
    int n = a.size();
    vector<int> d(n + 1, INF);
    // d[i] -> stores the min last element  for i-length subsequence
    d[0] = -INF;
    for (int i = 0; i < n; i++) {
        int j = upper_bound(d.begin(), d.end(), a[i]) - d.begin();
        if (d[j - 1] < a[i] && a[i] < d[j]) {
            d[j] = a[i];
        }
    }
    int ans = 0;
    for (int i = 0; i <= n; i++) {
        if (d[i] < INF) {
            ans = i;
        }
    }
    return ans;
}

int lds(vector<int> &a) {
    // longest decreasing subsequence
    int n = a.size();

    reverse(a.begin(), a.end());

    vector<int> d(n + 1, INF);
    // d[i] -> stores the min last element for i-length subsequence
    d[0] = -INF;
    for (int i = 0; i < n; i++) {
        int j = upper_bound(d.begin(), d.end(), a[i]) - d.begin();
        if (d[j - 1] < a[i] && a[i] < d[j]) {
            d[j] = a[i];
        }
    }
    int ans = 0;
    for (int i = 0; i <= n; i++) {
        if (d[i] < INF) {
            ans = i;
        }
    }
    return ans;
}

int l_non_is(vector<int> &a) {
    // longest non-increasing subsequence
    int n = a.size();

    reverse(a.begin(), a.end());

    vector<int> d(n + 1, INF);
    // d[i] -> stores the min last element for i-length subsequence
    d[0] = -INF;
    for (int i = 0; i < n; i++) {
        int j = upper_bound(d.begin(), d.end(), a[i]) - d.begin();
        if (d[j - 1] <= a[i] && a[i] <= d[j]) {
            d[j] = a[i];
        }
    }
    int ans = 0;
    for (int i = 0; i <= n; i++) {
        if (d[i] < INF) {
            ans = i;
        }
    }
    return ans;
}
// -------------------------------</LONGEST_SUBSEQUNECE>-----------------------------




// -------------------------------</COORDINATE COMPRESS>-----------------------------

auto coo_compress = [&](vector<int> &x) {
    // coordinates compresse x and over-writes on it !
    int sz = x.size(), id = 0;
    map<int, int> mp;

    for (int i = 0; i < sz; i++) {
        mp[x[i]] = 1;
    }

    for (auto x : mp) {
        mp[x.first] = id++;
    }

    for (int i = 0; i < sz; i++) {
        x[i] = mp[x[i]];
    }
};
// -------------------------------</COORDINATE COMPRESS>-----------------------------


// ---------------------------------<GEOMETRY>---------------------------------------
struct P {
    int x, y;
    void read() {
        cin >> x >> y;
    }

    P operator -=(P b) {
        return {x -= b.x, y -= b.y};
    }

    int operator *(P b) {
        return x * b.y - y * b.x;
    }
};

LD angle(P &p) {
    if (p.x == 0) {
        if (p.y == 0) {
            return INF;
        } else if (p.y > 0) {
            return PI / 2.0;
        } else {
            return -PI / 2.0;
        }
    }

    if (p.x > 0) {
        if (p.y > 0) {
            return atan2((LD)p.y, (LD)p.x);
        } else if (p.y == 0) {
            return 0;
        } else {
            return -atan2((LD)(-1 * p.y), (LD)p.x);
        }
    } else {
        if (p.y > 0) {
            return PI - atan2((LD)abs(p.y), (LD)abs(p.x));
        } else if (p.y == 0) {
            return -PI;
        } else {
            return -PI + atan2((LD)abs(p.y), (LD)abs(p.x));
        }
    }
}


// ---------------------------------</GEOMETRY>---------------------------------------





// -------------------------------<SEGTREE_ADDITIONAL_STUFF>-----------------------------


auto greater_than_equal_right = [&](int x, int id, int st_sz) {
    /*
    Instructions - *> Create MAX segtree on the vector[v] you want to search
                   *> returns if there exist a i such that i >= id ans v[i] >= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */

    target = x;
    auto i = st.max_right(id, f); // 'st' used !

    return (i == st_sz ? -1ll : i);
};




auto less_than_equal_right = [&](int x, int id, int st_sz) {
    /*
    Instructions - *> Create MAX segtree on the negative of the vector[v] you want to search
                   *> returns if there exist a i such that i >= id ans v[i] <= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */
    target = -x;
    auto i = st_neg.max_right(id, f); // 'st_neg' used !

    return (i == st_sz ? -1ll : i);
};



auto greater_than_equal_left = [&](int x, int id, int st_sz) {
    /*
    Instructions - *> Create MAX segtree on the reverse of the vector[v] you want to search
                   *> returns if there exist a i such that i <= id ans v[i] >= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */
    target = x;
    id = st_sz - id - 1;
    auto i = st_rev.max_right(id, f);// 'st_rev' used !

    return (i == st_sz ? -1ll : st_sz - 1 - i); // st_sz == n  -> size of the array??
};



auto less_than_equal_left = [&](int x, int id, int st_sz) {
    /*
    Instructions - *> Create MAX segtree on the reverse of negative of the vector[v] you want to search
                   *> returns if there exist a i such that i <= id ans v[i] <= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */
    target = -x;
    id = st_sz - id - 1;
    auto i = st_rev_neg.max_right(id, f);// 'st_rev_neg' used !
    return (i == st_sz ? -1ll : st_sz - 1 - i); // st_sz == n  -> size of the array ??
};



//--------------------------------





auto greater_than_equal_right = [&](int x, int id, SegmentTree<int> &segV) {
    /*
    Instructions - *> Create MAX segtree on the vector[v] you want to search
                   *> returns if there exist a i such that i >= id ans v[i] >= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */

    target = x;

    auto check1 = [&](int x) {
        return x >= target;
    };

    return segV.find_first(id, check1);
};




auto less_than_equal_right = [&](int x, int id, SegmentTree<int> &segV) {
    /*
    Instructions - *> Create MIN segtree on the vector[v] you want to search
                   *> returns if there exist a i such that i >= id ans v[i] <= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */
    target = x;

    auto check1 = [&](int x) {
        return x <= target;
    };
    return segV.find_first(id, check1);
};



auto greater_than_equal_left = [&](int x, int id, SegmentTree<int> &segV) {
    /*
    Instructions - *> Create MAX segtree on the reverse of the vector[v] you want to search
                   *> returns if there exist a i such that i <= id ans v[i] >= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */
    target = x;

    auto check1 = [&](int x) {
        return x >= target;
    };
    return segV.find_last(id, check1);
};



auto less_than_equal_left = [&](int x, int id, SegmentTree<int> &segV) {
    /*
    Instructions - *> Create MIN segtree on the reverse of negative of the vector[v] you want to search
                   *> returns if there exist a i such that i <= id ans v[i] <= x, O.W returns -1
                   *> CONDITION -> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    */
    target = x;
    auto check1 = [&](int x) {
        return x <= target;
    };
    return segV.find_last(id, check1);
};



// -------------------------------<SEGTREE_ADDITIONAL_STUFF>-----------------------------





// ----------------------------<MATRIX>-------------------------------------

const long long MOD2 = static_cast<long long>(MOD) * MOD;

struct Matrix {
    vector< vector<int> > mat;
    int n_rows, n_cols;

    Matrix() {}

    Matrix(vector< vector<int> > values): mat(values), n_rows(values.size()),
        n_cols(values[0].size()) {}

    static Matrix identity_matrix(int n) {
        vector< vector<int> > values(n, vector<int>(n, 0));
        for (int i = 0; i < n; i++)
            values[i][i] = 1;
        return values;
    }

    Matrix operator*(const Matrix &other) const {
        int n = n_rows, m = other.n_cols;
        vector< vector<int> > result(n_rows, vector<int>(n_cols, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                long long tmp = 0;
                for (int k = 0; k < n_cols; k++) {
                    tmp += mat[i][k] * 1ll * other.mat[k][j];
                    while (tmp >= MOD2)
                        tmp -= MOD2;
                }
                result[i][j] = tmp % MOD;
            }

        return move(Matrix(move(result)));
    }

    inline bool is_square() const {
        return n_rows == n_cols;
    }
};

// binary exponentiation, returns a^p
Matrix pw(Matrix a, int p) {
    Matrix result = Matrix::identity_matrix(a.n_cols);
    while (p > 0) {
        if (p & 1)
            result = a * result;
        a = a * a;
        p >>= 1;
    }
    return result;
}



// ----------------------------</MATRIX>-------------------------------------





// ----------------------------<LONG DOUBLE MATRIX>-------------------------------------


struct Matrix {
    vector< vector<long double> > mat;
    int n_rows, n_cols;

    Matrix() {}

    Matrix(vector< vector<long double> > values): mat(values), n_rows(values.size()),
        n_cols(values[0].size()) {}

    static Matrix identity_matrix(int n) {
        vector< vector< long double> > values(n, vector<long double>(n, 0.0));
        for (int i = 0; i < n; i++)
            values[i][i] = 1.0;
        return values;
    }

    Matrix operator*(const Matrix &other) const {
        int n = n_rows, m = other.n_cols;
        vector< vector<long double> > result(n_rows, vector<long double>(n_cols, 0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                long double tmp = 0;
                for (int k = 0; k < n_cols; k++) {
                    tmp += mat[i][k] * 1ll * other.mat[k][j];
                    while (tmp >= MOD2)
                        tmp -= MOD2;
                }
                result[i][j] = tmp;
            }

        return move(Matrix(move(result)));
    }

    inline bool is_square() const {
        return n_rows == n_cols;
    }
};

// binary exponentiation, returns a^p
Matrix pw(Matrix a, int p) {
    Matrix result = Matrix::identity_matrix(a.n_cols);
    while (p > 0) {
        if (p & 1)
            result = a * result;
        a = a * a;
        p >>= 1;
    }
    return result;
}


// ----------------------------</LONG DOUBLE MATRIX>-------------------------------------


// ----------------------------<PRIME CHECK>-------------------------------------


bool isprime(int n) {
    if (n < 2) {
        return false;
    }
    for (int x : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        if (n == x) {
            return true;
        }
        bool flag = true;
        int r = 1;
        int t = 1;
        while (r <= ((n - 1) >> __builtin_ctzll(n - 1))) {
            if (r & ((n - 1) >> __builtin_ctzll(n - 1))) {
                t = ((__int128)t * x) % n;
            }
            x = ((__int128)x * x) % n;
            r <<= 1;
        }
        if (t == 1 || t == n - 1) {
            flag = false;
        }
        for (r = 0; r < __builtin_ctzll(n - 1); r++) {
            t = ((__int128)t * t) % n;
            if (t == n - 1)
                flag = false;
        }
        if (flag)
            return false;
    }
    return true;
}


// ----------------------------</PRIME CHECK>-------------------------------------


struct mint {
    int vl;
    mint() {
        vl = 0;
    }
    mint(long long v) {
        if (v < 0) {
            v += ((-v + mod - 1) / mod) * mod;
        }
        vl = (int)(v % mod);
    }

    mint operator +(mint oth) {
        return mint((long long)vl + oth.vl);
    }

    mint operator -(mint oth) {
        return mint((long long)vl - oth.vl);
    }

    mint operator *(mint oth) {
        return mint((long long)vl * oth.vl);
    }

    mint powr(long long k) {
        if (k == 0) {
            return mint(1);
        }
        else if (k % 2) {
            return (*this) * powr(k - 1);
        }
        else {
            mint tmp = powr(k / 2);
            return tmp * tmp;
        }
    }

    mint inv() {
        temp = extended_gcd(vl, mod);
        return mint(temp.first);
    }

    mint operator /(mint oth) {
        return (*this) * oth.inv();
    }

    void operator +=(mint oth) {
        (*this) = (*this) + oth;
    }

    void operator -=(mint oth) {
        (*this) = (*this) - oth;
    }

    void operator *=(mint oth) {
        (*this) = (*this) * oth;
    }

    void operator /=(mint oth) {
        (*this) = (*this) / oth;
    }
};
ostream& operator<<(ostream& os, const mint& m)
{
    os << m.vl;
    return os;
}



