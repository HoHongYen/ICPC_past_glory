#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<ll> vll;
typedef complex<double> base;
const double PI = acos(-1);
const int MAXN = 19;
const int maxn = (1 << MAXN);
const int N = 1e5 + 10;
const int MOD = 998244353;
const int BLOCK = 4000;
const ll mod2 = 8ll * MOD * (ll)MOD;
base W[maxn], invW[maxn], P1[maxn], Q1[maxn];
int a[N], fact[N], ifact[N];
vector<ll> b, c, e, temp, Inv;

void precompute_powers()
{
    for (int i = 0; i < maxn / 2; i++)
    {
        double ang = (2 * PI * i) / maxn;
        double _cos = cos(ang), _sin = sin(ang);
        W[i] = base(_cos, _sin);
        invW[i] = base(_cos, -_sin);
    }
}
void fft(vector<base> &a, bool invert)
{
    int n = (int) a.size();

    for (int i = 1, j = 0; i < n; ++i)
    {
        int bit = n >> 1;
        for (; j >= bit; bit >>= 1)
            j -= bit;
        j += bit;
        if (i < j)
            swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1)
    {
        for (int i = 0; i < n; i += len)
        {
            int ind = 0, add = maxn / len;
            for (int j = 0; j < len / 2; ++j)
            {
                base u = a[i + j], v = (a[i + j + len / 2] * (invert ? invW[ind] : W[ind]));
                a[i + j] = (u + v);
                a[i + j + len / 2] = (u - v);
                ind += add;
            }
        }
    }
    if (invert)
        for (int i = 0; i < n; ++i)
            a[i] /= n;
}

void mul_big_mod(vll &a, vll &b)
{
    int n1 = a.size(), n2 = b.size();
    int final_size = a.size() + b.size() - 1;
    int n = 1;
    while (n < final_size)
        n <<= 1;
    vector<base> P(n), Q(n);
    int SQRTMOD = (int)sqrt(MOD) + 10;
    for (int i = 0; i < n1; i++)
        P[i] = base(a[i] % SQRTMOD, a[i] / SQRTMOD);
    for (int i = 0; i < n2; i++)
        Q[i] = base(b[i] % SQRTMOD, b[i] / SQRTMOD);
    fft(P, 0);
    fft(Q, 0);
    base A1, A2, B1, B2, X, Y;
    for (int i = 0; i < n; i++)
    {
        X = P[i];
        Y = conj(P[(n - i) % n]);
        A1 = (X + Y) * base(0.5, 0);
        A2 = (X - Y) * base(0, -0.5);
        X = Q[i];
        Y = conj(Q[(n - i) % n]);
        B1 = (X + Y) * base(0.5, 0);
        B2 = (X - Y) * base(0, -0.5);
        P1[i] = A1 * B1 + A2 * B2 * base(0, 1);
        Q1[i] = A1 * B2 + A2 * B1;
    }
    for (int i = 0; i < n; i++)
        P[i] = P1[i], Q[i] = Q1[i];
    fft(P, 1);
    fft(Q, 1);
    a.resize(final_size);
    for (int i = 0; i < final_size; i++)
    {
        ll x = (ll)(P[i].real() + 0.5);
        ll y = (ll)(P[i].imag() + 0.5) % MOD;
        ll z = (ll)(Q[i].real() + 0.5);
        a[i] = (x + ((y * SQRTMOD + z) % MOD) * SQRTMOD) % MOD;
    }
}

vll mul(vll a, vll b)
{
    mul_big_mod(a, b);
    return a;
}

int add(int x, int y)
{
    x += y;
    while (x >= MOD) x -= MOD;
    while (x < 0) x += MOD;
    return x;
}

int mul(int x, int y) { return (((ll)x) * y) % MOD; }

int power(int a, ll b)
{
    int x = 1 % MOD;
    while (b)
    {
        if (b & 1)
            x = mul(x, a);
        a = mul(a, a);
        b >>= 1;
    }
    return x;
}

int inv(int a) { return power(a, MOD - 2); }

void pre()
{
    fact[0] = 1;
    for (int i = 1; i < N; i++)
        fact[i] = mul(i, fact[i - 1]);
    ifact[N - 1] = inv(fact[N - 1]);
    for (int i = N - 2; i >= 0; i--)
        ifact[i] = mul(ifact[i + 1], i + 1);
}

int h(ll n, int k)
{
    n %= MOD;
    return mul(power(n + 1, k + 1), ifact[k + 1]);
}

vector<ll> getPoly(int l, int r)
{
    if (l == r)
        return {a[l], 1};
    int mid = (l + r) >> 1;
    return mul(getPoly(l, mid), getPoly(mid + 1, r));
}

ll get(vector<ll> &v, int pos) { return pos >= v.size() ? 0 : v[pos]; }

int main()
{
    pre();
    precompute_powers();
    int n, k;
    ll t = 0;
    cin >> n >> k;
    for (int i = 0; i <= k + 1; i++)
        Inv.push_back(ifact[i]);
    for (int i = 1; i <= n; i++)
    {
        a[i] = i;
        cin >> a[i];
    }
    e = getPoly(1, n);
    reverse(e.begin(), e.end());
    for (int i = 0; i <= n; i++)
        if (!(i & 1))
            e[i] = add(0, -e[i]);
    e.resize(k + 1);
    // calculate b values
    for (int i = 1; i <= k; i += BLOCK)
    {
        int beg = i, en = min(i + BLOCK - 1, k);
        for (int j = beg; j <= en; j++)
        {
            b.resize(j + 1);
            ll val = add(get(temp, j), mul(j, e[j]));

            for (int ind = beg; ind < j; ind++)
            {
                val += b[ind] * e[j - ind];
                if (val >= mod2)
                    val -= mod2;
            }
            b[j] = val % MOD;
        }
        if (en < k)
            temp = mul(b, e);
    }
    b[0] = n;
    temp.clear();
    for (int i = 0; i <= k; i += BLOCK)
    {
        int beg = i, en = min(i + BLOCK - 1, k);
        for (int j = beg; j <= en; j++)
        {
            c.resize(j + 1);
            ll val = add(h(t, j), -get(temp, j + 1));
            for (int ind = beg; ind < j; ind++)
            {
                val -= c[ind] * Inv[j + 1 - ind];
                if (val < 0)
                    val += mod2;
            }
            c[j] = val % MOD;
        }
        if (en < k)
            temp = mul(c, Inv);
    }
    for (int i = 0; i <= k; i++)
    {
        b[i] = mul(b[i], ifact[i]);
    }

    temp = mul(b, c);
    for (int i = 1; i <= k; i++)
        cout << mul(fact[i], temp[i]) << '\n';
    return 0;
}