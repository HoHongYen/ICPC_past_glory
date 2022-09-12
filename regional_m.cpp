#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
const int MOD = 1e9+7;

ll mul(ll a, ll b) {
    return (a*b) % MOD;
}

ll power(ll a, ll b) {
    ll res = 1;
    while (b) {
        if (b&1) res = mul(res, a);
        a = mul(a, a);
        b >>= 1;
    }
    return res;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    ll n, k;
    cin >> n >> k;
    cout << mul(power(2, k-1), (n-k+1));
    return 0;
}