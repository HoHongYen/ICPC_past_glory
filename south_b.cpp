#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
const int MOD = 1e9+7;
const int N = 1e7+5;
ll t, n;
ll fact[N], s1[N], s2[N];

ll add(ll a, ll b) {
    a += b;
    while (a >= MOD) a -= MOD;
    while (a < 0) a += MOD;
    return a;
}

ll mul(ll a, ll b) {
    return a*b % MOD;
}

void pre_compute() {
    fact[0] = 1;
    for (int i=1; i<N; i++) {
        fact[i] = mul(fact[i-1], i);
    }
    
    s1[0] = 1;
    for (int i=1; i<=N; i++) {
        s1[i] = add(s1[i-1], fact[i]);
    }
    
    s2[0]= 1;
    for (int i=1; i<=N; i++) {
        s2[i] = add(s2[i-1], s1[i]);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    cin >> t;
    pre_compute();
    while (t--) {
        cin >> n;
        cout << add(add(fact[n], -2*s1[n-1] + 1), s2[n-2]) << '\n';
    }
    return 0;
}