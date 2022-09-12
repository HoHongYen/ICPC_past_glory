#include <bits/stdc++.h>
using namespace std;

typedef long long ll;

ll power(int x, int y) {
    if (y == 0) return 1;
    ll temp = power(x, y/2);
    if (y&1) return temp*temp*x;
    return temp*temp;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int n;
    cin >> n;
    cout << power(2, n-1);
    return 0;
}