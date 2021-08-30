#include <iostream>
#include <map>
using namespace std;

int main() {
	map<int, int> coupleID;
	int n, ID_1, ID_2;
	cin >> n;
	while (n--) {
		cin >> ID_1 >> ID_2;
		coupleID.insert(pair<int, int>(ID_1, ID_2));
	}
	cin >> n;
	int** search;
	search = new int*[n];
	for (int i = 0; i < n; ++i) {
		search[i] = new int[2];
	}
	for (int i = 0; i < n; ++i) {
		cin >> search[i][1];
		search[i][0] = 0;
	}
	
	map<int, int>::iterator iter;
	for (int i = 0; i < n; ++i) {
		iter = coupleID.find(search[i][1]);

		if (iter != coupleID.end()) {
			for (int j = 0; j < n; ++j) {
				if (iter->second == search[j][1]) {
					search[i][0] = 1;
					search[j][0] = 1;
					break;
				}
			}
		}
	}

	for (int i = 0; i < n; ++i) {
		if (search[i][0] == 0) {
			cout << search[i][1] << " ";
		}
	}
	
	return 0;
}