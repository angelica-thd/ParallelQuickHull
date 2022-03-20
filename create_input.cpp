#include <fstream>
#include <random>
#include <iomanip>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
	int		n;
	double		max_x, max_y;
	ofstream	outputFile;
	random_device	rd;

	if (argc != 5) {
		cout << "Usage: " << argv[0] << " <number of points> <max x value> <max y value> <output file name>" << endl;
		exit(0);
	}

	n     = stoi(argv[1]);
	max_x = stod(argv[2]);
	max_y = stod(argv[3]);

	default_random_engine			eng_x(rd());
	uniform_real_distribution<double>	distr_x(0, max_x);
	default_random_engine			eng_y(rd());
	uniform_real_distribution<double>	distr_y(0, max_y);

	outputFile.open(argv[4]);
	for (auto i = 0; i < n; i++) {
		outputFile << fixed << setprecision(2) << distr_x(eng_x) << " " << distr_y(eng_y) << endl;
	}
	outputFile.close();

	return(0);
}

