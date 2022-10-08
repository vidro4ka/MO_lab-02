#include<vector>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>


double null(double isk) {
	if (isk == -0) {
		isk = 0;
		return isk;
	}
	return isk;
}
std::vector<std::vector<double>> input(std::string rezult, std::vector<std::vector<double>>& A, std::vector<double>& c,
	std::vector<double>& b, std::vector<std::string>& sign) {
	std::vector<std::vector<double>> A_test(A.size());
	int max = 0;
	if (c.size() >= A.size() && c.size() >= b.size()) {
		max = c.size();
	}
	if (A.size() >= c.size() && A.size() >= b.size()) {
		max = A.size();
	}
	if (b.size() >= A.size() && b.size() >= c.size()) {
		max = b.size();
	}
	for (size_t i = 0; i < A.size(); ++i) {
		A_test[i].resize(A[i].size() + 1);
		for (size_t j = 0; j < A[i].size(); ++j) {
			A_test[i][j] = A[i][j];
		}
	}
	std::cout << "F(x) = ";
	for (size_t i = 0; i < c.size() - 1; ++i) {
		std::cout << c[i] << "y" << i + 1 << " + ";
	}
	std::cout << c[c.size() - 1] << "y" << c.size();
	std::cout << " --> " << rezult;
	std::cout << std::endl;
	std::cout << "---Initial condition---\n";
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < abs(max) - 1; ++j) {
			std::cout << A[i][j] << "y" << j + 1 << " + ";
		}
		std::cout << A[i][A[i].size() - 1] << "y" << A.size() << " " << sign[i] << " " << b[i];
		std::cout << '\n';
	}
	std::cout << "\n";
	std::cout << "-----canonical form----\n";
	int ch = 1;
	for (size_t i = 0; i < A.size(); ++i) {
		if (sign[i] == "<=" || sign[i] == "<") {
			A[i].push_back(ch);
			sign[i] = "=";
		}
		if (sign[i] == ">=" || sign[i] == ">") {
			A[i].push_back(ch * (-1));
			sign[i] = "=";
		}
	}
	max++;
	for (size_t i = 0; i < A.size(); ++i) {
		if (A[i][0] < 0) {
			std::cout << " - " << -1 * A[i][0] << "y" << 1;
		}
		else {
			std::cout << A[i][0] << "y" << 1;
		}
		for (size_t j = 1; j < abs(max) - 1; ++j) {
			if (A[i][j] < 0) {
				std::cout << " - " << -1 * A[i][j] << "y" << j + 1;
			}
			else {
				std::cout << " + " << A[i][j] << "y" << j + 1;
			}
		}
		if (A[i][A[i].size() - 1] < 0) {
			std::cout << " - " << -1 * A[i][A[i].size() - 1] << "y" << A.size() + i + 1 << " " << sign[i] << " " << b[i];
		}
		else {
			std::cout << " + "<<A[i][A[i].size() - 1] << "y" << A.size() + i + 1 <<  " " << sign[i] << " " << b[i];
		}
		std::cout << '\n';
	}
	std::cout << "-----------------------\n";
	for (size_t i = 0; i < A.size(); ++i) {
		std::swap(A[i][0], A[i][max - 1]);
	}
	bool t = false;
	for (size_t i = 0; i < A.size(); ++i) {
		if (A[i][0] < 0) {
			t = true;
		}
		if (t) {
			std::cout << -1 * A[i][0] << "y" << A.size() + i + 1 << " " << sign[i] << " ";
		}
		else {
			std::cout << A[i][0] << "y" << A.size() + i + 1 << " " << sign[i] << " ";
		}

		A[i][1] = null(A[i][1]);
		if (t) {
			if (A[i][1] < 0) {
				std::cout << " + " << -1 * A[i][1] << "y" << 1 + 1;
			}
			else if (A[i][1] > 0) {
				std::cout << " - " << A[i][1] << "y" << 1 + 1;
			}
			else if (A[i][1] == 0) {
				std::cout << "   " << A[i][1] << "y" << 1 + 1;
			}
		}
		else {
			if (A[i][1] < 0) {
				std::cout << " + " << -1 * A[i][1] << "y" << 1 + 1;
			}
			else if (A[i][1] > 0) {
				std::cout << " - " << A[i][1] << "y" << 1 + 1;
			}
			else if (A[i][1] == 0) {
				std::cout << A[i][1] << "y" << 1 + 1;
			}
		}

		for (size_t j = 2; j < abs(max) - 1; ++j) {
			A[i][j] = null(A[i][j]);
			if (t) {
				if (A[i][j] < 0) {
					std::cout << " + " << -1 * A[i][j] << "y" << j + 1;
				}
				else if(A[i][j] >= 0) {
					std::cout << " - " << A[i][j] << "y" << j + 1;
				}
			}
			else {
				if (A[i][j] < 0) {
					std::cout << " + " << -1 * A[i][j] << "y" << j + 1;
				}
				else if (A[i][j] >= 0) {
					std::cout << " - " << A[i][j] << "y" << j + 1;
				}
			}

		}
		if (t) {
			if (A[i][A[i].size() - 1] >= 0) {
				std::cout << " + " << A[i][A[i].size() - 1] << "y" << A.size() - 2;
			}
			else {
				std::cout << " " << -1 * A[i][A[i].size() - 1] << "y" << A.size() - 2;
			}
			if (b[i] <= 0) {
				std::cout << " + " << -1 * b[i];
			}
			else {
				std::cout << " - " << b[i];
			}
			A_test[i][A[i].size() - 1] = -1 * b[i];
		}
		else {
			if (A[i][A[i].size() - 1] <= 0) {
				std::cout << " + " << A[i][A[i].size() - 1] << "y" << A.size() - 2;
			}
			else {
				std::cout << " " << -1 * A[i][A[i].size() - 1] << "y" << A.size() - 2;
			}
			if (b[i] >= 0) {
				std::cout << " + " << b[i];
			}
			else {
				if (b[i] < 0) {
					std::cout << " - " << -1 * b[i];
				}
			}
			A_test[i][A[i].size() - 1] = b[i];
		}
		t = false;
		std::cout << '\n';
	}
	for (size_t i = 0; i < A.size(); ++i) {
		for (size_t j = 0; j < A.size(); ++j) {
			if (A[i][0] >= 0 && A_test[i][j] != 0) {
				A_test[i][j] = -1 * A_test[i][j];
			}
		}
	}
	int n = A.size() + 1;
	int m = max;
	std::vector<std::vector<double>> sympl(n);
	for (size_t i = 0; i < n - 1; ++i) {
		sympl[i].resize(m);
		for (size_t j = 0; j < m; ++j) {
			sympl[i][j] = A_test[i][j];
		}
		sympl[i][m - 1] = A_test[i][A_test[i].size() - 1];
		std::cout << "\n";
	}
	sympl[n - 1].resize(m);
	for (size_t i = 0; i < m - 1; ++i) {
		sympl[n - 1][i] = 0;
	}
	for (size_t i = 0; i < c.size(); ++i) {
		if (rezult == "max") {
			sympl[n - 1][i] = c[i];
		}
		else if (rezult == "min") {
			sympl[n - 1][i] = -1 * c[i];
		}
	}
	int counter = 0;
	while (counter < m) {
		for (size_t i = 0; i < n; ++i) {
			std::swap(sympl[i][counter], sympl[i][m - 1]);
		}
		counter++;
	}
	for (size_t i = 0; i < n-1; ++i) {
		for (size_t j = 1; j < m; ++j) {
			if (sympl[i][j] != 0) {
				sympl[i][j] = -1 * sympl[i][j];
			}
		}
	}
	return sympl;
}
void printer(std::vector<std::vector<double>>& matr, const std::vector<std::string>& basis,
	const std::vector<std::string>& free) {
	std::cout << "\n\n\n\t ";
	for (size_t i = 0; i < free.size(); ++i) {
		std::cout << free[i] << "	 ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < basis.size(); ++i) {
		std::cout << basis[i] <<"	";
		for (size_t j = 0; j < matr[i].size(); ++j) {
			if (matr[i][j] == -0) {
				matr[i][j] = 0;
			}
			std::cout << std::setw(5) <<std::fixed << std::setprecision(2) << round(matr[i][j] * 100) / 100 << std::right  << "\t";
		}
		std::cout << "\n";
	}
}

int find_column(std::vector<double>& F) {
	int index_r_column = -1;
	for (size_t j = 1; j < F.size(); ++j) {
		if (F[j]>0) {
			index_r_column = j;
			break;
		}
	}
	return index_r_column;
}

int find_row(std::vector<std::vector<double>>& sympl, int index_r_column) {
	int index_r_row = -1;
	double min = INT64_MAX;
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (sympl[i][0] >= 0 && sympl[i][index_r_column] < 0) {
			continue;
		}
		double znach = sympl[i][0] / sympl[i][index_r_column];
		if (znach >= 0 && znach < min) {
			min = znach;
			index_r_row = i;
		}
	}
	return index_r_row;
}

void transformation(std::vector<std::vector<double>>& sympl, int index_r_row, int index_r_column) {
	double razr = sympl[index_r_row][index_r_column];
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (i != index_r_row) {
			for (size_t j = 0; j < sympl[i].size(); ++j) {
				if (j != index_r_column) {
					sympl[i][j] = sympl[i][j] - (sympl[index_r_row][j] * sympl[i][index_r_column] / razr);
				}
			}
		}
	}
	for (size_t i = 0; i < sympl.size(); ++i) {
		if (i != index_r_row) {
			sympl[i][index_r_column] = -1 * sympl[i][index_r_column] / razr;
		}
	}
	for (size_t j = 0; j < sympl[index_r_row].size(); ++j) {
		if (j != index_r_column) {
			sympl[index_r_row][j] = sympl[index_r_row][j] / razr;
		}
	}
	sympl[index_r_row][index_r_column] = 1 / razr;
}

void searher(std::vector<std::vector<double>>& matr, std::vector<std::string>& free, std::vector<std::string>& basis) {
	int index_r_column;
	for (size_t i = 0; i < matr.size() - 1; ++i) {
		if (matr[i][0] < 0) {
			bool otr = false;
			for (size_t j = 1; j < matr[i].size(); ++j) {
				if (matr[i][j] < 0) {
					index_r_column = j;
					size_t index_r_raw = find_row(matr, index_r_column);
					std::swap(free[index_r_column], basis[index_r_raw]);
					transformation(matr, index_r_raw, index_r_column);
					otr = true;
					break;
				}
			}
			if (otr) {
				printer(matr, basis, free);
				i = 0;
			}
			else {
				std::cout << "No solution";
				break;
			}
		}
	}
}

void printer_answer(std::string rezult, std::vector<std::vector<double>>& matr, std::vector<std::string>& free, std::vector<std::string>& basis) {
	std::cout << "\n\n\n---------ANSWER--------\n";
	if (rezult == "max") {
		std::cout << "F  = " << round(-1 * matr[matr.size() - 1][0] * 100) / 100 << std::endl;
	}
	else {
		std::cout << "F  = " << round(matr[matr.size() - 1][0] * 100) / 100 << std::endl;

	}
	for (size_t i = 0; i < matr.size() - 1; ++i) {
		std::cout << basis[i] << " = " << round(matr[i][0] * 100) / 100 << std::endl;
	}
	for (size_t j = 1; j < matr[0].size(); ++j) {
		std::cout << free[j] << " = 0" << std::endl;
	}
}

bool checking(std::vector<std::vector<double>>& matr, std::vector<std::vector<double>>& A, std::vector<double>& c,
	std::vector<double>& b, std::vector<std::string>& sign) {
	bool sign_chek;
	for (size_t i = 0; i < sign.size(); ++i) {
		if (sign[i] == "<=") {
			sign_chek = true;
		}
		if (sign[i] == ">=") {
			sign_chek = false;
		}
	}
	return sign_chek;
}
void method(std::vector<std::vector<double>>& matr, std::string rezult) {
	std::vector<std::string> basis(matr.size());
	std::vector<std::string> free(matr[0].size());
	free[0] = "sv";
	size_t number_of_free = 1;

	for (number_of_free = 1; number_of_free < free.size(); ++number_of_free) {
		free[number_of_free] = "y" + std::to_string(number_of_free);
	}
	for (size_t i = 0; i < basis.size() - 1; ++i) {
		basis[i] = "y" + std::to_string(number_of_free);
		++number_of_free;
	}
	basis[basis.size() - 1] = "F";
	printer(matr, basis, free);
	searher(matr, free, basis);
	int index_r_column = find_column(matr[matr.size() - 1]);
	while (index_r_column != -1) {
		int index_r_row = find_row(matr, index_r_column);
		if (index_r_row == -1) {
			std::cout << "Infinity number of sollution";
			break;
		}
		std::swap(basis[index_r_row], free[index_r_column]);
		transformation(matr, index_r_row, index_r_column);
		printer(matr, basis, free);
		index_r_column = find_column(matr[matr.size() - 1]);
	}
	printer_answer(rezult, matr, free, basis);
}

int main() {
	std::vector<std::vector<double>> x;
	std::vector<double> c = { 6,6,6 };
	std::vector<std::vector<double>> A = { {4,1,1},{1,2,0},{0, 0.5,4} };
	std::vector<double> b = { 5,3,8 };
	std::vector<std::string> sign = { "<=", "<=","<=" };
	std::string rezult = "max";
	/////
	std::vector<double> new_b(c.size());
	std::vector<double> new_c(b.size());
	for (size_t i = 0; i < c.size(); ++i) {
		new_b[i] = c[i];
	}
	for (size_t i = 0; i < b.size(); ++i) {
		new_c[i] = b[i];
	}
	if (rezult == "max") {
		rezult = "min";
	}
	else if (rezult == "min") {
		rezult = "max";
	}
	for (size_t i = 0; i < sign.size(); ++i) {
		if (sign[i] == "<=") {
			sign[i] = ">=";
		}
		else {
			if (sign[i] == "<") {
				sign[i] == ">";
			}
			else {
				if (sign[i] == ">=") {
					sign[i] = "<=";
				}
				else {
					if (sign[i] == ">") {
						sign[i] == "<";
					}
				}
			}
		}
	}
	std::vector<std::vector<double>> trans_A(A[0].size());
	for (size_t i = 0; i < A[0].size(); ++i) {
		trans_A[i].resize(A.size());
		for (size_t j = 0; j < A.size(); ++j) {
			trans_A[i][j] = A[j][i];
		}
	}
	///////////
	x = input(rezult, trans_A, new_c, new_b, sign);
	for (size_t i = 0; i < x.size(); ++i) {
		for (size_t j = 0; j < x[0].size(); ++j) {
			std::cout << x[i][j] << " ";
		}
		std::cout << '\n';
	}
	method(x, rezult);
	return 0;
}