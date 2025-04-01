//Primero: un programa que lea matrices

//Luego las matrices de reescala

//Probamos optimización. O bien ompmD o O1-O2-O3.

//Si lo hacemos en python lo hacemos con numba.

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

class Matrix {
private:
    std::vector<std::vector<int>> data;
    
public:
    Matrix(size_t rows, size_t cols) : data(rows, std::vector<int>(cols, 0)) {}
    
    size_t rows() const { return data.size(); }
    size_t cols() const { return data.empty() ? 0 : data[0].size(); }
    
    int& operator()(size_t row, size_t col) {
        if (row >= rows() || col >= cols()) {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        return data[row][col];
    }
    
    const int& operator()(size_t row, size_t col) const {
        if (row >= rows() || col >= cols()) {
            throw std::out_of_range("Matrix indices out of bounds");
        }
        return data[row][col];
    }
};

Matrix read_matrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    size_t rows, cols;
    file >> rows >> cols;
    if (file.fail()) {
        throw std::runtime_error("Invalid matrix dimensions in file: " + filename);
    }

    Matrix matrix(rows, cols);
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            file >> matrix(i, j);
            if (file.fail()) {
                throw std::runtime_error(
                    "Invalid matrix element at row " + std::to_string(i + 1) +
                    ", column " + std::to_string(j + 1));
            }
        }
    }

    return matrix;
}

void print_matrix(const Matrix& matrix) {
    size_t rows = matrix.rows();
    size_t cols = matrix.cols();
    
    std::cout << "Matrix dimensions: " << rows << "x" << cols << "\n";
    std::cout << "Matrix contents:\n";
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            std::cout << std::setw(4) << matrix(i, j);
        }
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <directory_path>\n";
        return 1;
    }

    try {
        fs::path directory(argv[1]);
        
        if (!fs::exists(directory)) {
            throw std::runtime_error("Directory does not exist: " + directory.string());
        }
        
        if (!fs::is_directory(directory)) {
            throw std::runtime_error("Path is not a directory: " + directory.string());
        }

        std::cout << "Reading matrices from directory: " << directory.string() << "\n";

        for (const auto& entry : fs::directory_iterator(directory)) {
            if (entry.is_regular_file()) {
                std::cout << "\nReading file: " << entry.path().filename().string() << "\n";
                
                try {
                    Matrix matrix = read_matrix(entry.path().string());
                    print_matrix(matrix);
                } catch (const std::exception& e) {
                    std::cerr << "Error processing file: " << entry.path().filename().string() 
                             << "\nError: " << e.what() << "\n";
                }
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}