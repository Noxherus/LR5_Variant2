#include "Matrix.h"

#include <algorithm>
#include <cmath>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <utility>

namespace {
constexpr std::size_t kMaxDimension = 10000;
constexpr std::size_t kMaxElements = 1000000;

void requireFinite(double value) {
    if (!std::isfinite(value)) {
        throw std::overflow_error("A matrix contains a non-finite value");
    }
}
}  // namespace

void Matrix::validateDimensions(std::size_t rows, std::size_t cols) {
    if ((rows == 0) != (cols == 0)) {
        throw std::invalid_argument("Both dimensions must be zero or positive");
    }
    if (rows > kMaxDimension || cols > kMaxDimension ||
        (rows != 0 && cols > kMaxElements / rows)) {
        throw std::length_error("Matrix dimensions are too large");
    }
}

Matrix::Matrix() : rows_(0), cols_(0) {}

Matrix::Matrix(std::size_t rows, std::size_t cols, double value)
    : rows_(0), cols_(0) {
    validateDimensions(rows, cols);
    requireFinite(value);
    data_.assign(rows, std::vector<double>(cols, value));
    rows_ = rows;
    cols_ = cols;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> values)
    : rows_(0), cols_(0) {
    if (values.size() == 0) return;
    const std::size_t cols = values.begin()->size();
    if (cols == 0) throw std::invalid_argument("Rows cannot be empty");
    validateDimensions(values.size(), cols);
    for (const auto& row : values) {
        if (row.size() != cols) throw std::invalid_argument("Rows have different lengths");
        data_.emplace_back(row);
    }
    rows_ = values.size();
    cols_ = cols;
    ensureFinite();
}

Matrix::Matrix(const Matrix& other)
    : data_(other.data_), rows_(other.rows_), cols_(other.cols_) {}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        Matrix temporary(other);
        swap(temporary);
    }
    return *this;
}

Matrix::~Matrix() = default;

std::size_t Matrix::getRows() const noexcept { return rows_; }
std::size_t Matrix::getCols() const noexcept { return cols_; }
bool Matrix::isEmpty() const noexcept { return rows_ == 0; }
double Matrix::get(std::size_t row, std::size_t col) const { return at(row, col); }
void Matrix::set(std::size_t row, std::size_t col, double value) {
    requireFinite(value);
    at(row, col) = value;
}

void Matrix::swap(Matrix& other) noexcept {
    data_.swap(other.data_);
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
}

void Matrix::requireIndex(std::size_t row, std::size_t col) const {
    if (row >= rows_ || col >= cols_) throw std::out_of_range("Matrix index is out of range");
}

double& Matrix::at(std::size_t row, std::size_t col) {
    requireIndex(row, col);
    return data_[row][col];
}
const double& Matrix::at(std::size_t row, std::size_t col) const {
    requireIndex(row, col);
    return data_[row][col];
}
Matrix::RowProxy Matrix::operator[](std::size_t row) {
    if (row >= rows_) throw std::out_of_range("Matrix row is out of range");
    return RowProxy(*this, row);
}
Matrix::ConstRowProxy Matrix::operator[](std::size_t row) const {
    if (row >= rows_) throw std::out_of_range("Matrix row is out of range");
    return ConstRowProxy(*this, row);
}
double& Matrix::RowProxy::operator[](std::size_t column) const { return owner_.at(row_, column); }
const double& Matrix::ConstRowProxy::operator[](std::size_t column) const {
    return owner_.at(row_, column);
}
double& Matrix::operator()(std::size_t row, std::size_t col) { return at(row, col); }
const double& Matrix::operator()(std::size_t row, std::size_t col) const { return at(row, col); }

void Matrix::ensureFinite() const {
    for (const auto& row : data_) {
        for (double v : row) requireFinite(v);
    }
}

void Matrix::requireSameShape(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_)
        throw std::invalid_argument("Matrix dimensions must match");
}

Matrix Matrix::transpose() const {
    Matrix result(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t j = 0; j < cols_; ++j)
            result.data_[j][i] = data_[i][j];
    return result;
}

Matrix Matrix::identity(std::size_t n) {
    if (n == 0) throw std::invalid_argument("Identity matrix must be nonempty");
    Matrix result(n, n);
    for (std::size_t i = 0; i < n; ++i) result.data_[i][i] = 1;
    return result;
}

Matrix Matrix::inverse() const {
    if (isEmpty() || rows_ != cols_)
        throw std::invalid_argument("Only nonempty square matrices have an inverse");
    ensureFinite();
    Matrix work(*this);
    Matrix result = identity(rows_);
    double largest = 0.0;
    for (const auto& row : data_)
        for (double x : row) largest = std::max(largest, std::abs(x));
    if (largest == 0.0) throw std::domain_error("Matrix is singular");
    const double tolerance = 1e-12 * largest;

    // Метод Гаусса-Жордана с частичным выбором главного элемента.
    for (std::size_t col = 0; col < rows_; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < rows_; ++row)
            if (std::abs(work.data_[row][col]) > std::abs(work.data_[pivot][col]))
                pivot = row;
        if (std::abs(work.data_[pivot][col]) <= tolerance)
            throw std::domain_error("Matrix is singular or numerically unstable");
        if (pivot != col) {
            work.data_[pivot].swap(work.data_[col]);
            result.data_[pivot].swap(result.data_[col]);
        }
        const double divisor = work.data_[col][col];
        for (std::size_t j = 0; j < rows_; ++j) {
            work.data_[col][j] /= divisor;
            result.data_[col][j] /= divisor;
        }
        for (std::size_t row = 0; row < rows_; ++row) {
            if (row == col) continue;
            const double factor = work.data_[row][col];
            for (std::size_t j = 0; j < rows_; ++j) {
                work.data_[row][j] -= factor * work.data_[col][j];
                result.data_[row][j] -= factor * result.data_[col][j];
            }
        }
    }
    result.ensureFinite();
    return result;
}

Matrix Matrix::operator+(const Matrix& other) const {
    requireSameShape(other);
    Matrix result(*this);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t j = 0; j < cols_; ++j)
            result.data_[i][j] += other.data_[i][j];
    result.ensureFinite();
    return result;
}
Matrix Matrix::operator-(const Matrix& other) const {
    requireSameShape(other);
    Matrix result(*this);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t j = 0; j < cols_; ++j)
            result.data_[i][j] -= other.data_[i][j];
    result.ensureFinite();
    return result;
}
Matrix Matrix::operator-() const { return *this * (-1.0); }

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_ || isEmpty() || other.isEmpty())
        throw std::invalid_argument("Incompatible matrix multiplication dimensions");
    Matrix result(rows_, other.cols_);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t k = 0; k < cols_; ++k)
            for (std::size_t j = 0; j < other.cols_; ++j)
                result.data_[i][j] += data_[i][k] * other.data_[k][j];
    result.ensureFinite();
    return result;
}
Matrix Matrix::operator*(double scalar) const {
    requireFinite(scalar);
    Matrix result(*this);
    for (auto& row : result.data_)
        for (double& element : row) element *= scalar;
    result.ensureFinite();
    return result;
}
Matrix operator*(double scalar, const Matrix& matrix) { return matrix * scalar; }

// Деление на матрицу определено как A / B = A * inverse(B).
Matrix Matrix::operator/(const Matrix& other) const { return *this * other.inverse(); }
Matrix Matrix::operator/(double scalar) const {
    requireFinite(scalar);
    if (scalar == 0.0) throw std::domain_error("Division by zero");
    Matrix result(*this);
    for (auto& row : result.data_)
        for (double& element : row) element /= scalar;
    result.ensureFinite();
    return result;
}

Matrix& Matrix::operator+=(const Matrix& other) { *this = *this + other; return *this; }
Matrix& Matrix::operator-=(const Matrix& other) { *this = *this - other; return *this; }
Matrix& Matrix::operator*=(const Matrix& other) { *this = *this * other; return *this; }
Matrix& Matrix::operator*=(double scalar) { *this = *this * scalar; return *this; }
Matrix& Matrix::operator/=(const Matrix& other) { *this = *this / other; return *this; }
Matrix& Matrix::operator/=(double scalar) { *this = *this / scalar; return *this; }

bool Matrix::operator==(const Matrix& other) const noexcept {
    return rows_ == other.rows_ && cols_ == other.cols_ && data_ == other.data_;
}
bool Matrix::operator!=(const Matrix& other) const noexcept { return !(*this == other); }

// Для матриц используется лексикографический порядок, не матричное неравенство.
bool Matrix::operator<(const Matrix& other) const noexcept {
    if (rows_ != other.rows_) return rows_ < other.rows_;
    if (cols_ != other.cols_) return cols_ < other.cols_;
    return data_ < other.data_;
}
bool Matrix::operator>(const Matrix& other) const noexcept { return other < *this; }
bool Matrix::operator<=(const Matrix& other) const noexcept { return !(other < *this); }
bool Matrix::operator>=(const Matrix& other) const noexcept { return !(*this < other); }

Matrix& Matrix::operator++() {
    for (auto& row : data_) for (double& x : row) ++x;
    ensureFinite();
    return *this;
}
Matrix Matrix::operator++(int) { Matrix before(*this); ++(*this); return before; }
Matrix& Matrix::operator--() {
    for (auto& row : data_) for (double& x : row) --x;
    ensureFinite();
    return *this;
}
Matrix Matrix::operator--(int) { Matrix before(*this); --(*this); return before; }

Matrix::operator bool() const noexcept { return !isEmpty(); }
Matrix::operator double() const {
    double norm = 0;
    for (const auto& row : data_)
        for (double x : row) norm = std::hypot(norm, x);
    requireFinite(norm);
    return norm;
}
Matrix::operator std::vector<std::vector<double>>() const { return data_; }

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    os << matrix.rows_ << 'x' << matrix.cols_ << '\n';
    for (const auto& row : matrix.data_) {
        for (std::size_t j = 0; j < row.size(); ++j) {
            if (j) os << ' ';
            os << row[j];
        }
        os << '\n';
    }
    return os;
}

// Ввод: число строк, столбцов, затем элементы построчно.
std::istream& operator>>(std::istream& is, Matrix& matrix) {
    long long rows = 0, cols = 0;
    if (!(is >> rows >> cols)) return is;
    if (rows <= 0 || cols <= 0 || rows > static_cast<long long>(kMaxDimension) ||
        cols > static_cast<long long>(kMaxDimension) ||
        rows > static_cast<long long>(kMaxElements) / cols) {
        is.setstate(std::ios::failbit);
        return is;
    }
    Matrix temp(static_cast<std::size_t>(rows), static_cast<std::size_t>(cols));
    for (std::size_t i = 0; i < temp.rows_; ++i) {
        for (std::size_t j = 0; j < temp.cols_; ++j) {
            if (!(is >> temp.data_[i][j]) || !std::isfinite(temp.data_[i][j])) {
                is.setstate(std::ios::failbit);
                return is;
            }
        }
    }
    matrix.swap(temp);
    return is;
}
