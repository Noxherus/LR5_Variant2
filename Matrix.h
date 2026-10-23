#pragma once

#include <cstddef>
#include <iosfwd>
#include <initializer_list>
#include <vector>

// Класс матрицы с безопасным доступом и перегруженными операторами.
class Matrix {
public:
    class RowProxy {
    public:
        double& operator[](std::size_t column) const;
    private:
        friend class Matrix;
        RowProxy(Matrix& owner, std::size_t row) : owner_(owner), row_(row) {}
        Matrix& owner_;
        std::size_t row_;
    };

    class ConstRowProxy {
    public:
        const double& operator[](std::size_t column) const;
    private:
        friend class Matrix;
        ConstRowProxy(const Matrix& owner, std::size_t row) : owner_(owner), row_(row) {}
        const Matrix& owner_;
        std::size_t row_;
    };

    Matrix();
    Matrix(std::size_t rows, std::size_t columns, double value = 0.0);
    Matrix(std::initializer_list<std::initializer_list<double>> values);
    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    ~Matrix();

    std::size_t getRows() const noexcept;
    std::size_t getCols() const noexcept;
    bool isEmpty() const noexcept;
    double get(std::size_t row, std::size_t col) const;
    void set(std::size_t row, std::size_t col, double value);
    void swap(Matrix& other) noexcept;

    double& at(std::size_t row, std::size_t col);
    const double& at(std::size_t row, std::size_t col) const;
    RowProxy operator[](std::size_t row);
    ConstRowProxy operator[](std::size_t row) const;
    double& operator()(std::size_t row, std::size_t col);
    const double& operator()(std::size_t row, std::size_t col) const;

    Matrix transpose() const;
    Matrix inverse() const;
    static Matrix identity(std::size_t n);

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator-() const;
    Matrix operator*(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix operator/(const Matrix& other) const;
    Matrix operator/(double scalar) const;

    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix& operator*=(const Matrix& other);
    Matrix& operator*=(double scalar);
    Matrix& operator/=(const Matrix& other);
    Matrix& operator/=(double scalar);

    bool operator==(const Matrix& other) const noexcept;
    bool operator!=(const Matrix& other) const noexcept;
    bool operator<(const Matrix& other) const noexcept;
    bool operator>(const Matrix& other) const noexcept;
    bool operator<=(const Matrix& other) const noexcept;
    bool operator>=(const Matrix& other) const noexcept;

    Matrix& operator++();
    Matrix operator++(int);
    Matrix& operator--();
    Matrix operator--(int);

    // Явные преобразования: непустая матрица, норма Фробениуса, копия данных.
    explicit operator bool() const noexcept;
    explicit operator double() const;
    explicit operator std::vector<std::vector<double>>() const;

    friend Matrix operator*(double scalar, const Matrix& matrix);
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
    friend std::istream& operator>>(std::istream& is, Matrix& matrix);

private:
    std::vector<std::vector<double>> data_;
    std::size_t rows_;
    std::size_t cols_;

    static void validateDimensions(std::size_t rows, std::size_t cols);
    void requireSameShape(const Matrix& other) const;
    void requireIndex(std::size_t row, std::size_t col) const;
    void ensureFinite() const;
};
