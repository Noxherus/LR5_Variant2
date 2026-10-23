#include "Matrix.h"

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error("Self-test failed: " + message);
}
void printMatrix(const std::string& title, const Matrix& value) {
    std::cout << title << ":\n" << value << '\n';
}
template <typename Function>
void demonstrateError(const std::string& title, Function action) {
    try {
        action();
        throw std::runtime_error("Expected exception was not thrown");
    } catch (const std::invalid_argument& error) {
        std::cout << title << ": " << error.what() << '\n';
    } catch (const std::out_of_range& error) {
        std::cout << title << ": " << error.what() << '\n';
    } catch (const std::domain_error& error) {
        std::cout << title << ": " << error.what() << '\n';
    }
}
}

int main() {
    try {
        std::cout << std::boolalpha << std::fixed << std::setprecision(3);
        std::cout << "LABORATORY WORK 5 - VARIANT 2 - ADVANCED LEVEL\n";
        std::cout << "Operator overloading and friend functions (Matrix)\n\n";

        Matrix a{{1, 2}, {3, 4}};
        Matrix b{{2, 0}, {1, 2}};
        Matrix empty;
        Matrix zeroes(2, 2);
        Matrix copy(a);                   // Конструктор копирования.
        Matrix assigned;
        assigned = b;                     // Оператор присваивания.
        Matrix& alias = assigned;
        assigned = alias;                 // Безопасное самоприсваивание.
        check(copy == a && assigned == b, "copy / assignment");
        check(empty.isEmpty() && zeroes.getRows() == 2, "constructors");

        printMatrix("Matrix A", a);
        printMatrix("Matrix B", b);
        printMatrix("A + B", a + b);
        printMatrix("A - B", a - b);
        printMatrix("A * B (matrix multiplication)", a * b);
        printMatrix("A / B (A * inverse(B))", a / b);
        printMatrix("Transpose of A", a.transpose());
        printMatrix("2 * A (friend function)", 2.0 * a);
        check((a * b)[0][0] == 4 && (a * b)[1][1] == 8, "product");
        check((a.transpose())[0][1] == 3, "transpose");

        // Составные операции и инкремент/декремент.
        Matrix arithmetic = a;
        arithmetic += b;
        arithmetic -= b;
        arithmetic *= Matrix::identity(2);
        arithmetic /= Matrix::identity(2);
        arithmetic *= 2;
        arithmetic /= 2;
        check(arithmetic == a, "compound assignment");
        Matrix old = arithmetic++;
        check(old == a && arithmetic[0][0] == 2, "postfix increment");
        ++arithmetic;
        --arithmetic;
        arithmetic--;
        check(arithmetic == a, "increment / decrement");
        std::cout << "Compound operators and ++/--: passed\n\n";

        // Для сравнения выбран лексикографический порядок.
        std::cout << "A == B: " << (a == b) << '\n';
        std::cout << "A != B: " << (a != b) << '\n';
        std::cout << "A <  B: " << (a < b) << '\n';
        std::cout << "A >  B: " << (a > b) << '\n';
        std::cout << "A <= B: " << (a <= b) << '\n';
        std::cout << "A >= B: " << (a >= b) << "\n\n";

        Matrix edited(a);
        edited[0][1] = 7;                 // Доступ через оператор [].
        edited(1, 0) = 8;                 // Доступ через функтор ().
        edited.set(1, 1, 9);             // Сеттер с проверкой.
        check(edited[0][1] == 7 && edited(1, 0) == 8 && edited.get(1, 1) == 9,
              "index / functor / accessors");
        printMatrix("Modified by [] and ()", edited);

        // Явные преобразования.
        const double norm = static_cast<double>(a);
        const bool nonempty = static_cast<bool>(a);
        const auto raw = static_cast<std::vector<std::vector<double>>>(a);
        check(std::abs(norm - std::sqrt(30.0)) < 1e-10, "Frobenius norm");
        check(nonempty && raw[1][1] == 4, "explicit conversions");
        std::cout << "Explicit double conversion (Frobenius norm): " << norm << '\n';
        std::cout << "Explicit bool conversion (nonempty): " << nonempty << "\n\n";

        // Чтение через дружественный оператор >> без ввода с клавиатуры.
        std::istringstream input("2 2 4 3 2 1");
        Matrix parsed;
        input >> parsed;
        check(static_cast<bool>(input) && parsed(1, 1) == 1, "stream input");
        printMatrix("Read using operator >>", parsed);

        // Умножение прямоугольных матриц.
        Matrix rectangularA{{1, 2, 3}, {4, 5, 6}};
        Matrix rectangularB{{7, 8}, {9, 10}, {11, 12}};
        Matrix product = rectangularA * rectangularB;
        check(product.getRows() == 2 && product.getCols() == 2 &&
              product[0][0] == 58 && product[1][1] == 154, "rectangular multiplication");
        printMatrix("Rectangular 2x3 * 3x2", product);

        demonstrateError("Invalid index handled", [&] { (void)a[99][0]; });
        demonstrateError("Invalid dimensions handled", [&] { (void)(a + rectangularA); });
        demonstrateError("Singular matrix handled", [&] {
            Matrix singular{{1, 2}, {2, 4}};
            (void)(a / singular);
        });
        demonstrateError("Division by zero handled", [&] { (void)(a / 0.0); });

        std::cout << "\nProgram completed successfully.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Unexpected error: " << error.what() << '\n';
        return 1;
    }
}
