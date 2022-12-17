#include <iostream>
#include "Matrix.h"
using namespace std;

void input(Matrix &, int, int);
int main()
{
    int rows, columns;
    cout << "Please enter matrix's row and column[row column]: ";
    cin >> rows >> columns;
    Matrix m = Matrix(rows, columns, 0);
    input(m, rows, columns);
    cout << m << endl;
    while (1)
    {
        int r;
        int c;
        int option;
        cout << "1) Transpose 2) Add 3) Multiply 4) Determinant "
             << "5) Adjoint 6) Inverse 7) cofactor 8) exit: ";
        cin >> option;
        if (option == 1)
        {
            Matrix t = m.Transpose();
            cout << t << endl;
        }
        else if (option == 2)
        {
            cout << "[Please enter another "
                 << rows << " * " << columns
                 << " matrix you want to multiply]\n\n";
            cout << "Please enter matrix's row and column[row column]: ";
            cin >> r >> c;
            Matrix a = Matrix(r, c, 0);
            input(a, r, c);
            cout << a << endl;
            Matrix C = m.Add(a);
            cout << C << endl;
        }
        else if (option == 3)
        {
            cout << "[Please enter another "
                 << rows << " *  X "
                 << "matrix you want to multiply]\n\n";
            cout << "Please enter matrix's row and column[row column]: ";
            cin >> r >> c;
            Matrix a = Matrix(r, c, 0);
            input(a, r, c);
            Matrix d = m.Mult(a);
            cout << d << endl;
        }
        else if (option == 4)
        {
            int d = m.Determinant();
            cout << "Determinant: "
                 << d << endl;
        }
        else if (option == 5)
        {
            cout << "Adjoint matrix: " << endl;
            Matrix e = m.Adjoint();
            cout << e << endl;
        }
        else if (option == 6)
        {
            cout << "Inverse matrix: " << endl;
            Matrix f = m.Inverse();
            cout << f << endl;
        }
        else if (option == 7)
        {
            cout << "Cofactor matrix: " << endl;
            Matrix g = m.Cofactor();
            cout << g << endl;
        }
        else if (option == 8)
        {
            cout << "Exit the program!" << endl;
            break;
        }
    }
    return 0;
}

void input(Matrix &m, int rows, int columns)
{
    int row;
    int column;
    do
    {
        int value;
        cout << "Enter non-zero elements of the matrix\n"
             << "[eg. row column value, 999 999 999 to end]: ";
        cin >> row >> column >> value;
        if (row == 999 && column == 999 && value == 999)
        {
            break;
        }
        else if ((column >= columns && column < 0) || (row >= rows && row < 0))
        {
            cout << "out of boundary" << endl;
            continue;
        }
        m.setElement(value, row, column);
    } while (m.getTerms() < rows * columns);
}