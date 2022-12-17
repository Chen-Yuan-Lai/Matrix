#ifndef MATRIX_H
#define MATRIX_H
using namespace std;

class Matrix;

class Term
{
    friend Matrix;

public:
    int operator>(const Term &s);
    int operator<(const Term &s);
    int operator==(const Term &s);
    friend ostream &operator<<(ostream &os, Matrix &b);

private:
    int row, col;
    double value;
};

class Matrix
{
public:
    Matrix(int r, int c, int t);
    Matrix(Matrix &s);
    ~Matrix();
    void setElement(const float sum, const int r, const int c);
    Matrix Transpose();
    Matrix Add(Matrix b);
    Matrix Mult(Matrix b);
    Matrix minor(int q, int p);
    int Determinant();
    Matrix Adjoint();
    Matrix Inverse();
    Matrix Cofactor();
    void ChangeSize1D(const int newSize);
    int getTerms() const { return terms; }
    friend ostream &operator<<(ostream &os, Matrix &b);
    friend istream &operator>>(istream &is, Matrix &b);

private:
    int rows, cols, terms, capacity;
    Term *smArray;
};

// define custom exception type
class e1 : public exception
{
public:
    e1(string msg) : message(msg) {}
    // override std:: what()
    string content()
    {
        return message;
    }

private:
    string message;
};

#endif