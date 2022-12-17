#include <iostream>
#include "Matrix.h"
using namespace std;

int Term::operator>(const Term &s)
{ // override Term's <
    if (row > s.row)
    {
        return 1;
    }
    else if (row == s.row)
    {
        if (col > s.col)
        {
            return 1;
        }
        return 0;
    }
    return 0;
}

int Term::operator<(const Term &s)
{ // override Term's  >
    if (row < s.row)
    {
        return 1;
    }
    else if (row == s.row)
    {
        if (col < s.col)
        {
            return 1;
        }
        return 0;
    }
    return 0;
}

int Term::operator==(const Term &s)
{ // override Term's ==
    if (this == &s)
    {
        return 1;
    }
    if (col == s.col && row == s.row)
    {
        return 1;
    }

    return 0;
}
Matrix::Matrix(int r, int c, int t) : rows(r), cols(c), capacity(t)
{ // A set of triples, <row, column, value>, where row and column are non-negtive
  // integers and form a unique combination; value is also and integer.
    terms = 0;
    if (t != 0)
    {
        smArray = new Term[capacity];
    }
    else
    {
        capacity = 1;
        smArray = new Term[1];
    }
}

Matrix::Matrix(Matrix &s)
{
    rows = s.rows;
    cols = s.cols;
    capacity = s.capacity;
    terms = s.terms;
    smArray = new Term[capacity];
    copy(s.smArray, s.smArray + terms, smArray);
}

Matrix::~Matrix()
{
    delete[] smArray;
}

void Matrix::ChangeSize1D(const int newSize)
{ // Change the size of smArray to newSize
    if (newSize < terms)
        throw e1("New size must be >= number of terms");

    Term *temp = new Term[newSize]; // new array
    copy(smArray, smArray + terms, temp);
    delete[] smArray; // deallocate old memory
    smArray = temp;
    capacity = newSize;
}

void Matrix::setElement(const float sum, const int r, const int c)
{ // If sum or element != 0, then is along with its and column opsition are stored as the last term in *this
    if (sum != 0)
    {
        if (terms == capacity)
        {
            ChangeSize1D(2 * capacity); // double size
        }
        smArray[terms].row = r;
        smArray[terms].col = c;
        smArray[terms++].value = sum;
    }
}

Matrix Matrix::Transpose()
{ // The transpose of a(*this) is placed in b and is found in O(terms + column) time.
    Matrix b(rows, cols, terms);

    if (terms > 0)
    {
        int *rowSize = new int[cols + 1];
        int *rowStart = new int[cols + 1];
        fill(rowSize, rowSize + cols + 1, 0); // initialize
        // Compute rowSize[i] = number of terms in row i of b
        for (int i = 0; i < terms; i++)
        {
            rowSize[smArray[i].col]++;
        }
        // rowStart[i] = starting position of row i in b

        rowStart[0] = 0;
        for (int i = 1; i < cols + 1; i++)
        {
            rowStart[i] = rowStart[i - 1] + rowSize[i - 1];
        }

        for (int i = 0; i < terms; i++)
        { // copy from *this to b

            int j = rowStart[smArray[i].col];
            b.smArray[j].row = smArray[i].col;
            b.smArray[j].col = smArray[i].row;
            b.smArray[j].value = smArray[i].value;
            b.terms++;
            rowStart[smArray[i].col]++;
        } // end of for

        delete[] rowSize;
        delete[] rowStart;
    }
    return b;
}

Matrix Matrix::Add(Matrix b)
{
    if (rows != b.rows || cols != b.cols)
        throw e1("Rows or cols of two matrixs are not equal");

    Matrix c(cols, rows, 0);
    int aPos = 0, bPos = 0;
    while ((aPos < terms) && (bPos < terms))
    {
        if (smArray[aPos] == b.smArray[bPos])
        {
            float t = smArray[aPos].value + b.smArray[bPos].value;
            c.setElement(t, smArray[aPos].row, smArray[aPos].col);
            aPos++;
            bPos++;
        }
        else if (smArray[aPos] < b.smArray[bPos])
        {
            c.setElement(smArray[aPos].value, smArray[aPos].row, smArray[aPos].col);
            aPos++;
        }
        else
        {
            c.setElement(smArray[bPos].value, smArray[bPos].row, smArray[bPos].col);
            bPos++;
        }
    }
    // add in remaining terms of *this
    for (; aPos < terms; aPos++)
    {
        c.setElement(smArray[aPos].value, smArray[aPos].row, smArray[aPos].col);
    }
    // add in remaining terms of b matrix
    for (; bPos < b.terms; bPos++)
    {
        c.setElement(smArray[bPos].value, smArray[bPos].row, smArray[bPos].col);
    }
    return c;
}

Matrix Matrix::Mult(Matrix b)
{ // Return the product of the sparse matrices *his and b
    if (cols != b.rows)
        throw e1("Incompatible matrices");

    Matrix bXpose = b.Transpose();
    Matrix d(rows, cols, 0);
    int currRowIndex = 0,
        currRowBegin = 0,
        currRowA = smArray[0].row;
    // Set boundary conditions
    if (terms == capacity)
    {
        ChangeSize1D(terms + 1);
    }
    bXpose.ChangeSize1D(bXpose.terms + 1);
    smArray[terms].row = rows;
    bXpose.smArray[b.terms].row = b.cols;
    bXpose.smArray[b.terms].col = -1;
    int sum = 0;
    while (currRowIndex < terms)
    { // generate row currentRowA of d
        int currColB = bXpose.smArray[0].row;
        int currColIndex = 0;
        while (currColIndex <= b.terms)
        { // multiply row currRowA of *this by column crrColB of b
            if (smArray[currRowIndex].row != currRowA)
            { // end of row currRowA
                d.setElement(sum, currRowA, currColB);
                sum = 0; // reset sum
                currRowIndex = currRowBegin;
                // advance to next column
                while (bXpose.smArray[currColIndex].row == currColB)
                {
                    currColIndex++;
                }
                currColB = bXpose.smArray[currColIndex].row;
            }
            else if (bXpose.smArray[currColIndex].row != currColB)
            { // end of column currColB of b
                d.setElement(sum, currRowA, currColB);
                sum = 0; // reset sum
                // set multiply row currRowA with next column
                currColIndex = currRowBegin;
                currColB = bXpose.smArray[currColIndex].row;
            }
            else
            {
                if (smArray[currRowIndex].col < bXpose.smArray[currColIndex].col)
                {
                    currRowIndex++; // advance to next term in row
                }
                else if (smArray[currRowIndex].col == bXpose.smArray[currColIndex].col)
                { // add sum
                    sum += smArray[currRowIndex].value * bXpose.smArray[currColIndex].value;
                    currRowIndex++;
                    currColIndex++;
                }
                else
                {
                    currColIndex++; // next term in currColB
                }
            }
        }                                             // end of while(currColIndex <= b.terms)
        while (smArray[currRowIndex].row == currRowA) // advance to next row
        {
            currRowIndex++;
        }
        currRowBegin = currRowIndex;
        currRowA = smArray[currRowIndex].row;
    } // end of while(currRowIndex , terms)
    return d;
}

Matrix Matrix::minor(int q, int p)
{
    Matrix d(rows - 1, cols - 1, 0);
    int index = 0;
    int k = 0;
    for (int i = 0; i < cols * rows; i++)
    {
        if (index >= (rows - 1) * (cols - 1) || k >= terms)
        {
            break;
        }
        else if (smArray[k].row == i / rows && smArray[k].col == i % cols)
        {
            if (smArray[k].row != q && smArray[k].col != p)
            {
                d.setElement(smArray[k].value, index / (rows - 1), index % (cols - 1));
                k++;
                index++;
            }
            else
            {
                k++;
            }
        }
        else
        {
            if ((i / rows) != q && (i % cols) != p)
                index++;
        }
    }
    return d;
}

int Matrix::Determinant()
{
    if (cols != rows)
        throw e1("not a square matrix");
    int D = 0;
    if (cols == 1)
    {
        if (terms != 0)
        {
            return smArray[0].value;
        }
        return 0;
    }
    int sign = 1; // To store sign multiplier

    int currIndex = 0;
    while (currIndex < terms)
    {
        if (smArray[currIndex].row == 0)
        {
            Matrix d = minor(0, smArray[currIndex].col);
            if (smArray[currIndex].col % 2 == 1)
            {
                sign = -sign;
            }
            D += sign * smArray[currIndex].value * d.Determinant();
            sign = 1;
        }
        currIndex++;
    }
    return D;
}

Matrix Matrix::Cofactor()
{
    if (cols != rows)
        throw e1("not a square matrix");
    Matrix d(rows, cols, 0);
    int sign = 1;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if ((i + j) % 2 == 1)
            {
                sign = -sign;
            }
            Matrix temp = minor(i, j);
            int co = sign * temp.Determinant();
            if (co != 0)
            {
                d.setElement(co, i, j);
            }
            sign = 1;
        }
    }
    return d;
}

Matrix Matrix::Adjoint()
{
    Matrix temp = Cofactor().Transpose();
    Matrix d(rows, cols, 0);
    for (int i = 0; i < temp.terms; i++)
        d.setElement(temp.smArray[i].value, temp.smArray[i].row, temp.smArray[i].col);
    return d;
}

Matrix Matrix::Inverse()
{
    if (cols != rows)
        throw e1("not a square matrix");

    Matrix temp = Adjoint();
    // Matrix d(temp.rows, temp.cols, 0);
    int D = Determinant();
    for (int i = 0; i < terms; i++)
    {
        // d.setElement((temp.smArray[i].value / D), temp.smArray[i].row, temp.smArray[i].col);
        temp.smArray[i].value = float(temp.smArray[i].value / D);
        // cout << d.smArray[i].value << endl;
    }
    return temp;
}

ostream &operator<<(ostream &os, Matrix &b)
{
    os << "\n"
       << "matrix dimensions:" << b.rows << " X " << b.cols << "\n"
       << "\n"
       << "      \t";

    for (int i = 0; i < b.cols; i++)
    {
        os << "col[" << i << "]\t";
    }
    os << "\n\n";

    int k = 0;
    for (int i = 0; i < b.rows; i++)
    {
        os << "row[" << i << "]\t";
        for (int j = 0; j < b.cols; j++)
        {
            // whether the element in the matrix is a none-zero term
            // cout << i << " " << j << "  " << k << endl;
            if (b.smArray[k].row == i && b.smArray[k].col == j && k < b.terms)
            {
                os << b.smArray[k].value << "\t";
                k++;
            }
            else
            {
                os << "0\t";
            }
        }
        os << "\n\n";
    }

    return os;
}

std::istream &operator>>(std::istream &is, Matrix &b)
{
    int r, c, v;
    is >> v >> c >> r;
    b.setElement(r, c, v);
    return is;
}