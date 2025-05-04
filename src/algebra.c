#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    Matrix t;
    if(a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    t = create_matrix(a.rows, a.cols);
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < a.cols; j++)
        {
            t.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return t;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    Matrix t;
    if(a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    t = create_matrix(a.rows, a.cols);
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < a.cols; j++)
        {
            t.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return t;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    Matrix t;
    t = create_matrix(a.rows, b.cols);
    if(a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<b.cols;j++)
        {
            t.data[i][j] = 0;
            for(int k=0;k<a.cols;k++)
            {
                t.data[i][j] += a.data[i][k] * b.data[k][j];

            }
        }
    }
    return t;
}

Matrix scale_matrix(Matrix a, double k)
{
    for(int i = 0; i < a.rows; i++)
    {
        for(int j = 0; j < a.cols; j++)
        {
            a.data[i][j] = a.data[i][j]*k;
        }
    }
    return a;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix t;
    t = create_matrix(a.cols, a.rows);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
        {
            t.data[j][i] = a.data[i][j];
        }
    }
    return t;
}

double det_matrix(Matrix a)
{
    if(a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    if(a.rows == 1)
        return a.data[0][0];
    double k=0;
    for(int i=0;i<a.cols;i++)
    {
        Matrix t = create_matrix(a.rows-1, a.cols-1);
        for(int j=1;j<a.rows;j++)
        {
            for(int kk=0;kk<a.cols;kk++)
            {
                if(kk<i)
                    t.data[j-1][kk] = a.data[j][kk];
                else if(kk>i)
                    t.data[j-1][kk-1] = a.data[j][kk];
            }
        }
        k += pow(-1, i) * a.data[0][i] * det_matrix(t);
    }
    return k;
}

Matrix inv_matrix(Matrix a)
{
    Matrix t;
    if(a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    if(det_matrix(a) == 0)
    {
        printf("Error: det(a) == 0.\n");
        return create_matrix(0, 0);
    }
    t = create_matrix(a.rows, a.cols);
    for(int i=0;i<a.rows;i++)
    {
        for(int j=0;j<a.cols;j++)
        {
            Matrix m = create_matrix(a.rows-1, a.cols-1);
            for(int k=0;k<a.rows;k++)
            {
                for(int kk=0;kk<a.cols;kk++)
                {
                    if(k<i && kk<j)
                        m.data[k][kk] = a.data[k][kk];
                    else if(k<i && kk>j)
                        m.data[k][kk-1] = a.data[k][kk];
                    else if(k>i && kk<j)
                        m.data[k-1][kk] = a.data[k][kk];
                    else if(k>i && kk>j)
                        m.data[k-1][kk-1] = a.data[k][kk];
                }
            }
            t.data[i][j] = pow(-1,i+j) * det_matrix(m);
        }
    }
    t = scale_matrix(t, 1/det_matrix(a));
    return t;
}

int rank_matrix(Matrix a)
{
    if(a.rows > a.cols)
        a = transpose_matrix(a);
    int rank = 0;
    for(int i=0;i<a.rows;i++)
    {
        if(a.data[i][i]!=0)
        {
            rank++;
            for(int j=i+1;j<a.rows;j++)
            {
                double k = a.data[j][i] / a.data[i][i];
                for(int kk=0;kk<a.cols;kk++)
                {
                    a.data[j][kk] -= k * a.data[i][kk];
                }
            }
        }
        else
        {
            int k;
            for(k=i+1;k<a.rows;k++)
            {
                if(a.data[k][i]!=0)
                    break;
            }
            if(k == a.rows)
                continue;
            for(int j=i;j<a.cols;j++)
            {
                double t = a.data[i][j];
                a.data[i][j] = a.data[k][j];
                a.data[k][j] = t;
            }
            rank++;
            for(int j=i+1;j<a.rows;j++)
            {
                double k = a.data[j][i] / a.data[i][i];
                for(int kk=0;kk<a.cols;kk++)
                {
                    a.data[j][kk] -= k * a.data[i][kk];
                }
            }
        }
    }
    return rank;
}

double trace_matrix(Matrix a)
{
    double sum = 0;
    if(a.rows != a.cols)
    {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    for(int i=0;i<a.rows;i++)
    {
        sum += a.data[i][i];
    }
    return sum;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}