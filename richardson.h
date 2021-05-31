#ifndef richardson
#define richardson

#include <bits/stdc++.h>

template <class T>
class rich_exp
{
public:
    // richardson method constructors
    rich_exp(long double (*func)(long double obj));

    //differentiation at point x0
    long double diff(T x0);

    //integration from l to r
    long double intg(T l, T r);

private:
    long double h_init = 0.000000001;            //intial step size
    int max_itr = 10;                            //limiting iterations
    long double epsilon = 0.000000000001;        //convergence criterion
    long double (*func)(long double obj);        //mathematical function
    long long power(long long x1, long long y1); //power function using binary exponentiation
};

template <class T>
rich_exp<T>::rich_exp(long double (*func)(long double obj))
{
    this->func = func;
}

template <class T>
long double rich_exp<T>::diff(T x0)
{
    long double mat[max_itr][max_itr]; //recursion matrix
    long double h_arr[max_itr];        //step size array
    memset(mat, 0.0, sizeof(mat));
    memset(h_arr, 0.0, sizeof(h_arr));
    h_arr[0] = h_init;
    for (int i = 0; i < max_itr - 1; i++)
    {
        h_arr[i + 1] = h_arr[i] / 2.0;
    }

    //central difference formula
    for (int i = 0; i < max_itr; i++)
    {
        mat[i][0] = (func(x0 + h_arr[i]) - func(x0 - h_arr[i])) / (2.0 * h_arr[i]);
    }

    int pos = 0;
    for (int i = 1; i < max_itr; i++)
    {
        for (int j = i; j < max_itr; j++)
        {
            //recurrence relation
            mat[j][i] = ((long double)power(4, i) * mat[j][i - 1] - mat[j - 1][i - 1]) / (long double)(power(4, i) - 1.0);
        }
        long double rel_error = abs(mat[i][i] - mat[i - 1][i - 1]) / abs(mat[i - 1][i - 1]);
        if (rel_error < epsilon)
        {
            //numerically converged
            pos = i;
            break;
        }
    }
    return mat[pos][pos];
}

template <class T>
long double rich_exp<T>::intg(T l, T r)
{
    long double mat[max_itr][max_itr]; //recursion matrix
    long double h_arr[max_itr];        //step size array
    memset(mat, 0.0, sizeof(mat));
    memset(h_arr, 0.0, sizeof(h_arr));
    long double h_init_temp = 2;
    h_arr[0] = h_init_temp;
    for (int i = 0; i < max_itr - 1; i++)
    {
        h_arr[i + 1] = h_arr[i] / 2.0;
    }

    //trapezoidal rule
    for (int i = 0; i < max_itr; i++)
    {
        long double ans = func(l) + func(r);
        int np = (r - l) / h_arr[i];
        for (int j = 1; j <= np - 1; j++)
        {
            ans = ans + 2.0 * (func(l + j * h_arr[i]));
        }
        ans = ans * h_arr[i];
        ans = ans / 2.0;
        mat[i][0] = ans;
    }

    int pos = -1;
    for (int i = 1; i < max_itr; i++)
    {
        for (int j = i; j < max_itr; j++)
        {
            //recurrence relation
            mat[j][i] = ((long double)power(4, i) * mat[j][i - 1] - mat[j - 1][i - 1]) / (long double)(power(4, i) - 1.0);
        }
        long double rel_error = abs(mat[i][i] - mat[i - 1][i - 1]) / abs(mat[i - 1][i - 1]);
    }

    return mat[max_itr - 1][max_itr - 1];
}

// power function using binary exponentiation
template <class T>
long long rich_exp<T>::power(long long x, long long y)
{
    long long res = 1;
    x = x;
    if (x == 0)
        return 0;
    while (y > 0)
    {
        if (y & 1)
            res = (res * x);
        y = y >> 1;
        x = (x * x);
    }
    return res;
}

#endif //richardson
