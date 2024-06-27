﻿#include "solver.h"

#include <fstream>
#include <iostream>

/// <summary>LU(sq) разложение</summary>
void solver::LU()
{
    //копирование-инициализация
    for (int i = 0; i < di_n.size(); i++)
    {
        double sd = 0; //переменные суммирования
        int i0 = ig[i];
        int i1 = ig[i + 1];
        for (int k = i0; k < i1; k++)
        {
            int j = jg[k];
            double sl = 0;
            double su = 0;
            int j0 = ig[j];
            int j1 = ig[j + 1];
            int ki = i0;
            int kj = j0;
            while (ki < k && kj < j1)
            {
                int jl = jg[ki];
                int ju = jg[kj];
                if (jl == ju)
                {
                    sl += ggu_n[kj] * ggl_n[ki];
                    su += ggl_n[kj] * ggu_n[ki];
                    ki++;
                    kj++;
                }
                else if (jl < ju) ki++;
                else kj++;
            }
            ggl_n[k] -= sl;
            ggu_n[k] = (ggu_n[k] - su) / di_n[j];
            sd += ggu_n[k] * ggl_n[k];
        }
        di_n[i] -= sd;
    }
}

/// <summary>y = L^-1 * b</summary>
void solver::step_direct(
    std::vector<double>& ggl,
    std::vector<double>& di,
    std::vector<double>& res,
    std::vector<double>& b
)
{
    for (int i = 0; i < res.size(); i++)
    {
        double sum = 0;
        for (int j = ig[i]; j < ig[i + 1]; j++)
            sum += ggl[j] * res[jg[j]];
        res[i] = (b[i] - sum) / di[i];
    }
}

/// <summary>x = U^-1 * y</summary>
void solver::step_reverse(
    std::vector<double>& ggu,
    std::vector<double>& di,
    std::vector<double>& res,
    std::vector<double>& b
)
{
    res = b;
    for (int i = res.size() - 1; i >= 0; i--)
    {
        for (int j = ig[i + 1] - 1; j >= ig[i]; j--)
            res[jg[j]] -= res[i] * ggu[j];
    }
}

std::vector<double> solver::mult(const std::vector<double>& v)
{
    std::vector<double> res(v.size(), 0);
    for (int i = 0; i < v.size(); i++)
    {
        res[i] = di[i] * v[i];
        for (int j = ig[i]; j < ig[i + 1]; j++)
        {
            res[i] += gg[j] * v[jg[j]];
            res[jg[j]] += gg[j] * v[i];
        }
    }
    return res;
}


/// <summary>Стабилизированный метод бисопряженных градиентов</summary>
std::vector<double> solver::BCGSTAB(
    std::vector<double>& _gg,
    std::vector<double>& _di,
    std::vector<double>& F,
    std::vector<int>& _ig,
    std::vector<int>& _jg
)
{
    di = _di;
    gg = _gg;
    di_n = di;
    ggl_n = gg;
    ggu_n = gg;
    ig = _ig;
    jg = _jg;
    LU();
    int size = F.size();
    q.resize(size);
    std::vector<double> temp1(size);
    std::vector<double> temp2(size);
    std::vector<double> z_sl(size);
    std::vector<double> y_sl(size);
    std::vector<double> r0(size);
    r.resize(size);
    p.resize(size);

    // инициализация
    temp1 = mult(q);
    temp2 = F - temp1;
    step_direct(ggl_n, di_n, r0, temp2);
    z_sl = r0;
    r = r0;

    //iteration
    double nev = sqrt(r0 * r0);
    int iter = 0;
    for (; iter < maxIter && nev >= eps; iter++)
    {
        step_reverse(ggu_n, di_n, temp1, z_sl);
        temp1 = mult(temp1);
        step_direct(ggl_n, di_n, temp1, temp1);

        double r_r0 = r * r0;
        double alpha = r_r0 / (temp1 * r0);

        p = r - alpha * temp1;

        step_reverse(ggu_n, di_n, temp2, p);
        temp2 = mult(temp2);
        step_direct(ggl_n, di_n, temp2, temp2);

        double gamma = (p * temp2) / (temp2 * temp2);

        y_sl = y_sl + alpha * z_sl + gamma * p;

        r = p - gamma * temp2;

        double beta = alpha / (gamma * r_r0);
        beta *= (r * r0);

        z_sl = r + beta * z_sl - beta * gamma * temp1;

        step_reverse(ggu_n, di, q, y_sl);

        nev = sqrt(r * r);
    }
    std::cout << "iterations: " << iter << '\n';

    r.clear();
    p.clear();
    z_sl.clear();
    y_sl.clear();
    temp1.clear();
    temp2.clear();
    return q;
}

double solver::check_inaccuracy(std::vector<double>& q, const char* str, int n)
{
    std::ifstream in(str);
    std::vector<double> pardiso_solve = std::vector<double>();

    if (in.is_open())
    {
        double data;
        while (in >> data)
        {
            pardiso_solve.push_back(data);
        }
    }
    in.close();

    double s1 = 0, s2 = 0;
    for (size_t i = 0; i < n; i++)
    {
        s1 += (pardiso_solve[i] - q[i]) * (pardiso_solve[i] - q[i]);
        s2 += pardiso_solve[i] * pardiso_solve[i];
    }
    return s1 / s2;
}