#include "grid.h"

#include <fstream>

void read_gg_from_file(std::vector<double>& gg)
{
    std::ifstream in("Au.txt");
    gg = std::vector<double>();
    
    if (in.is_open())
    {
        double data;
        while (in >> data)
        {
            gg.push_back(data);
        }
    }
    in.close();
}

void read_ig_from_file(std::vector<int>& ig)
{
    std::ifstream in("Ia.txt");
    ig = std::vector<int>();
    
    if (in.is_open())
    {
        int data;
        while (in >> data)
        {
            ig.push_back(data);
        }
    }
    in.close();
}

void read_jg_from_file(std::vector<int>& jg)
{
    std::ifstream in("Ja.txt");
    jg = std::vector<int>();
    
    if (in.is_open())
    {
        int data;
        while (in >> data)
        {
            jg.push_back(data);
        }
    }
    in.close();
}

void read_di_from_file(std::vector<double>& di)
{
    std::ifstream in("Di.txt");
    di = std::vector<double>();
    
    if (in.is_open())
    {
        double data;
        while (in >> data)
        {
            di.push_back(data);
        }
    }
    in.close();
}

void read_b_from_file(std::vector<double>& F)
{
    std::ifstream in("B.txt");
    F = std::vector<double>();
    
    if (in.is_open())
    {
        double data;
        while (in >> data)
        {
            F.push_back(data);
        }
    }
    in.close();
}

void read_n_from_file(int& N)
{
    std::ifstream in("N.txt");

    if (in.is_open())
    {
        int data;
        while (in >> data)
        {
            N = data;
        }
    }
    in.close();
}