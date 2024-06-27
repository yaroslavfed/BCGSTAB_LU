#include <fstream>
#include <iostream>
#include <vector>

#include "solver.h"
#include "grid.h"
#include "memory_counter.h"

int main(int argc, char* argv[])
{
    std::vector<double> di, gg, F;
    std::vector<int> ig, jg;
    int N;

    read_gg_from_file(gg);
    read_ig_from_file(ig);
    read_jg_from_file(jg);
    read_di_from_file(di);
    read_b_from_file(F);
    read_n_from_file(N);

    solver s;

    //первый запуск нужен для прогрева, т.е чтобы инициализировать внутренние структуры API
    //(так как API счетчиков производительности сам выделяет память при первом запуске)
    memory_counter::print_memory();

    //второй запуск, соответственно, дает достоверные значения
    memory_counter::print_memory();

    clock_t start = clock();
    std::vector<double> q = s.BCGSTAB(gg, di, F, ig, jg);
    clock_t end = clock();

    std::ofstream out;
    out.open("output_158046.txt");
    if (out.is_open())
    {
        for (const auto value : q)
        {
            out << std::uppercase << std::scientific << value << '\n';
        }
    }
    out.close(); 
    std::cout << "File has been written" << '\n';


    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    printf("The time: %f seconds\n", seconds);
    return 0;
}
