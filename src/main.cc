/** \file

    \brief Parabolička jednadžba diskretizirana konformnom metodom KE u prostoru
          i implicitnom Eulerovom metodom u vremenu.
    */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include "driver.hh"

// Upotreba:
// main [n [tend [dt]]]
//    n = broj profinjenja (0 default)
//    dt = vremenski korak (0.1 default)
//    tend = krajnje vrijeme simulacije (2.0 default)
int main(int argc, char **argv)
{
     Dune::MPIHelper::instance(argc, argv);
     // Konstrukcija mreže
     const int dim = 2;
     Dune::FieldVector<double, dim> L(1.0);
     std::array<int, dim> N={10,10};
     Dune::YaspGrid<dim> grid(L, N);
     int level = 0;
     if(argc>1)
         level = std::stoi(argv[1]);
     grid.globalRefine(level);
     using GV = Dune::YaspGrid<dim>::LeafGridView;
     const GV & gv = grid.leafGridView();
     double tend = 2.0;
     double dt = 0.1;
     if(argc > 2)
         tend = std::stod(argv[2]);
     if(argc > 3)
         dt = std::stod(argv[3]);

     driver(gv, dt, tend);
     return 0;
}
