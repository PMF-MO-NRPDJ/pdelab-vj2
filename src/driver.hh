#pragma once

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include <dune/pdelab/constraints/conforming.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>

#include "bctype.hh"
#include "space_operator.hh"
#include "time_operator.hh"


/** Upravljačka rutina koja koordinira sav posao osim konstrukcije
 *  mreže. 
 *  @tparam GV = Leaf grid view tip 
 *  @param gv = leaf grid view 
 *  @param dt = vremenski korak
 *  @param tend = krajnje vrijeme simulacije
 *  */
template<typename GV>
void driver(const GV& gv, double dt, double tend)
{
    // 1. Izbor konačnog elementa.
    const int k=1;
    using FEM =  Dune::PDELab::QkLocalFiniteElementMap<GV, double, double, k>;
    FEM fem(gv);
    // 2. Tip ograničenja
    using CON = Dune::PDELab::ConformingDirichletConstraints;
    // 3. Tip vektora i matrica
    using VBE = Dune::PDELab::ISTL::VectorBackend<>;
    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    // 4. Prostor konačnih elemenata
    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
    GFS gfs(gv,fem);
    // 5. Konstruiraj spremnika za ograničenja u prostoru konačnih elemenata.
    using  CC = typename GFS::template ConstraintsContainer<double>::Type;
    CC cc;
    // 6.  Asembliranje Dirichleteovih ograničenja.
    DirichletBdry bctype;     // Korisnička klasa koja detektora Dirichletov rubni uvjet.
//    bctype.setTime(0.0);    // Ako b.c. ovisi o vremenu postavi vrijeme, time = 0.0
    Dune::PDELab::constraints(bctype, gfs, cc); // Definiraj Dirichletove vrhove
    // 7. Konstrukcija mrežnog operatora.
    // Lokalni operator za prostorni dio
    using SLOP = ...
    SLOP ...;        // prostorni lokalni operator
    // Lokalni operator za vremenski dio
    using TLOP = ...;
    TLOP ...;               // vremenski lokalni operator
    MBE mbe(9);
    // Grid operatori koji odgovaraju vremenskom i prostornom dijelu
    using GO0 = ....;
    GO0 go0(...);  // prostorni GO
    using GO1 = ...;
    GO1 go1(...);  // vremenski GO
    // Potpuni grid operator se dobiva jednokoračnom vremenskom metodom
    using IGO = ... ;
    IGO igo(...);           // grid operator jednokoračne metode
    // 8. Konstruirati vektor stupnjeva slobode
    using U =  Dune::PDELab::Backend::Vector<GFS, double>;
    U uold(gfs,0.0);             // rješenje u t=t^n, početni uvjet je 0
    // Vektor stupnjeva slobode može se dobiti i na ovaj način:
    //  using U = typename IGO::Traits::Domain;
    // 9. Instancirati klasu koja daje Dirichletov rubni uvjet i ubaciti Dirichletove vrijednosti
    //    u vektor rješenja (stupnjeva slobode).
    using G = BCExtension<GV>;    // rubni uvjet i početni uvjet
    G bcond(gv);                 // rubni uvjet sada ovisi o vremenu
    bcond.setTime(0.0);          // time = 0.0
    Dune::PDELab::interpolate(bcond,gfs,uold);  // početni uvjet je dan ovdje
    // 10. Izbor linearnog rješavača i prekondicionera.
    using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR;
    LS ls(5000,false);
    // 11. Instanciranje linearnog ili nelinearnog rješavača i rješavanje problema
    // Zadaća je linearna
    using PDESOLVER = Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,U>;
    PDESOLVER pdesolver(igo,ls,1e-10);
    // 12. Izbor parametara vremenske diskretizacije
    using TDM = ...;
    TDM method;                      // koeficijenti vremenske diskretizacije
    // Druge mogućnosti:
    //  using TDM = Dune::PDELab::ImplicitEulerParameter<double>; // implicitni Euler
    //  using TDM = Dune::PDELab::OneStepThetaParameter<double>; // theta metoda, theta parametar uzima konstruktor
    //  using TDM = Dune::PDELab::Alexander2Parameter<double>;
    // 13. Tip vremenske diskretizacije - jednokoračna metoda
    using OSM = ...;
    OSM osm(...);   // jednokoračna metoda za rješenje sustava
    osm.setVerbosityLevel(1);

    //14. Grafički prikaz - kreiranje diskretne mrežne funkcije
    using DGF = Dune::PDELab::DiscreteGridFunction<GFS,U>;
    DGF udgf(gfs,uold);
    // Za Q2 elemente trebamo koristiti subsampling.
    using VTKW = Dune::SubsamplingVTKWriter<GV>;
    VTKW vtkwriter(gv, Dune::RefinementIntervals{k});
    vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf,"solution"));
    // Za ispisivanje vremenskog niza podataka koristimo VTKSequenceWriter<>
    Dune::VTKSequenceWriter<GV> writer(std::make_shared<VTKW>(vtkwriter), "out");
    writer.write(0.0); // ispiši podatke u trenutku t=0.

    // 15. Vremenska petlja
    double time = 0.0;      // početni trenutak
    U unew(gfs,0.0);        // sljedeći vremenski sloj -jednokoračna metoda
    while (time < tend) {
        // postavi novo vrijeme u BC klasu -- za Dirichletovu granicu koja se mijenja u vremenu
//        bctype.setTime(time+dt);                   // izračunaj Dirichletov r.u
//        cc.clear();                                // u ovom vremenskom trenutku
//        Dune::PDELab::constraints(bctype,gfs,cc);

         // ....

        
        // Kontrola vremenskog koraka ovisno o broju linearnih iteracija.
        // Parametre namjestiti eksperimentalno
        //  if(noIter < 50) dt *= 1.2;
        //  if(noIter > 1000) dt /= 1.2;
    }
}

