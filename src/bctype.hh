#pragma once

#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/common/function.hh>

#include <cmath>

/*  Selekcija Dirichletove granice */
class DirichletBdry : public Dune::PDELab::DirichletConstraintsParameters
{
    // Dirichletova granica se eventualno može micati u vremenu.
  double time;
public:
  template<typename I>
  bool isDirichlet(const I & intersection,   
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
    auto xg = intersection.geometry().global( coord );
    
    if( xg[0]< 1E-6 )
      return false; // x=0 nije Dirichletova granica
    return true;    // sve ostalo je
  }

  //! Postavi vrijeme. Nije nužno ako Dirichletova granica ne ovisi o vremenu.
//  void setTime (double t) { time = t; }

};


// Klasa obilježja. Ovo je samo pokrata koja čini kod čitljivijim.
template<typename GV>
using ScalarTraits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1>>;

/*  Dirichletov rubni uvjet proširen na čitavu domenu.
 *  U t=0 daje inicijalni uvjet!
 */
template<typename GV>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<ScalarTraits<GV>, BCExtension<GV> >
{
  const GV& gv;
  double  time;
public:
  using Traits = ScalarTraits<GV>;

  // Sačuvaj gridview
  BCExtension (const GV& gv_) : gv(gv_) {}

  inline void evaluate(const typename Traits::ElementType& e,
                       const typename Traits::DomainType& xlocal,
                       typename Traits::RangeType& y) const
  {
    //auto x = e.geometry().global(xlocal);
    y = 0.0;
    return;
  }

  // Referenca na gridview
  inline const GV& getGridView () {return gv;}

  // Postavljanje vremena. PDELab očekuje tu metodu.
  void setTime (double t) {time = t;}
};

// Funkcija koja daje fluks na granici x=0 i ovisi samo o vremenu.
double flux(double time){
  if(time < 0.25)
      return -2*std::sin(2*M_PI*time);
  return 0.0;
}
