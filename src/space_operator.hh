#include<dune/pdelab/localoperator/idefault.hh>

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>

#include<dune/pdelab/finiteelement/localbasiscache.hh>


/** Lokalni operator za zadaću:
 *
 *   - \Delta u + a*u = f   u \Omega
 *                  u = g   na \Gamma_D\subseteq\partial\Omega
 *  -\nabla u \cdot n = h   na \Gamma_N = \partial\Omega\setminus\Gamma_D
 *
 *
 * \tparam BCType - klasa koja određuje tip rubnog uvjeta.
 */
template <typename BCType, typename FEM>
class DiffusionLOP :
  public Dune::PDELab::NumericalJacobianApplyVolume<DiffusionLOP<BCType,FEM> >,
  public Dune::PDELab::NumericalJacobianVolume<DiffusionLOP<BCType,FEM> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary<DiffusionLOP<BCType,FEM> >,
  public Dune::PDELab::NumericalJacobianBoundary<DiffusionLOP<BCType,FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>         
{
public:

  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };
  enum { doAlphaBoundary = true };

  DiffusionLOP(BCType& bctype_, // boundary cond.type
                         unsigned int intorder_=2) :
    bctype( bctype_ ), intorder( intorder_ )
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    //  lfsu == lfsv
    const int dim = EG::Geometry::coorddimension;
    using Gradient = Dune::FieldVector<double,dim>;

    auto gt = eg.geometry().type();
    const auto & rule = Dune::QuadratureRules<double,dim>::rule(gt,intorder);

    for (const auto & ip : rule)
      {
        // Izračunaj bazne funkcije
        auto& phi = cache.evaluateFunction(ip.position(), lfsu.finiteElement().localBasis());

        // rješenje
        double u=0.0;
        for (size_t i=0; i<lfsu.size(); ++i)
          u += x(lfsu,i)*phi[i];

        // Gradijenti baznih funkcija
        auto& gradphihat = cache.evaluateJacobian(ip.position(), lfsu.finiteElement().localBasis());

        // grad g_K^{-t}
        const auto & jac = eg.geometry().jacobianInverseTransposed(ip.position());
        // Gradijenti baznih funkcija na fizičkom elementu.
        std::vector<Gradient> gradphi(lfsu.size());
        for (size_t i=0; i<lfsu.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i]);

        // Gradijent rješenja u
        Gradient gradu(0.0);
        for (size_t i=0; i<lfsu.size(); ++i)
          gradu.axpy(x(lfsu,i),gradphi[i]);

        // evaluate parameters;
        // auto globalpos = eg.geometry().global(ip.position());
        double f = 0;
        double a = 0;

        // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
        double factor = ip.weight()*eg.geometry().integrationElement(ip.position());
        for (size_t i=0; i<lfsv.size(); ++i)
          r.accumulate(lfsv, i, (gradu*gradphi[i] + a*u*phi[i] - f*phi[i]) * factor);
      }
  }

  // boundary integral
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    // lfsu_s == lfsv_s
    // dimensions
    const int dim = IG::coorddimension;

    Dune::GeometryType gtface = ig.geometryInInside().type();
    const auto & rule = Dune::QuadratureRules<double,dim-1>::rule(gtface,intorder);

    for (auto const & ip : rule)
      {
        // Zanemari Dirichletovu granicu.
        if ( bctype.isDirichlet( ig, ip.position() ) )
          continue;

        // Pozicija integracijske točke u lokalnim koordinatama elementa.
        auto local = ig.geometryInInside().global(ip.position());

        // bazne funkcije u integracijskoj točki
        auto& phi = cache.evaluateFunction(local,lfsu_s.finiteElement().localBasis());
        // rješenje u
        double u=0.0;
        for (size_t i=0; i<lfsu_s.size(); ++i)
          u += x_s(lfsu_s,i)*phi[i];

        // Globalna pozicija integracijske točke
        auto globalpos = ig.geometry().global(ip.position());
        // Neumannov fluks
        double h = ...;

        double factor = ip.weight()*ig.geometry().integrationElement(ip.position());
        for (size_t i=0; i<lfsv_s.size(); ++i)
          r_s.accumulate(lfsv_s,i, h*phi[i]*factor);
      }
  }

private:
  BCType& bctype;
  unsigned int intorder;
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
