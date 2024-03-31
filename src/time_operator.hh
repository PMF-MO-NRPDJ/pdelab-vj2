#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

/** Vremenski lokalni operator. Ovdje je to samo skalarni produkt.
 *
 * \f{align*}{
                \int_\Omega uv dx
 * \f}
 */
template <typename FEM>
class TimeLOP
  : public Dune::PDELab::NumericalJacobianApplyVolume<TimeLOP<FEM>>,
    public Dune::PDELab::NumericalJacobianVolume<TimeLOP<FEM>>,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  enum { doPatternVolume = true };
  enum { doAlphaVolume = true };

  TimeLOP (unsigned int intorder_=2)
    : intorder(intorder_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    const int dim = EG::Geometry::mydimension;
    // ...
  }
private:
  unsigned int intorder;
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};
