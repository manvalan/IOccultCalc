/**
 * @file astdyn_propagation_helper.h
 * @brief Helper functions per propagazione usando ITALOccultLibrary/AstDyn
 */

 #ifndef IOCCULTCALC_ASTDYN_PROPAGATION_HELPER_H
 #define IOCCULTCALC_ASTDYN_PROPAGATION_HELPER_H
 
 #include "astdyn_interface.h"
 #include "orbital_elements.h"
 #include "types.h"
 #include <vector>
 #include <memory>
 #include <tuple>
 
 namespace ioccultcalc {
 
 class AstDynPropagationHelper {
 public:
     static AstDynPropagationHelper& getInstance();
     
     AstDySElements propagate(const AstDySElements& elements, double target_mjd_tdb);
     std::pair<double, double> propagateToRADec(
         const EquinoctialElements& elements,
         const JulianDate& target_jd);
     std::pair<double, double> getRADec(
         const AstDySElements& elements,
         double target_mjd_tdb);
     std::vector<std::tuple<double, double, double>> propagateRange(
         const AstDySElements& elements,
         double start_mjd_tdb,
         double end_mjd_tdb,
         double step_days = 1.0);
     
     static AstDySElements convertFromEquinoctial(const EquinoctialElements& eq);
     static AstDySElements convertFromOrbital(const OrbitalElements& orb);
     
     void setTolerance(double tolerance);
     void setPerturbations(bool planets, bool asteroids, bool relativity);
     
 private:
     AstDynPropagationHelper();
     ~AstDynPropagationHelper() = default;
     AstDynPropagationHelper(const AstDynPropagationHelper&) = delete;
     AstDynPropagationHelper& operator=(const AstDynPropagationHelper&) = delete;
     
     std::unique_ptr<AstDynPropagator> propagator_;
     double tolerance_;
 };
 
 namespace astdyn_utils {
     AstDySElements quickPropagate(const AstDySElements& elements, double target_mjd_tdb);
     std::pair<double, double> quickRADec(const AstDySElements& elements, double target_mjd_tdb);
     AstDySElements toAstDyS(const EquinoctialElements& eq);
     AstDySElements toAstDyS(const OrbitalElements& orb);
 }
 
 } // namespace ioccultcalc
 
 #endif // IOCCULTCALC_ASTDYN_PROPAGATION_HELPER_H