
#ifndef CELLPOLARITYODESYSTEM_HPP_
#define CELLPOLARITYODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Delta-Notch ODE system described by Collier et al,
 * "Pattern formation by lateral inhibition with feedback: a mathematical
 * model of delta-notch intercellular signalling" (Journal of Theoretical
 * Biology 183:429-446, 1996).
 */
class CellPolarityOdeSystem : public AbstractOdeSystem
{
private:

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

    ///\todo extract model parameters as member variables

public:

    /**
     * Default constructor.
     *
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    CellPolarityOdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~CellPolarityOdeSystem();

    /**
     * Compute the RHS of the  Nissen et al. polarity ODE system
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using  Collier et al. system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellPolarityOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a DeltaNotchOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const CellPolarityOdeSystem * t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a DeltaNotchOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, CellPolarityOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)CellPolarityOdeSystem(state_variables);
}
}
} // namespace ...

#endif /*CELLPOLARITYODESYSTEM_HPP_*/
