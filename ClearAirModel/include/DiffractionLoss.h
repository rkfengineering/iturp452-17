#ifndef DIFFRACTION_LOSS_H
#define DIFFRACTION_LOSS_H

#include "ProfilePath.h"
#include "Enumerations.h"

namespace DiffractionLoss {

    double delta_bullington(const ProfilePath& path, const double hts, const double hrs, const double hstd, 
            const double hsrd, const double ae, double freqGHz,
            const double omega, const Enumerations::PolarizationType pol);

    struct DiffResults{
        double Ldp;
        double LD50;
    };

    DiffResults diffLoss(const ProfilePath& path, const double hts, const double hrs, const double hstd, 
            const double hsrd, double freqGHz, const double omega, const double p, const double ae, 
            const double ab, const Enumerations::PolarizationType pol);

}//end namespace DiffractionLoss

#endif /* DIFFRACTION_LOSS_H */