#include "Particle.h"
#include "ParticleType.h"
#include "ResonanceType.h"
#include <cmath>
#include <cstdlib> //for RAND_MAX
#include <iostream>

int Particle::fNParticleType;
ParticleType *Particle::fParticleType[fMaxNumParticleType];

int Particle::FindParticle(const char *Name) {
  for (int i = 0; i < fNParticleType; i++) {
    if (Name == fParticleType[i]->GetName()) {
      return i;
    }
  }
  return -1;
}

Particle::Particle(const char *Name, double Px = 0., double Py = 0.,double Pz = 0.)
    : fPx(Px), fPy(Py), fPz(Pz) {
  const int Index = FindParticle(Name);
  if (Index == -1) {
    std::cout << "No correspondece" << '\n';
  }
  fIndex = Index;
}

Particle::Particle() = default; 

int Particle::GetIndex() const { return fIndex; }
double Particle::GetPx() const { return fPx; }
double Particle::GetPy() const { return fPy; }
double Particle::GetPz() const { return fPz; }
double Particle::GetMass() const { return fParticleType[fIndex]->GetMass(); }
int Particle::GetCharge() const { return fParticleType[fIndex]->GetCharge(); }
double Particle::GetNParticleType() { return fNParticleType; }

double Particle::TotalP2() const { return fPx * fPx + fPy * fPy + fPz * fPz; }
double Particle::GetEnergy() const {
  double Energy = sqrt(GetMass() * GetMass() + TotalP2());
  return Energy;
}

double Particle::InvMass(Particle &p) const {
  double EnergySum2 = pow(GetEnergy() + p.GetEnergy(), 2);
  double Pxtot= fPx + p.GetPx();
  double Pytot= fPy + p.GetPy();
  double Pztot= fPz + p.GetPz();
  double PSum2 = Pxtot * Pxtot + Pytot * Pytot + Pztot * Pztot;
  double InvMass = sqrt(EnergySum2 - PSum2);
  return InvMass;
}

void Particle::AddParticleType(const char *Name, double Mass, int Charge,
                               double Width) {

  if (FindParticle(Name) != -1) {
    std::cout << "Already present" << '\n';
  }
  if (fNParticleType == fMaxNumParticleType) {
    std::cout << "Too many particles" << '\n';
  }
  if (Width == 0) {
    //fNParticleType++;
    fParticleType[fNParticleType++] = new ParticleType(Name, Mass, Charge);
  } else {
    //fNParticleType++;
    fParticleType[fNParticleType++] =
        new ResonanceType(Name, Mass, Charge, Width);
  }
}

void Particle::SetIndex(int Index) {
  if (Index < 0 || Index >= fMaxNumParticleType || Index >= fNParticleType) {
    std::cout << "Invalid index" << '\n';
  }
  fIndex = Index;
}
void Particle::SetIndex(const char *Name) { SetIndex(FindParticle(Name)); }

void Particle::SetP(double Px, double Py, double Pz) {
  fPx = Px;
  fPy = Py;
  fPz = Pz;
}

void Particle::PrintParticleType() {
  for (int i = 0; i < fNParticleType; i++) {
    fParticleType[i]->Print();
  }
}

void Particle::Print() const {
  std::cout << "Index: " << fIndex << '\n'
            << "Name: " << fParticleType[fIndex]->GetName() << '\n'
            << "P: (" << fPx << ", " << fPy << ", " << fPz << ")" << '\n';
}


/*-----------------------------------------------------------------------------------------------*/
// parte 2 da online

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetMass() == 0.0) {
    std::cout << "Decayment cannot be preformed if mass is zero" << '\n';
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  if (fIndex > -1) { // add width effect
    // gaussian random numbers

    float x1, x2, w, y1, y2;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;

    massMot += fParticleType[fIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cout << "Decayment cannot be preformed because mass is too low in "
                 "this channel"
              << '\n';
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
            pout * cos(theta));
  dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
            -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz) {
  double energy = GetEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}