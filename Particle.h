#ifndef PARTICLE_H
#define PARTICLE_H

#include "ParticleType.h"
#include "ResonanceType.h"

class Particle {

public:
  Particle(const char *Name, double Px, double Py, double Pz);
  Particle();

  int GetIndex() const;
  double GetPx() const;
  double GetPy() const;
  double GetPz() const;
  static double GetNParticleType();

  double GetMass() const;
  int GetCharge() const;
  double TotalP2() const;
  double GetEnergy() const;
  double InvMass(Particle &p) const;

  static void AddParticleType(const char *Name, double Mass, int Charge,double Width = 0.);

  void SetIndex(int Index);
  void SetIndex(const char *Name);
  void SetP(double Px, double Py, double Pz);

  static void PrintParticleType();
  void Print() const;

  int Decay2body(Particle &dau1, Particle &dau2) const; // parte 2

private:
  static const int fMaxNumParticleType = 10; // va messo nel public??? NOOOO
  static ParticleType *fParticleType[fMaxNumParticleType];

  static int fNParticleType;

  int fIndex;

  double fPx;
  double fPy;
  double fPz;

  static int FindParticle(const char *Name);

  void Boost(double bx, double by, double bz); // parte 2
};

#endif