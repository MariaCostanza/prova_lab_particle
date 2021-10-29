#include "ResonanceType.h"
#include "ParticleType.h"

#include <iostream>

ResonanceType::ResonanceType(const char *Name, double Mass, int Charge, double Width)
    : ParticleType (Name, Mass, Charge), fWidth(Width) {}

double ResonanceType::GetWidth() const { return fWidth; }

void ResonanceType::Print() const{
    ParticleType::Print();
    std::cout << "Width: "<<fWidth<<'\n';
}