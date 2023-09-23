#ifndef _MODMULT_H
#define _MODMULT_H

#include "bignbr.h"

// These used to be globals. Now calling GetMontgomeryParams* creates them.
struct MontgomeryParams {
  limb MontgomeryMultN[MAX_LEN];
  limb MontgomeryMultR1[MAX_LEN];
  limb MontgomeryMultR2[MAX_LEN];
  int NumberLengthR1;
  // Indicates that the modulus is a power of two.
  int powerOf2Exponent;
};

MontgomeryParams GetMontgomeryParams(int modulus_length, const limb *modulus);
MontgomeryParams GetMontgomeryParamsPowerOf2(int powerOf2,
                                             // computed from the power of 2
                                             int *modulus_length);

// product <- factor1 * factor2 mod modulus
void ModMult(const limb *factor1, const limb *factor2,
             int modulus_length, const limb *modulus,
             limb *product);

void ModPowBaseInt(const MontgomeryParams &params,
                   int modulus_length, const limb *modulus,
                   int base, const limb *exp, int nbrGroupsExp, limb *power);
void ModPow(const MontgomeryParams &params,
            int modulus_length, const limb *modulus,
            const limb *base, const limb *exp, int nbrGroupsExp, limb *power);

void BigIntGeneralModularDivision(const BigInteger *Num, const BigInteger *Den,
                                  const BigInteger *mod, BigInteger *quotient);

// The interface to this is pretty bad: mod and modulus represent the
// same number, and modulus is modified.
void BigIntModularDivision(const MontgomeryParams &params,
                           int modulus_length, limb *modulus,
                           const BigInteger* Num, const BigInteger* Den,
                           const BigInteger* mod, BigInteger* quotient);
void BigIntModularPower(const MontgomeryParams &params,
                        int modulus_length, const limb *modulus,
                        const BigInteger* base, const BigInteger* exponent,
                        BigInteger* power);

void AddBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *TestNbr,
                   int number_length);
void SubtBigNbrModN(const limb *Nbr1, const limb *Nbr2, limb *Sum, const limb *TestNbr,
                    int number_length);

void ComputeInversePower2(const limb *value, /*@out@*/limb *result,
                          int number_length);

#endif