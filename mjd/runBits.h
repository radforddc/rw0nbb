/*
 *  runBits.h -- define run bit descriptions
 *  David Radford   Nov 2016
 */

#ifndef _RUN_BITS_H
#define _RUN_BITS_H


char runBitDesc[32][128] =
  {
    "DAQExpertMode",
    "bb-decay: 'good' data taking.",
    "Calibration-Prototype: Calibration source present on Prototype.",
    "Calibration-Module 1: Calibration source present on Module 1.",
    
    "Calibration-Module 2: Calibration source present on Module 2.",
    "Co60: Use of Cobalt-60 source.",
    "Th228: Use of Th-228 source.",
    "Partial shield: Part of poly shield is not present.",
    
    "Prototype offline: Prototype is offline.",
    "Module 1 offline: Module 1 is offline.",
    "Module 2 offline: Module 2 is offline.",
    "Radon purge offline: Radon purge is offline or abnormal.",
    
    "Machine shop:  Work ongoing in machine shop.",
    "Disruptive work: Work in Detector Room might be disruptive to data taking",
    "Blank Monolith (East side): Blank monolith is in east shield spot.",
    "Blank Monolith (South side): Blank monolith is in south shield spot.",
    
    "Transition (Source): Moving a source in or out a module.",
    "Nonstandard: Unverified change has been made to the software and/or hardware.",
    "Pulser-Cal: Pulser calibration runs.",
    "Electronics-Cal: Electronics calibrations runs.",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    ""
  };

#endif /*#ifndef _RUN_BITS_H */
