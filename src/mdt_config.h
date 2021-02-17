/** \file mdt_config.h     Utility macros used by all MDT headers.
 *
 *             Part of MDT, Copyright(c) 1989-2021 Andrej Sali
 */

#ifndef __MDT_CONFIG_H
#define __MDT_CONFIG_H

/* Provide macros to mark functions and classes as exported from a DLL/.so */
#ifdef _MSC_VER
  #ifdef MDT_EXPORTS
    #define MDTDLLEXPORT __declspec(dllexport)
  #else
    #define MDTDLLEXPORT __declspec(dllimport)
  #endif
  #define MDTDLLLOCAL
#else
  #ifdef GCC_VISIBILITY
    #define MDTDLLEXPORT __attribute__ ((visibility("default")))
    #define MDTDLLLOCAL __attribute__ ((visibility("hidden")))
  #else
    #define MDTDLLEXPORT
    #define MDTDLLLOCAL
  #endif
#endif

#endif  /* __MDT_CONFIG_H */
