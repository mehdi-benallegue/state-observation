#pragma once

// Handle portable symbol export.
// Defining manually which symbol should be exported is required
// under Windows whether MinGW or MSVC is used.
//
// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#  define STATE_OBSERVATION_DLLIMPORT __declspec(dllimport)
#  define STATE_OBSERVATION_DLLEXPORT __declspec(dllexport)
#  define STATE_OBSERVATION_DLLLOCAL
#else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#    define STATE_OBSERVATION_DLLIMPORT __attribute__((visibility("default")))
#    define STATE_OBSERVATION_DLLEXPORT __attribute__((visibility("default")))
#    define STATE_OBSERVATION_DLLLOCAL __attribute__((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#    define STATE_OBSERVATION_DLLIMPORT
#    define STATE_OBSERVATION_DLLEXPORT
#    define STATE_OBSERVATION_DLLLOCAL
#  endif // __GNUC__ >= 4
#endif // defined _WIN32 || defined __CYGWIN__

#ifdef STATE_OBSERVATION_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define STATE_OBSERVATION_DLLAPI
#  define STATE_OBSERVATION_LOCAL
#else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef STATE_OBSERVATION_EXPORTS
#    define STATE_OBSERVATION_DLLAPI STATE_OBSERVATION_DLLEXPORT
#  else
#    define STATE_OBSERVATION_DLLAPI STATE_OBSERVATION_DLLIMPORT
#  endif // STATE_OBSERVATION_EXPORTS
#  define STATE_OBSERVATION_LOCAL STATE_OBSERVATION_DLLLOCAL
#endif // STATE_OBSERVATION_STATIC
