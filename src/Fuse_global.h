#ifndef FUSE_GLOBAL_H
#define FUSE_GLOBAL_H

#ifdef _WIN32
#  define EXPORTIT __declspec( dllexport )
#  define IMPORTIT __declspec( dllimport )
#  define CALL_CONV __cdecl
#else
#  define EXPORTIT
#  define IMPORTIT
#  define CALL_CONV
#endif

#if defined(FUSE_LIBRARY)
#  define FUSESHARED_EXPORT EXPORTIT
#else
#  define FUSESHARED_EXPORT IMPORTIT
#endif

#define R_CALL extern "C" FUSESHARED_EXPORT

#endif // FUSE_GLOBAL_H
