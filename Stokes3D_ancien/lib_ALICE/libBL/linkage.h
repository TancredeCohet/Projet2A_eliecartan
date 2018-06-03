#ifndef __LIBNL_LINKAGE__
#define __LIBNL_LINKAGE__

#ifdef WIN32

#ifdef BL_EXPORTS
#define LIBNL_API __declspec( dllexport )
#else
#define LIBNL_API __declspec( dllimport )
#endif

#else

#define LIBNL_API 
#endif

#endif
