
PROJECT(BL)

SET( libBL_SOURCES   sparse_matrix.cpp OPEN_NL_C_API.cpp )

IF(WIN32)
   SET(libBL_SOURCES ${libBL_SOURCES} libBL.h linkage.h oNL_types.h OPEN_NL_C_API.h sparse_matrix.h geex_defs.h )
ENDIF(WIN32)

ADD_LIBRARY( BL SHARED ${libBL_SOURCES} )


