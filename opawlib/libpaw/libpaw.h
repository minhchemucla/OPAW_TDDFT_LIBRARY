/* libpaw.h */

/*
 * This file is part of the libPAW library.
 * It has to be customized according to the host code.
 * For the time being there are 2 known host codes:
 * ABINIT (www.abinit.org) and BigDFT (bigdft.org).
 */

/*
 * Copyright (C) 2014-2018 ABINIT Group (MT)
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 */

/* config.h should contain all preprocessing directives */
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/*Note!!!*/
/*#define HAVE_LIBXC */
/*#define HAVE_FC_ISO_C_BINDING */
/*#define LIBPAW_HAVE_LIBXC */
/*#define LIBPAW_ISO_C_BINDING */
/* =============================
 * ========= DEFAULT ===========
 * ============================= */

/* Constants and defs */
#  define USE_DEFS use m_libpaw_defs

/* MPI wrappers */
#  define USE_MPI_WRAPPERS use m_libpaw_mpi

/* Messages, errors */
#  define USE_MSG_HANDLING use m_libpaw_tools, only : wrtout => libpaw_wrtout, libpaw_msg_hndl
#  define MSG_COMMENT(msg) call libpaw_msg_hndl(msg,"COMMENT","PERS")
#  define MSG_WARNING(msg) call libpaw_msg_hndl(msg,"WARNING","PERS")
#  define MSG_ERROR(msg)   call libpaw_msg_hndl(msg,"ERROR"  ,"PERS")
#  define MSG_BUG(msg)     call libpaw_msg_hndl(msg,"BUG"    ,"PERS")
#  undef  HAVE_YAML

/* Allocation/deallocation */
#  define USE_MEMORY_PROFILING
/* Use this to allocate/deallocate basic-type arrays with sizes */
#  define LIBPAW_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate/deallocate basic-type pointers with sizes */
#  define LIBPAW_POINTER_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_POINTER_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate/deallocate user-defined-type arrays with sizes */
#  define LIBPAW_DATATYPE_ALLOCATE(ARR,SIZE) allocate(ARR SIZE)
#  define LIBPAW_DATATYPE_DEALLOCATE(ARR) deallocate(ARR)
/* Use this to allocate user-defined-type arrays with explicit bounds */
#  define LIBPAW_BOUND1_ALLOCATE(ARR,BND1) allocate(ARR(BND1))
#  define LIBPAW_BOUND2_ALLOCATE(ARR,BND1,BND2) allocate(ARR(BND1,BND2))
#  define BOUNDS(LBND,UBND) LBND : UBND

/* libXC support */
/*#  undef LIBPAW_HAVE_LIBXC*/

/* Netcdf support */
#  undef LIBPAW_HAVE_NETCDF

/* FoX support */
#  undef LIBPAW_HAVE_FOX

/* F2008 support */
#  define LIBPAW_CONTIGUOUS

/* =============================
 * =========== END =============
 * ============================= */



/* =============================
 * ===== COMMON DEFINITIONS ====
 * ============================= */

/* Error handlers for netcdf; libpaw_netcdf_check is defined in m_libpaw_tools */
#ifndef NCF_CHECK
#  define NCF_CHECK(ncerr) if (ncerr/=nf90_noerr) call libpaw_netcdf_check(ncerr, "ncf_check")
#endif

