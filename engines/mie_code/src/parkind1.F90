MODULE PARKIND1

! +----------------------------------------------------------------------------------------+
! | engines/mie_code/src/parkind1.f90                                                      |
! |                                                                                        |
! | (C) Copyright 2018- ECMWF.                                                             |
! |                                                                                        |
! | This software is licensed under the terms of the Apache Licence Version 2.0            |
! | which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.                   |
! |                                                                                        |
! | In applying this licence, ECMWF does not waive the privileges and immunities           |
! | granted to it by virtue of its status as an intergovernmental organisation             |
! | nor does it submit to any jurisdiction.                                                |
! |                                                                                        |
! | Author:                                                                                |
! |    Alessio Bozzo.       ECMWF                                                          |
! |    Ramiro Checa-Garcia. ECMWF                                                          |
! |                                                                                        |
! | Modifications:                                                                         |
! |    10-May-2018   Alessio Bozzo    1st. version                                         |
! |                                                                                        |
! +----------------------------------------------------------------------------------------+


!
!     *** Define usual kinds for strong typing ***
!
IMPLICIT NONE
SAVE
!
!     Integer Kinds
!     -------------
!
INTEGER, PARAMETER :: JPIT = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: JPIS = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: JPIB = SELECTED_INT_KIND(12)

!Special integer type to be used for sensative adress calculations
!should be *8 for a machine with 8byte adressing for optimum performance
#ifdef ADDRESS64
INTEGER, PARAMETER :: JPIA = JPIB
#else
INTEGER, PARAMETER :: JPIA = JPIM
#endif

!
!     Real Kinds
!     ----------
!
INTEGER, PARAMETER :: JPRT = SELECTED_REAL_KIND(2,1)
INTEGER, PARAMETER :: JPRS = SELECTED_REAL_KIND(4,2)
INTEGER, PARAMETER :: JPRM = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER :: JPRB = SELECTED_REAL_KIND(13,300)
!
END MODULE PARKIND1
