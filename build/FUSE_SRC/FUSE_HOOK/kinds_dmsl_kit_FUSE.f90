!******************************************************************
! (C) Copyright 2000-2010  ---  Dmitri Kavetski  ---  All rights reserved
! NB: CUSTOMIZED VERSION FOR FUSE SUITE OF MARTYN CLARK
!******************************************************************
module kinds_dmsl_kit_FUSE
! Purpose: a. Defines global numeric kinds for DMSL;
!          b. Contains machine precision information;
!          c. Contains global information for DMSL library support.
! This module is typically made globally available.
! ---
! Programmer: Dmitri Kavetski.
! 2000 - 2004
! Civil, Environmental Engineering and Surveying
! University of Newcastle, Callaghan, NSW 2308, Australia.
! 2004 - 2007
! Department of Civil and Environmental Engineering
! Princeton University, Princeton, NJ 08544, USA.
! 2007 - current
! Civil, Environmental Engineering and Surveying
! University of Newcastle, Callaghan, NSW 2308, Australia.
! ---
! Comments:
! 1. The log[] function (2b) may not compile on some compilers.
! 2. The complex constants may not compile on some compilers.
! 3. If the compiler prevents the direct definitions below, hardcode them
!    with at least 40 significant digits of precision.
implicit none
public
! ---
! (1) Parameterised numeric data types
! (a) Available precision (CVF reals: 4=single, 8=double, 16=quad; CVF integers: 4=short, 8=long)
integer,     parameter::srk=selected_real_kind(p=4)   ! single precision
integer,     parameter::drk=selected_real_kind(p=8)   ! double precision
integer,     parameter::qrk=selected_real_kind(p=16)  ! quadruple precision
integer,     parameter::sik=selected_int_kind(r=4)    ! short integer
integer,     parameter::lik=selected_int_kind(r=8)    ! long integer (NB: NR-90 uses r=9)
! (b) Selected global precision in all DMSL units
integer,     parameter::mrk=drk                       ! global real kind
integer,     parameter::mik=lik                       ! global integer kind
real(mrk),   parameter::protoRe=1._mrk                ! prototype of real(mrk) number
integer(mik),parameter::protoInt=1_mik                ! prototype of integer(mik) number
integer,     parameter::mck=kind((1._mrk,1._mrk))     ! global complex kind
integer,     parameter::mlk=kind(.true.)              ! global logical kind
complex(mck),parameter::protoCmx=((1._mrk,1._mrk))    ! prototype of complex(mck) number
logical(mlk),parameter::protoLog=.true.               ! prototype of logical(mlk) number
! (c) Compiler-specific info [best kept up to date, I guess ...]
integer(mik),parameter::mrkBy=mrk                     ! number of bytes to store protoRe
integer(mik),parameter::mikBy=mik                     ! number of bytes to store protoInt
integer(mik),parameter::mckBy=2*mrk                   ! number of bytes to store protoCmx
integer(mik),parameter::mlkBy=4                       ! number of bytes to store protoLog
! NB:
! On CVF compiler: mrk and mik also denote the number of bytes used to store the value,
!                  mlk requires 4 bytes storage
! single precision      =  32-bit   (4 bytes)
! double precision      =  64-bit   (8 bytes)
! quadruple precisison  = 128-bits (16 bytes)
! ---
! (2) Machine precision information
! (a) Intrinsix
real(mrk),   parameter::tinyRe=tiny(protoRe)          ! smallest real on machine
real(mrk),   parameter::epsRe= epsilon(protoRe)       ! normalised machine accuracy
real(mrk),   parameter::hugeRe=huge(protoRe)          ! largest real on machine
integer(mik),parameter::hugeInt= huge(protoInt)       ! largest integer on machine
real(mrk),   parameter::hugeIntR=real(hugeInt,mrk)    ! largest integer (real format)
! real(mrk),   parameter::hugeIntR=2.14748364700000E+009_mrk      ! Salford Software FTN95
! complex(mck),parameter::tinyC=(tinyRe,tinyRe)         ! smallest complex on machine
! complex(mck),parameter::epsC= (epsRe,epsRe)           ! complex machine precision
! complex(mck),parameter::hugeC=(hugeRe,hugeRe)         ! largest complex on machine
! (b) Functions of machine precision
integer(mik),parameter::minExpRei=minexponent(protoRe)  ! min exponent (int)  in machine base (usually radix=2)
real(mrk),   parameter::minExpRer=real(minExpRei,mrk)   ! min exponent (real) in machine base (usually radix=2)
! real(mrk),   parameter::minExpRer=-1.02100000000000E+003_mrk    ! Salford Software FTN95
integer(mik),parameter::maxExpRei=maxexponent(protoRe)  ! max exponent (int)  in machine base (usually radix=2)
real(mrk),   parameter::maxExpRer=real(maxExpRei,mrk)   ! max exponent (real) in machine base (usually radix=2)
! real(mrk),   parameter::maxExpRer=+1.02400000000000E+003_mrk    ! Salford Software FTN95
real(mrk),   parameter::radixRer=real(radix(protoRe),mrk)       ! radix expressed as real
! real(mrk),   parameter::radixRer=2.00000000000000E+000_mrk      ! Salford Software FTN95
! real(mrk),   parameter::nDecDigitsRe=-log10(epsRe)            ! number of decimal digits
!real(mrk),   parameter::nDecDigitsRe=-log(epsRe)/log(10._mrk)   ! number of decimal digits
! real(mrk),   parameter::nDecDigitsRe=1.56535597745270E+001_mrk  ! Salford Software FTN95
!real(mrk),   parameter::lnEpsRe=log(epsRe)                      ! ln[] of machine precision
real(mrk),   parameter::lnEpsRe=3.60436533891172E+001_mrk       ! Salford Software FTN95
!real(mrk),   parameter::lunflw=minExpRer*log(radixRer)                ! =log(tinyRe)  ! ln[] of smallest real
! real(mrk),   parameter::lunflw=-7.07703271351704E+002_mrk       ! Salford Software FTN95
!real(mrk),   parameter::lovflw=(1._mrk-epsRe)*maxExpRer*log(radixRer) ! =log(hugeRe)  ! ln[] of largest real, safe to exponentiate
! real(mrk),   parameter::lovflw=+7.09782712893384E+002_mrk     ! Salford Software FTN95
! ---
! (3) Parameterised machine settings
integer(mik),parameter::keyboardUnit=5                ! keyboard unit (default input)
integer(mik),parameter::screenUnit=6                  ! screen unit   (default output)
! ---
! (4) Library support features
integer(mik),parameter::DMSL_vernum=417
character(*),parameter::DMSL_authorName="Dmitri Kavetski"
character(*),parameter::DMSL_authorEmail="dmitri.kavetski@newcastle.edu.au"
! ---
! (5) Special DMSL values, conventionally used to flag un-initialised variables
real(mrk),   parameter::undefRN=-999999999._mrk       ! flag for undefined real numbers
real(mrk),   parameter::undefRNH=-0.5_mrk*hugeRe      ! huge flag for undefined real numbers
integer(mik),parameter::undefIN=-999999999            ! flag for undefined integer numbers
integer(mik),parameter::undefINH=-hugeInt/2     ! huge flag for undefined integer numbers
logical(mlk),parameter::undefLG=.false.               ! flag for undefined logicals
! complex(mck),parameter::undefCZ=(undefRN,undefRN)    ! flag for undefined complex numbers
complex(mck),parameter::undefCZ=(-999999999._mrk,-999999999._mrk) ! flag for undefined complex numbers
! complex(mck),parameter::undefCZH=(undefRNH,undefRNH)  ! huge flag for undefined complex numbers
! complex(mck),parameter::undefCZH=cmplx(undefRNH,undefRNH,kind=mck)  ! huge flag for undefined complex numbers
character(*),parameter::undefCH="undefined"          ! flag for undefined character strings
! ---
! (6) DMSL-wide registered settings
integer(mik),parameter::iyes=1,ino=0 ! integer flags for true/false
! ---
endmodule kinds_dmsl_kit_FUSE
!******************************************************************
! module makeKinds_dmsl_kit
! implicit none
! contains
! !----------------------------------------------------
! !----------------------------------------------------
! endmodule makeKinds_dmsl_kit
!******************************************************************
