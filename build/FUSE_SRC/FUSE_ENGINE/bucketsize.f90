SUBROUTINE BUCKETSIZE()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the maximum water holding capacity of the different reservoirs
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- bucket sizes stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE multiparam                                        ! model parameters
IMPLICIT NONE
! ---------------------------------------------------------------------------------------
! derive maximum tension water in each layer
DPARAM%MAXTENS_1  =     MPARAM%FRACTEN  * MPARAM%MAXWATR_1
DPARAM%MAXTENS_2  =     MPARAM%FRACTEN  * MPARAM%MAXWATR_2
! derive maximum free water in each layer
DPARAM%MAXFREE_1  = (1._sp-MPARAM%FRACTEN) * MPARAM%MAXWATR_1
DPARAM%MAXFREE_2  = (1._sp-MPARAM%FRACTEN) * MPARAM%MAXWATR_2
! derive capacities of the recharge and lower zone (ONLY USED if upper tension is divided in two)
DPARAM%MAXTENS_1A =     MPARAM%FRCHZNE  * DPARAM%MAXTENS_1
DPARAM%MAXTENS_1B = (1._sp-MPARAM%FRCHZNE) * DPARAM%MAXTENS_1
! derive capacities of the primary and secondary parallel baseflow reservoirs
DPARAM%MAXFREE_2A =     MPARAM%FPRIMQB  * DPARAM%MAXFREE_2
DPARAM%MAXFREE_2B = (1._sp-MPARAM%FPRIMQB) * DPARAM%MAXFREE_2
! ---------------------------------------------------------------------------------------
END SUBROUTINE BUCKETSIZE
