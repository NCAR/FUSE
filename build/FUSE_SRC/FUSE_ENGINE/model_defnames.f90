module model_defnames
! Purpose: Contains routines for alternating between char <-> int names
! Programmers: David McInerney and Dmitri Kavetski (University of Adelaide)
USE nrtype 
implicit none
! parameterised descriptions
integer(I4B), parameter ::  iopt_additive_e = 1001, & 
                            iopt_multiplc_e = 1002, &
                            iopt_tension1_1 = 2001, &
                            iopt_tension2_1 = 2002, &
                            iopt_onestate_1 = 2003, & 
                            iopt_tens2pll_2 = 3001, &
                            iopt_unlimfrc_2 = 3002, & 
                            iopt_unlimpow_2 = 3003, &
                            iopt_fixedsiz_2 = 3004, &
                            iopt_topmdexp_2 = 3005, &
                            iopt_arno_x_vic = 4001, &
                            iopt_prms_varnt = 4002, &
                            iopt_tmdl_param = 4003, &
                            iopt_perc_f2sat = 5001, &
                            iopt_perc_w2sat = 5002, &
                            iopt_perc_lower = 5003, &
                            iopt_sequential = 6001, &
                            iopt_rootweight = 6002, &
                            iopt_intflwnone = 7001, &
                            iopt_intflwsome = 7002, &
                            iopt_rout_gamma = 8001, &
                            iopt_no_routing = 8002, &
                            iopt_no_snowmod = 8501, &
                            iopt_temp_index = 8502
! ---
integer(I4B), parameter ::  iopt_TENS1A = 9001, &
                            iopt_TENS1B = 9002, &
                            iopt_TENS_1 = 9003, &
                            iopt_FREE_1 = 9004, &
                            iopt_WATR_1 = 9005, &
                            iopt_TENS_2 = 9006, &
                            iopt_FREE2A = 9007, &
                            iopt_FREE2B = 9008, &
                            iopt_WATR_2 = 9009
! ------------------------------------------
contains
! ------------------------------------------
elemental function desc_str2int(name)result(res)
! Purpose: Converts a string description into its corresponding integer value.
implicit none
! dummies
character(*), intent(in) :: name
integer(I4B) :: res
! Start procedure here
selectcase(name)
case("additive_e"); res = iopt_additive_e
case("multiplc_e"); res = iopt_multiplc_e
case("tension1_1"); res = iopt_tension1_1
case("tension2_1"); res = iopt_tension2_1
case("onestate_1"); res = iopt_onestate_1
case("tens2pll_2"); res = iopt_tens2pll_2
case("unlimfrc_2"); res = iopt_unlimfrc_2
case("unlimpow_2"); res = iopt_unlimpow_2
case("fixedsiz_2"); res = iopt_fixedsiz_2
case("arno_x_vic"); res = iopt_arno_x_vic
case("prms_varnt"); res = iopt_prms_varnt
case("tmdl_param"); res = iopt_tmdl_param
case("perc_f2sat"); res = iopt_perc_f2sat
case("perc_w2sat"); res = iopt_perc_w2sat
case("perc_lower"); res = iopt_perc_lower
case("sequential"); res = iopt_sequential
case("rootweight"); res = iopt_rootweight
case("intflwnone"); res = iopt_intflwnone
case("intflwsome"); res = iopt_intflwsome
case("rout_gamma"); res = iopt_rout_gamma
case("no_routing"); res = iopt_no_routing
case("no_snowmod"); res = iopt_no_snowmod
case("temp_index"); res = iopt_temp_index
case("TENS1B");  res = iopt_TENS1B
case("TENS_1");  res = iopt_TENS_1
case("FREE_1");  res = iopt_FREE_1
case("WATR_1");  res = iopt_WATR_1
case("TENS_2");  res = iopt_TENS_2
case("FREE2A");  res = iopt_FREE2A
case("FREE2B");  res = iopt_FREE2B
case("WATR_2");  res = iopt_WATR_2
case default;    res = -999
endselect
! End procedure here
endfunction desc_str2int
! ------------------------------------------
elemental function desc_int2str(intVal)result(res)
! Purpose: Converts an integer description into corresponding string value
implicit none
! dummies
integer(I4B), intent(in) :: intVal
character(10) :: res
! Start procedure here
selectcase(intVal)
case(iopt_TENS1B); res = "TENS1B"
case(iopt_TENS_1); res = "TENS_1"
case(iopt_FREE_1); res = "FREE_1"
case(iopt_WATR_1); res = "WATR_1"
case(iopt_TENS_2); res = "TENS_2"
case(iopt_FREE2A); res = "FREE2A"
case(iopt_FREE2B); res = "FREE2B"
case(iopt_WATR_2); res = "WATR_2"
case(iopt_additive_e); res = "additive_e" 
case(iopt_multiplc_e); res = "multiplc_e"
case(iopt_tension1_1); res = "tension1_1" 
case(iopt_tension2_1); res = "tension2_1"
case(iopt_onestate_1); res = "onestate_1"
case(iopt_tens2pll_2); res = "tens2pll_2"
case(iopt_unlimfrc_2); res = "unlimfrc_2"
case(iopt_unlimpow_2); res = "unlimpow_2"
case(iopt_fixedsiz_2); res = "fixedsiz_2"
case(iopt_arno_x_vic); res = "arno_x_vic"
case(iopt_prms_varnt); res = "prms_varnt"
case(iopt_tmdl_param); res = "tmdl_param"
case(iopt_perc_f2sat); res = "perc_f2sat"
case(iopt_perc_w2sat); res = "perc_w2sat"
case(iopt_perc_lower); res = "perc_lower"
case(iopt_sequential); res = "sequential"
case(iopt_rootweight); res = "rootweight"
case(iopt_intflwnone); res = "intflwnone"
case(iopt_intflwsome); res = "intflwsome"
case(iopt_rout_gamma); res = "rout_gamma"
case(iopt_no_routing); res = "no_routing"
case(iopt_no_snowmod); res = "no_snowmod"
case(iopt_temp_index); res = "temp_index"
case default;      res = "UNDFND"
endselect
! End procedure here
endfunction desc_int2str
! ------------------------------------------
endmodule model_defnames
