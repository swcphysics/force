! Bingwei 02/03/2011
! - added PC_hbarc_cube
! Bingwei 12/16/2010
! Physical constants
! * The default units are MeV and MeV^(-1), unless noted otherwise.

module nneft_phyconst

  use nneft_type
  implicit none

  ! PC = Physical Constant
  ! PC_default_* : the real-world values of various constant
  ! Essential constants
  real(NER), parameter  ::                  &
  & PC_default_gA = 1.29_NER,               &
  & PC_default_mN = 939.0_NER ,             &
  & PC_default_fpi = 92.4_NER ,             &
  & PC_default_mpi = 138.0_NER ,            &
  & PC_default_hbarc = 197.3286_NER

  ! Values of ci's
  real(NER), parameter  ::                  &
  ! Epelbaum et al, 1999
  ! & PC_Epel99_c1 = -0.00081_NER,          &
  ! & PC_Epel99_c2 = 0.00328_NER,           &
  ! & PC_Epel99_c3 = -0.00470_NER,          &
  ! & PC_Epel99_c4 = 0.0034_NER,            &
  & PC_Epel99_c1 = -0.81E-3_NER,            &
  & PC_Epel99_c2 = 3.28E-3_NER,             &
  & PC_Epel99_c3 = -4.7E-3_NER,             &
  & PC_Epel99_c4 = 3.4E-3_NER,              &
  ! Krebs' Q^2 delta-less fit
  ! & PC_Kreb_c1 = -0.00057_NER,          &
  ! & PC_Kreb_c2 = 0.00284_NER,           &
  ! & PC_Kreb_c3 = -0.00387_NER,          &
  ! & PC_Kreb_c4 = 0.00289_NER,           &
  & PC_Kreb_c1 = -0.57E-3_NER,              &
  & PC_Kreb_c2 = 2.84E-3_NER,               &
  & PC_Kreb_c3 = -3.87E-3_NER,              &
  & PC_Kreb_c4 = 2.89E-3_NER,               &
  ! Valderrama
  ! & PC_Pavon_c1 = -0.00081_NER,             &
  ! & PC_Pavon_c2 = 0.00328_NER,              &
  ! & PC_Pavon_c3 = -0.0034_NER,              &
  ! & PC_Pavon_c4 = 0.0034_NER,               &
  & PC_Pavon_c1 = -0.81E-3_NER,             &
  & PC_Pavon_c2 = 3.28E-3_NER,              &
  & PC_Pavon_c3 = -3.4E-3_NER,              &
  & PC_Pavon_c4 = 3.4E-3_NER,               &
  ! Roy Steiner equation, Table II of PRC 96, 024004
  & PC_RSE2_c1 = -0.74E-3_NER,               &
  & PC_RSE2_c3 = -3.61E-3_NER,               &
  & PC_RSE2_c4 = 2.44E-3_NER,                &
  & PC_RSE3_c1 = -1.07E-3_NER,               &
  & PC_RSE3_c3 = -5.32E-3_NER,               &
  & PC_RSE3_c4 = 3.56E-3_NER,                &
  & PC_RSE4_c1 = -1.10E-3_NER,               &
  & PC_RSE4_c3 = -5.54E-3_NER,               &
  & PC_RSE4_c4 = 4.17E-3_NER,                &
  & PC_default_c1 = PC_Pavon_c1,            &
  & PC_default_c2 = PC_Pavon_c2,            &
  & PC_default_c3 = PC_Pavon_c3,            &
  & PC_default_c4 = PC_Pavon_c4

  ! Delta-isobar related parameters
  real(NER), parameter  ::                  &
  & PC_default_c1_del = -0.00057_NER ,      &
  & PC_default_c2_del = -0.00192_NER ,      &
  & PC_default_c3_del = -0.00088_NER ,      &
  & PC_default_c4_del =  0.0005_NER ,       &
  & PC_default_hA = 1.75_NER,               &
  & PC_default_b3 = 0.0_NER,                &
  & PC_default_b8 = 0.0_NER,                &
  & PC_default_cr = 1.0_NER,                &
  & PC_default_delta = 293.0_NER

  ! Misc
  real(NER), parameter  ::                  &
  & PC_default_Lambda_inner = 1000.0_NER

  real(NER), parameter   ::                               &
  & PC_default_gA_sqr = PC_default_gA * PC_default_gA,             &
  & PC_default_mpi_sqr = PC_default_mpi * PC_default_mpi,          &
  & PC_default_fpi_sqr = PC_default_fpi * PC_default_fpi,          &
  & PC_default_gA_4th = PC_default_gA_sqr * PC_default_gA_sqr,     &
  & PC_default_mpi_4th = PC_default_mpi_sqr * PC_default_mpi_sqr,  &
  & PC_default_fpi_4th = PC_default_fpi_sqr * PC_default_fpi_sqr,  &
  & PC_default_hbarc_cube = PC_default_hbarc * PC_default_hbarc * PC_default_hbarc

  ! Variables for physical constants
  real(NER) ::                    &
  & PC_gA = PC_default_gA,           &
  & PC_mN = PC_default_mN,           &
  & PC_fpi = PC_default_fpi,         &
  & PC_mpi = PC_default_mpi,         &
  & PC_hbarc = PC_default_hbarc,     &
  & PC_c1_del = PC_default_c1_del,   &
  & PC_c2_del = PC_default_c2_del,   &
  & PC_c3_del = PC_default_c3_del,   &
  & PC_c4_del = PC_default_c4_del,   &
  & PC_c1 = PC_default_c1,           &
  & PC_c2 = PC_default_c2,           &
  & PC_c3 = PC_default_c3,           &
  & PC_c4 = PC_default_c4,           &
  & PC_hA = PC_default_hA,           &
  & PC_b3 = PC_default_b3,           &
  & PC_b8 = PC_default_b8,           &
  & PC_cr = PC_default_cr,           &
  & PC_delta = PC_default_delta,     &
  & PC_Lambda_inner = PC_default_Lambda_inner

  ! real(NER) ::                    &
  ! & PC_gA_sqr = PC_default_gA_sqr,   &
  ! & PC_mpi_sqr = PC_default_mpi_sqr, &
  ! & PC_fpi_sqr = PC_default_fpi_sqr, &
  ! & PC_gA_4th = PC_default_gA_4th,   &
  ! & PC_mpi_4th = PC_default_mpi_4th, &
  ! & PC_fpi_4th = PC_default_fpi_4th, &
  ! & PC_hbarc_cube = PC_default_hbarc_cube

  type :: struct_pc

    real(NER) ::                    &
    & gA = PC_default_gA,           &
    & mN = PC_default_mN,           &
    & fpi = PC_default_fpi,         &
    & mpi = PC_default_mpi,         &
    & hbarc = PC_default_hbarc,     &
    & c1_del = PC_default_c1_del,   &
    & c2_del = PC_default_c2_del,   &
    & c3_del = PC_default_c3_del,   &
    & c4_del = PC_default_c4_del,   &
    & c1 = PC_default_c1,           &
    & c2 = PC_default_c2,           &
    & c3 = PC_default_c3,           &
    & c4 = PC_default_c4,           &
    & hA = PC_default_hA,           &
    & b3 = PC_default_b3,           &
    & b8 = PC_default_b8,           &
    & cr = PC_default_cr,           &
    & delta = PC_default_delta,     &
    & Lambda_inner = PC_default_Lambda_inner

  end type struct_pc

  type  :: strct_cis

    real(NER) ::            &
    & c1 = PC_default_c1,   &
    & c3 = PC_default_c3,   &
    & c4 = PC_default_c4

  end type strct_cis

  ! type (strct_cis)  ::

end module nneft_phyconst
