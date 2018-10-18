module V_retarGV0_pwd
	use V_retar_pwd
	implicit none

 contains


  function V_retarGV0_j0j(j,p_,p,E)
      real(NER) :: V_retarGV0_j0j
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
 

     if(j<0) then
      V_retarGV0_j0j=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*V_retar_j0j(j,p_,k)*OPE_j0j(j,k,p)*wghtslg_pwo(ii)

       end do

       V_retarGV0_j0j=regsum
    
  end function V_retarGV0_j0j



  function V_retarGV0_j1j(j,p_,p,E)
      real(NER) :: V_retarGV0_j1j
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
  

     if(j<0) then
      V_retarGV0_j1j=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*V_retar_j1j(j,p_,k)*OPE_j1j(j,k,p)*wghtslg_pwo(ii)

       end do

       V_retarGV0_j1j=regsum

  end function V_retarGV0_j1j



    function V_retarGV0_jpp(j,p_,p,E)
      real(NER) :: V_retarGV0_jpp
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     

     if(j<0) then
      V_retarGV0_jpp=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*V_retar_jpp(j,p_,k)*OPE_jpp(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*V_retar_jmp(j,p_,k)*OPE_jpm(j,k,p)*wghtslg_pwo(ii) 

       end do

      V_retarGV0_jpp=regsum

  end function V_retarGV0_jpp



  function V_retarGV0_jmm(j,p_,p,E)
      real(NER) :: V_retarGV0_jmm
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     
     if(j<0) then
      V_retarGV0_jmm=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*V_retar_jpm(j,k,p_)*OPE_jpm(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*V_retar_jmm(j,p_,k)*OPE_jmm(j,k,p)*wghtslg_pwo(ii) 

       end do

      V_retarGV0_jmm=regsum

  end function V_retarGV0_jmm



   function V_retarGV0_jpm(j,p_,p,E)
      real(NER) :: V_retarGV0_jpm
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     

     if(j<0) then
      V_retarGV0_jpm=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*V_retar_jpm(j,p_,k)*OPE_jpp(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*V_retar_jmm(j,p_,k)*OPE_jpm(j,k,p)*wghtslg_pwo(ii) 

       end do

      V_retarGV0_jpm=regsum

  end function V_retarGV0_jpm



   function V_retarGV0_jmp(j,p_,p,E)
      real(NER) :: V_retarGV0_jmp
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     

     if(j<0) then
      V_retarGV0_jmp=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*V_retar_jpp(j,k,p_)*OPE_jpm(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*V_retar_jmp(j,p_,k)*OPE_jmm(j,k,p)*wghtslg_pwo(ii) 

       end do

      V_retarGV0_jmp=regsum

  end function V_retarGV0_jmp
end module V_retarGV0_pwd

