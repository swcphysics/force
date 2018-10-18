module V0GV_retar_pwd
	use V_retar_pwd
	implicit none

 contains


  function V0GV_retar_j0j(j,p_,p,E)
      real(NER) :: V0GV_retar_j0j
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
 

     if(j<0) then
      V0GV_retar_j0j=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*OPE_j0j(j,p_,k)*V_retar_j0j(j,k,p)*wghtslg_pwo(ii)

       end do

       V0GV_retar_j0j=regsum
    
  end function V0GV_retar_j0j



  function V0GV_retar_j1j(j,p_,p,E)
      real(NER) :: V0GV_retar_j1j
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
  

     if(j<0) then
      V0GV_retar_j1j=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*OPE_j1j(j,p_,k)*V_retar_j1j(j,k,p)*wghtslg_pwo(ii)

       end do

       V0GV_retar_j1j=regsum

  end function V0GV_retar_j1j



    function V0GV_retar_jpp(j,p_,p,E)
      real(NER) :: V0GV_retar_jpp
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     

     if(j<0) then
      V0GV_retar_jpp=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*OPE_jpp(j,p_,k)*V_retar_jpp(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*OPE_jpm(j,p_,k)*V_retar_jmp(j,k,p)*wghtslg_pwo(ii) 

       end do

      V0GV_retar_jpp=regsum

  end function V0GV_retar_jpp



  function V0GV_retar_jmm(j,p_,p,E)
      real(NER) :: V0GV_retar_jmm
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     
     if(j<0) then
      V0GV_retar_jmm=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*OPE_jpm(j,k,p_)*V_retar_jpm(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*OPE_jmm(j,p_,k)*V_retar_jmm(j,k,p)*wghtslg_pwo(ii) 

       end do

      V0GV_retar_jmm=regsum

  end function V0GV_retar_jmm



   function V0GV_retar_jpm(j,p_,p,E)
      real(NER) :: V0GV_retar_jpm
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     

     if(j<0) then
      V0GV_retar_jpm=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*OPE_jpp(j,p_,k)*V_retar_jpm(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*OPE_jpm(j,p_,k)*V_retar_jmm(j,k,p)*wghtslg_pwo(ii) 

       end do

      V0GV_retar_jpm=regsum

  end function V0GV_retar_jpm



   function V0GV_retar_jmp(j,p_,p,E)
      real(NER) :: V0GV_retar_jmp
     real(NER),intent(in) :: p_,p,E
     integer,intent(in) :: j
     integer :: ii
     real(NER) :: regsum,k
     

     if(j<0) then
      V0GV_retar_jmp=0.0_NER
      return
     end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         k=mshlg_pwo(ii)
        
          regsum=regsum+k*k/(E-k*k/PC_mN)*OPE_jpm(j,k,p_)*V_retar_jpp(j,k,p)*wghtslg_pwo(ii)   &
         &   +k*k/(E-k*k/PC_mN)*OPE_jmm(j,p_,k)*V_retar_jmp(j,k,p)*wghtslg_pwo(ii) 

       end do

      V0GV_retar_jmp=regsum

  end function V0GV_retar_jmp
end module V0GV_retar_pwd

