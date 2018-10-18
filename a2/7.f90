! in the module the potential is V_retardation


module V_retar_pwd

  
  use ope_pwd
  implicit none


  contains
 function WT1abc(q0sqr,qsqr)
    real(NER) ::  WT1abc
    real(NER),intent(in) :: q0sqr,qsqr
    real(NER) :: factor,PC_mpi_sqr
 
    factor=(PC_gA/PC_fpi)*(PC_gA/PC_fpi)*0.25_NER
    PC_mpi_sqr = PC_mpi * PC_mpi

    WT1abc=factor*q0sqr/((qsqr +PC_mpi_sqr)*(qsqr +PC_mpi_sqr))
 end function WT1abc 


 function WY1abc(q0,qsqr)
    real(NER) ::  WY1abc,z
    real(NER),intent(in) :: q0,qsqr
    real(NER) :: factor,PC_mpi_sqr,m_N
    factor=(PC_gA/PC_fpi)*(PC_gA/PC_fpi)*0.25_NER
    PC_mpi_sqr = PC_mpi * PC_mpi
    m_N=PC_mN

    WY1abc=factor*q0/(2.0_NER*m_N)/(qsqr+PC_mpi_sqr )
 end function WY1abc



 function V_retar_j0j(j,k,p)
    real(NER) :: V_retar_j0j
    real(NER),intent(in) :: k,p
    integer,intent(in) :: j
    real(NER) :: relation_factor
    integer :: ii
    real(NER) :: regsum,z,q0sqr,qsqr
    relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
    q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))
     

    if(j<0) then
      V_retar_j0j=0.0_NER
      return
    end if

      call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
         z=mshlg_pwo(ii)
         qsqr=k*k+p*p-2.0_NER*k*p*z
          regsum=regsum+2.0_NER*PI_NE*WT1abc(q0sqr,qsqr)*(-(k*k+p*p)+2.0_NER*k*p*z)*lgndr_table_pwo(j, ii)*wghtslg_pwo(ii)
       end do

   if(mod(j,2)==0) then
      
      regsum=-3_NER*regsum
    
    end if

    V_retar_j0j=relation_factor*regsum  
  end function V_retar_j0j



 function V_retar_j1j(j,k,p)
    real(NER) :: V_retar_j1j
    real(NER),intent(in) :: k,p
    integer,intent(in) :: j
    real(NER) :: relation_factor
    integer :: ii
    real(NER) :: regsum,z,q0sqr,qsqr
    relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
    q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))

    if(j<0) then
      V_retar_j1j=0.0_NER
      return
    end if

     call initmshlg_pwo()
       regsum=0.0_NER

     do ii=1,Nlg_pwo,1
          z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
        regsum=regsum+2.0_NER*PI_NE*WT1abc(q0sqr,qsqr)*((k*k+p*p)*lgndr_table_pwo(j, ii)&
        &-2.0_NER/(2.0_NER*j+1.0_NER)*j*k*p*lgndr_table_pwo(j+1, ii)&
       &-2.0_NER/(2.0_NER*j+1.0_NER)*(j+1.0_NER)*k*p*lgndr_table_pwo(j-1, ii))*wghtslg_pwo(ii)
      end do

    if(mod(j,2)==1) then
      
      regsum=-3_NER*regsum
    end if

     V_retar_j1j=relation_factor*regsum
 end function V_retar_j1j



 function V_retar_jpp(j,k,p)
   real(NER) :: V_retar_jpp
     real(NER),intent(in) :: k,p
      integer,intent(in) :: j
   real(NER) :: relation_factor
        integer :: ii
     real(NER) :: regsum,z,q0sqr,qsqr
      relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
   q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))

        if(j<0) then
        V_retar_jpp=0.0_NER
       return
       end if

    call initmshlg_pwo()
       regsum=0.0_NER
         do ii=1,Nlg_pwo,1
       z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
      regsum=regsum+2.0_NER*PI_NE*WT1abc(q0sqr,qsqr)*&
       &1.0_NER/(2.0_NER*j+1.0_NER)*(-(k*k+p*p)*lgndr_table_pwo(j+1, ii)&
     &+2.0_NER*k*p*lgndr_table_pwo(j, ii) )*wghtslg_pwo(ii)
      end do

    if(mod(j,2)==1) then
         
      regsum=-3_NER*regsum
    
    end if

    V_retar_jpp=relation_factor*regsum
 end function V_retar_jpp



 function V_retar_jmm(j,k,p)
    real(NER) :: V_retar_jmm
    real(NER),intent(in) :: k,p
       integer,intent(in) :: j
     real(NER) :: relation_factor
   integer :: ii
       real(NER) :: regsum,z,q0sqr,qsqr
     relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
     q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))

    if(j<0) then
       V_retar_jmm=0.0_NER
      return
   end if

       call initmshlg_pwo()
      regsum=0.0_NER
       do ii=1,Nlg_pwo,1
      z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
         regsum=regsum+2.0_NER*PI_NE*WT1abc(q0sqr,qsqr)*&
      &1.0_NER/(2.0_NER*j+1.0_NER)*((k*k+p*p)*lgndr_table_pwo(j-1, ii)&
     &-2.0_NER*k*p*lgndr_table_pwo(j, ii) )*wghtslg_pwo(ii)
      end do

    if(mod(j,2)==1) then
        
        regsum=-3_NER*regsum
   
    end if

      V_retar_jmm=relation_factor*regsum
 end function V_retar_jmm



 function V_retar_jpm(j,k,p)
    real(NER) :: V_retar_jpm
   real(NER),intent(in) :: k,p
       integer,intent(in) :: j
        real(NER) :: relation_factor
       integer :: ii
       real(NER) :: regsum,z,q0sqr,qsqr,q0
       relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
       q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))
        q0=(p*p-k*k)/(2.0_NER*PC_mN)

   if(j<0) then
       V_retar_jpm=0.0_NER
         return
   
    end if

       call initmshlg_pwo()
       regsum=0.0_NER
        do ii=1,Nlg_pwo,1
       z=mshlg_pwo(ii)
       qsqr=k*k+p*p-2.0_NER*k*p*z
         regsum=regsum+(-2.0_NER*PI_NE)*&
       &2.0_NER/(2.0_NER*j+1.0_NER)*sqrt(j*(j+1.0_NER))*WT1abc(q0sqr,qsqr)*((p*p*lgndr_table_pwo(j+1, ii)&
       &+k*k*lgndr_table_pwo(j-1, ii)-2.0_NER*k*p*lgndr_table_pwo(j, ii) )&
       &- WY1abc(q0,qsqr)*2.0_NER*k*p*(lgndr_table_pwo(j+1, ii)-lgndr_table_pwo(j-1, ii)))*wghtslg_pwo(ii)
       end do

    
    if(mod(j,2)==1) then
        
        regsum=-3_NER*regsum
    
    end if

       V_retar_jpm=relation_factor*regsum
  end function V_retar_jpm



  function V_retar_jmp(j,k,p)
   real(NER) :: V_retar_jmp
       real(NER),intent(in) :: k,p
       integer,intent(in) :: j
     real(NER) :: relation_factor
      integer :: ii
         real(NER) :: regsum,z,q0sqr,qsqr,q0
       relation_factor=PC_mN/(8.0_NER*PI_NE*PI_NE*PI_NE)
       q0sqr=((p*p-k*k)/(2.0_NER*PC_mN))*((p*p-k*k)/(2.0_NER*PC_mN))
    q0=(p*p-k*k)/(2.0_NER*PC_mN)

       if(j<0) then
       V_retar_jmp=0.0_NER
         return
        end if

      call initmshlg_pwo()
    regsum=0.0_NER
      do ii=1,Nlg_pwo,1
        z=mshlg_pwo(ii)
         qsqr=k*k+p*p-2.0_NER*k*p*z
        regsum=regsum+(-2.0_NER*PI_NE)*&
      &2.0_NER/(2.0_NER*j+1.0_NER)*sqrt(j*(j+1.0_NER))*(WT1abc(q0sqr,qsqr)*(k*k*lgndr_table_pwo(j+1, ii)&
     &+p*p*lgndr_table_pwo(j-1, ii)-2.0_NER*k*p*lgndr_table_pwo(j, ii) )&
      &+ WY1abc(q0,qsqr)*2.0_NER*k*p*(lgndr_table_pwo(j+1, ii)-lgndr_table_pwo(j-1, ii)))*wghtslg_pwo(ii)
      end do

   if(mod(j,2)==1) then
       
       regsum=-3_NER*regsum
    

    end if
      
       V_retar_jmp=relation_factor*regsum
  end function V_retar_jmp



end module V_retar_pwd