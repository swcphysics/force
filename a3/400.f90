program D1
use V0GV_retar_pwd
use V_retarGV0_pwd
implicit none
real(NER) :: p_,p,E
integer :: j
read(*,*) j,p_,p,E
write(*,*) V0GV_retar_j0j(j,p_,p,E),V_retarGV0_j0j(j,p,p_,E)
end program D1
