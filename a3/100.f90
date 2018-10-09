program D1
use V_retar_pwd
implicit none
real(NER) :: k,p
integer :: j
read(*,*) j,k,p
write(*,*) V_retar_j0j(j,k,p)
end program D1
