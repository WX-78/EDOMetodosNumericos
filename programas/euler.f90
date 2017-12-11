module constantes
	implicit none
	real ::  omg = 1.4, psi0=1.0, z0=0.0
end module constantes

module derivadas ! Um modulo para definir as funcoes derivadas
! Esse modulo deve ser usado na subrotina que ira resolver seu 
! sistema de EDOs
	implicit none
	contains
	function f(s,z) ! Aqui sao definidas as funcoes derivadas
		use constantes
		real :: s ! variavel independente
		real :: z(2) ! vetor com variaveis dependentes
		real :: f(2) ! cada f corresponde ao lado direito da EDO 
		             ! dz/dx = f
		f(1) = z(2)  ! v
		f(2) = -z(1)*omg**2 
	end function f
end module derivadas

module euler 
	implicit none
	contains
	subroutine metodo(h,t,y0,y) 
		use derivadas
		implicit none
		real :: y(:,:)      ! variaveis dependentes e derivadas
		real :: h           ! tamanho do passo
		real :: t(:)        ! vetor de valores da variavel 
		                    ! independente
		real :: y0(:)       ! condicoes iniciais
		integer :: Npassos, Neqs !No. de passos e No. de equacoes
		integer :: i
		real, allocatable, dimension(:) :: k1,k2
		Npassos = size(t)
		Neqs = size(y0)
		allocate ( k1(Neqs), k2(Neqs))
		y(:,1)=y0 ! valores iniciais
		do i=1, Npassos-1
			y(:,i+1) = y(:,i) + h*f(t(i),y(:,i))
		end do
	end subroutine metodo
end module euler 

program passing
	use constantes
	use euler 
	implicit none
	real,parameter :: h=0.12 !alterar aqui o tamanho do passo
	integer, parameter :: Npts=2000+1
	real :: yinit(2)
	real :: y(2,Npts)
	integer :: i
	real :: t(Npts)
	t(1) = 0.0
	do i = 2,Npts
		t(i) = t(i-1)+h
	end do
	yinit(1) = psi0
	yinit(2) = z0
	call metodo(h,t,yinit,y)
	write(*,*) "#   t              x           dx/dt"
	do i=1,Npts
		write(*,*) t(i), y(1,i), y(2,i)
	end do
end program passing

