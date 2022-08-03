! ---------------------------------------------------------------------
! Copyright (C) 2014-2022 Universidad de Las Palmas de Gran Canaria:
!                         Jacob D.R. Bordon
!                         Guillermo M. Alamo
!                         Juan J. Aznarez
!                         Orlando Maeso.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ---------------------------------------------------------------------

!! @author Juan Jose Aznarez (juanjose.aznarez@ulpgc.es), Jacob David Rodriguez Bordon (jacobdavid.rodriguezbordon@ulpgc.es)
!!
!! @version 2.0
!!
!! <b> This module implements some incident fields for elastic solids. </b>
module fbem_harela_incident_field

  ! Fortran 2003 standard
  use iso_fortran_env

  ! fbem modules
  use fbem_numerical
  use fbem_string_handling

  ! No implicit variables are allowed in the module
  implicit none

  ! By default all are private
  private

  ! Subroutines
  public :: fbem_harela_incident_plane_wave_multilayered_amplitudes
  public :: fbem_harela_incident_plane_wave_vertical_multilayered
  public :: fbem_harela_incident_plane_wave
  public :: fbem_harela_c_rayleigh
  public :: harela_calculate_rayleigh_wavenumber
  public :: fbem_ela_properties

  ! Parameters
  integer, parameter, public :: fbem_harela_sh_wave      =1 ! indice asociado al tipo de onda SH
  integer, parameter, public :: fbem_harela_p_wave       =2 ! indice asociado al tipo de onda P
  integer, parameter, public :: fbem_harela_sv_wave      =3 ! indice asociado al tipo de onda SV
  integer, parameter, public :: fbem_harela_rayleigh_wave=4 ! indice asociado al tipo de onda Rayleigh

contains

  !! Calculate waves amplitudes
  subroutine fbem_harela_incident_plane_wave_multilayered_amplitudes(omega,np,nl,xl,lambdal,mul,rhol,type_of_wave,A)
    implicit none
    ! I/O variables
    real(kind=real64)    :: omega
    integer              :: np ! +-1,+-2,+-3
    integer              :: nl
    real(kind=real64)    :: xl(nl)
    complex(kind=real64) :: lambdal(nl)
    complex(kind=real64) :: mul(nl)
    real(kind=real64)    :: rhol(nl)
    integer              :: type_of_wave
    complex(kind=real64) :: A(nl,2)
    ! Local variables
    integer               :: j
    real(kind=real64)     :: h(nl-1)
    complex(kind=real64)  :: cimag, e(nl), k(nl), R
    !
    ! Initialization
    !
    cimag=(0.d0,1.d0)
    ! Check np
    if ((np.gt.3).or.(np.eq.0).or.(np.lt.-3)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'invalid value of np')
    end if
    ! Calculate layer thickness
    do j=1,nl-1
      h(j)=xl(j+1)-xl(j)
    end do
    ! Check layer thickness (indicate invalid layer coordinates)
    if (np.lt.0) then
      if (minval(h).le.0.d0) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'calculated layer thickness less or equal to zero.')
      end if
    else
      if (maxval(h).ge.0.d0) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,'calculated layer thickness less or equal to zero.')
      end if
    end if
    ! Assign wavenumber and elastic constant
    select case (type_of_wave)
      case (fbem_harela_sh_wave,fbem_harela_sv_wave)
        do j=1,nl
          e(j)=mul(j)
          k(j)=omega/sqrt(mul(j)/rhol(j))
        end do
      case (fbem_harela_p_wave)
        do j=1,nl
          e(j)=lambdal(j)
          k(j)=omega/sqrt((lambdal(j)+2.d0*mul(j))/rhol(j))
        end do
    end select
    ! Calculate amplitudes
    A(1,1)=0.5d0
    A(1,2)=0.5d0
    do j=1,nl-1
      R=e(j)*k(j)/(e(j+1)*k(j+1))
      A(j+1,1)=0.5d0*((1.d0+R)*A(j,1)*exp(cimag*k(j)*h(j))+(1.d0-R)*A(j,2)*exp(-cimag*k(j)*h(j)))
      A(j+1,2)=0.5d0*((1.d0-R)*A(j,1)*exp(cimag*k(j)*h(j))+(1.d0+R)*A(j,2)*exp(-cimag*k(j)*h(j)))
    end do
  end subroutine fbem_harela_incident_plane_wave_multilayered_amplitudes

  subroutine fbem_harela_incident_plane_wave_vertical_multilayered(d,omega,np,nl,xl,lambdal,mul,rhol,alfa,type_of_wave,A,x,n,ui,ti)
    implicit none
    ! I/O variables
    integer              :: d
    real(kind=real64)    :: omega
    integer              :: np ! +-1,+-2,+-3
    integer              :: nl
    real(kind=real64)    :: xl(nl)
    complex(kind=real64) :: lambdal(nl)
    complex(kind=real64) :: mul(nl)
    real(kind=real64)    :: rhol(nl)
    real(kind=real64)    :: alfa
    integer              :: type_of_wave
    complex(kind=real64) :: A(nl,2)
    real(kind=real64)    :: x(3)
    real(kind=real64)    :: n(3)
    complex(kind=real64) :: ui(3)
    complex(kind=real64) :: ti(3)
    ! Local variables
    integer               :: j, kl
    complex(kind=real64)  :: cimag
    complex(kind=real64)  :: lambda, mu, ks, kp, Ap, Am, expp, expm
    real(kind=real64)     :: tol,rho, h, dx
    complex(kind=real64)  :: dui(3,3), eps(6), sig(6)
    !
    ! Initialization
    !
    ! 2D problem:
    !   x2 (y) in 2d == x3 (z) in 3d
    !   x1 (x) in 2d == x2 (y) in 3d
    !
    if (d.eq.2) then
      x(3)=x(2)
      n(3)=n(2)
      x(2)=x(1)
      n(2)=n(1)
      x(1)=0.d0
      n(1)=0.d0
    end if
    tol=1.d-6
    cimag=(0.d0,1.d0)
    ! Check np
    if ((np.ne.3).and.(np.ne.-3)) then
      call fbem_error_message(error_unit,0,__FILE__,__LINE__,'only np=3 or np=-3 is available')
    end if
    ! Check layer thickness (indicate invalid layer coordinates)
    if (np.gt.0) then
      do j=1,nl-1
        if ((xl(j)-xl(j+1)).le.0.d0) then
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'layer thickness less or equal to zero.')
        end if
      end do
    else
      do j=1,nl-1
        if ((xl(j+1)-xl(j)).le.0.d0) then
          call fbem_error_message(error_unit,0,__FILE__,__LINE__,'layer thickness less or equal to zero.')
        end if
      end do
    end if
    !
    ! Find in which layer is the evaluation point
    !
    if (np.gt.0) then
      if (x(3).gt.(xl(1)+tol)) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
        'trying to evaluate incident field at a point not belonging to the layered half-space.')
      else
        do kl=nl,1,-1
          if (x(3).le.(xl(kl)+tol)) then
            lambda=lambdal(kl)
            mu=mul(kl)
            rho=rhol(kl)
            ks=omega/sqrt(mu/rho)
            kp=omega/sqrt((lambda+2.d0*mu)/rho)
            Ap=A(kl,1)
            Am=A(kl,2)
            exit
          end if
        end do
      end if
    else
      if (x(3).lt.(xl(1)-tol)) then
        call fbem_error_message(error_unit,0,__FILE__,__LINE__,&
        'trying to evaluate incident field at a point not belonging to the layered half-space.')
      else
        do kl=nl,1,-1
          if (x(3).ge.(xl(kl)-tol)) then
            lambda=lambdal(kl)
            mu=mul(kl)
            rho=rhol(kl)
            ks=omega/sqrt(mu/rho)
            kp=omega/sqrt((lambda+2.d0*mu)/rho)
            Ap=A(kl,1)
            Am=A(kl,2)
            exit
          end if
        end do
      end if
    end if
    !
    ! CALCULATE DISPLACEMENTS AND DISPLACEMENT DERIVATIVE
    !
    ui=0.d0
    dui=0.d0
    dx=x(3)-xl(kl)
    select case (type_of_wave)
      !
      ! SH WAVE (symmetric with respect to plane XZ)
      !
      case (fbem_harela_sh_wave)
        if (nint(alfa).ge.0) then
          expp=exp( cimag*ks*dx)
          expm=exp(-cimag*ks*dx)
          ui(1)=Ap*expp+Am*expm
          dui(1,3)=cimag*ks*(Ap*expp-Am*expm)
        end if
      !
      ! SV WAVE (antisymmetric with respect to plane XZ)
      !
      case (fbem_harela_sv_wave)
        if (nint(alfa).le.0) then
          expp=exp( cimag*ks*dx)
          expm=exp(-cimag*ks*dx)
          ui(2)=Ap*expp+Am*expm
          dui(2,3)=cimag*ks*(Ap*expp-Am*expm)
        end if
      !
      ! P WAVE (symmetric with respect to plane XZ)
      !
      case (fbem_harela_p_wave)
        if (nint(alfa).ge.0) then
          expp=exp( cimag*kp*dx)
          expm=exp(-cimag*kp*dx)
          ui(3)=Ap*expp+Am*expm
          dui(3,3)=cimag*kp*(Ap*expp-Am*expm)
        end if
    end select
    !
    ! Tensor de deformaciones
    !
    eps(1)=dui(1,1)                  ! epsilon 11
    eps(2)=0.5d0*(dui(1,2)+dui(2,1)) ! epsilon 12
    eps(3)=0.5d0*(dui(1,3)+dui(3,1)) ! epsilon 13
    eps(4)=dui(2,2)                  ! epsilon 22
    eps(5)=0.5d0*(dui(2,3)+dui(3,2)) ! epsilon 23
    eps(6)=dui(3,3)                  ! epsilon 33
    !
    ! Tensor de tensiones
    !
    sig(1)=2.d0*mu*eps(1)+lambda*(eps(1)+eps(4)+eps(6)) ! sigma 11
    sig(2)=2.d0*mu*eps(2)                               ! sigma 12
    sig(3)=2.d0*mu*eps(3)                               ! sigma 13
    sig(4)=2.d0*mu*eps(4)+lambda*(eps(1)+eps(4)+eps(6)) ! sigma 22
    sig(5)=2.d0*mu*eps(5)                               ! sigma 23
    sig(6)=2.d0*mu*eps(6)+lambda*(eps(1)+eps(4)+eps(6)) ! sigma 33
    !
    ! Vector tension
    !
    ti(1)=sig(1)*n(1)+sig(2)*n(2)+sig(3)*n(3)
    ti(2)=sig(2)*n(1)+sig(4)*n(2)+sig(5)*n(3)
    ti(3)=sig(3)*n(1)+sig(5)*n(2)+sig(6)*n(3)
    !
    !  2D problem:
    !    x2 [y] en 2d == x3 [z] en 3d
    !    x1 [x] en 2d == x2 [y] en 3d
    !
    if (d.eq.2) then
      ui(1)=ui(2)
      ti(1)=ti(2)
      ui(2)=ui(3)
      ti(2)=ti(3)
    end if
  end subroutine fbem_harela_incident_plane_wave_vertical_multilayered

  function fbem_harela_c_rayleigh(cs,nu)
    implicit none
    complex*16 cs
    real*8     nu
    complex*16 fbem_harela_c_rayleigh
    real*8     rgammas
    call obtenerc(nu,rgammas)
    fbem_harela_c_rayleigh=rgammas*cs
  end function fbem_harela_c_rayleigh

  ! --------------------------------------------------------------------------------------------------------------------------------
  !
  !     Calcula campo de movimientos y tensores de tension de un campo
  !     incidente de ondas planas en un semiespacio viscoelastico,
  !     incluyendo posible simetria respecto al plano XZ.
  !
  !     Variables de entrada:
  !       - d           : dimension del problem (2D o 3D)
  !       - k1,k2       : numeros de onda de la onda P y S, respectivamente
  !       - mu,lambda,nu: propiedades del medio viscoelastico
  !       - alfa        : control de la simetria:
  !                         alfa=-1. => solo la parte antisimetrica con respecto al plano XZ (eje x en 2D)
  !                         alfa= 0. => sin simetria (campo incidente completo)
  !                         alfa= 1. => solo la parte simetrica con respecto al plano XZ (eje x en 2D)
  !       - type_of_wave: tipo de onda:
  !                         type_of_wave=1 => Onda SH
  !                         type_of_wave=2 => Onda P
  !                         type_of_wave=3 => Onda SV
  !                         type_of_wave=4 => Onda de Rayleigh
  !       - fi          : sea la proyeccion del vector director de la onda plana
  !                       incidente sobre el plano XY, FI es el angulo entre esta y
  !                       el eje y (es el angulo pi/2-angulo azimutal de las coord.
  !                       esfericas). En problemas 2D debe ser nulo.
  !       - theta       : angulo del vector director de la onda plana con
  !                       respecto al plano XY (es el angulo pi/2-angulo polar
  !                       de las coord. esfericas)
  !       - x(3)        : coordenadas del punto sobre el que calcular el campo
  !                       incidente
  !       - z_fs        : coordenada de posicion de la superficie libre
  !       - n(3)        : vector normal del punto
  !
  !     Variables de salida:
  !       - ui(3) : vector desplazamiento del campo incidente
  !       - ti(3) : vector tension del campo incidente
  !
  !     Autor:
  !       Juan Jose Aznarez Gonzalez
  !       Instituto Universitario SIANI
  !       Universidad de Las Palmas de Gran Canaria
  !
  ! --------------------------------------------------------------------------------------------------------------------------------

  subroutine fbem_harela_incident_plane_wave(d,k1,k2,mu,lambda,nu,space,alfa,type_of_wave,fi,theta,x,z_fs,n,ui,ti)
    implicit none
    ! variables de entrada y salida
    integer    d
    complex*16 k1
    complex*16 k2
    complex*16 mu
    complex*16 lambda
    real*8     nu
    integer    space
    real*8     alfa
    integer    type_of_wave
    real*8     theta
    real*8     fi
    real*8     x(3)
    real*8     z_fs
    real*8     n(3)
    complex*16 ui(3)
    complex*16 ti(3)
    ! variables locales
    complex*16 c0,s0,c1,s1,c2,s2                ! cosenos y senos de angulos
    complex*16 s20,s22,c20,c22,s40              ! cosenos y senos en complejo
    complex*16 dd(3,3),pp(3,3),dg(3,3),pg(3,3)  ! vectores directores de la onda
    complex*16 aa(3)                            ! amplitudes de las ondas implicadas
    complex*16 cimag                            ! numero imaginario
    complex*16 zcs,zcp,zz(3)                    ! numero de onda de las ondas implicadas
    complex*16 exz(3),eys(3),eya(3)
    integer    non                              ! numero de ondas implicadas
    complex*16 du(3,3)                          ! derivadas de desplazamientos (u_{i,j})
    complex*16 eps(6)                           ! tensor de deformaciones (eps_{ij})
    complex*16 sigma(6)                         ! tensor de tensiones (sigma_{ij})
    complex*16 zcr                              ! terminos relacionados con la onda de rayleigh
    real*8     xr,yr,zr                         ! coordenadas relativas
    real*8     rk
    real*8     rgammas,gammas,gammap
    integer    i,j                              ! contadores
    !
    ! Inicializacion
    !
    ! Si el problema es 2d:
    !   x2 (y) en 2d == x3 (z) en 3d
    !   x1 (x) en 2d == x2 (y) en 3d
    if (d.eq.2) then
      x(3)=x(2)
      n(3)=n(2)
      x(2)=x(1)
      n(2)=n(1)
      x(1)=0.d0
      n(1)=0.d0
    end if
    ! Unidad imaginaria
    cimag=(0.d0,1.d0)
    ! Inicializacion variables
    ui=(0.d0,0.d0)
    ti=(0.d0,0.d0)
    sigma=(0.d0,0.d0)
    ! Coordenadas relativas
    xr=x(1)
    yr=x(2)
    zr=x(3)-z_fs
    ! Numeros de onda multiplicados por i
    zcp=cimag*k1
    zcs=cimag*k2
    ! Constante
    rk=dsqrt((1.0d0-2.0d0*nu)/(2.0d0*(1.0d0-nu)))
    ! Inicializacion de los vectores directores y derivadas del desplazamiento
    pp=(0.d0,0.d0)
    dd=(0.d0,0.d0)
    pg=(0.d0,0.d0)
    dg=(0.d0,0.d0)
    du=(0.d0,0.d0)
    !
    ! -------
    ! Onda SH
    ! -------
    !
    if (type_of_wave.eq.fbem_harela_sh_wave) then
      non=2
      c0=dcos(theta)
      s0=dsin(theta)
      c1=dcos(theta)
      s1=dsin(theta)
      !
      ! vectores directores p y d
      !
      ! Onda incidente SH
      pp(2,1)=c0
      pp(3,1)=s0
      dd(1,1)=(1.d0,0.d0)
      ! Onda reflejada SH
      pp(2,2)=c1
      pp(3,2)=-s1
      dd(1,2)=(1.d0,0.d0)
      ! Amplitud onda incidente (SH)
      aa(1)=(1.d0,0.d0)
      ! Amplitud onda reflejada (SH)
      aa(2)=(1.d0,0.d0)
      ! no tiene sentido
      aa(3)=(0.d0,0.d0)
      ! numeros de onda (i*k) asociados
      zz(1)=zcs
      zz(2)=zcs
      zz(3)=(0.d0,0.d0)
    end if
    !
    ! ------
    ! Onda P
    ! ------
    !
    if (type_of_wave.eq.fbem_harela_p_wave) then
      non=3
      c0=dcos(theta)
      s0=dsin(theta)
      c1=dcos(theta)
      s1=dsin(theta)
      c2=rk*c0
      s2=zsqrt(1.d0-c2*c2)
      s20=2.d0*s0*c0
      s22=2.d0*s2*c2
      c22=c2*c2-s2*s2
      !
      ! Vectores directores p y d
      !
      ! onda incidente P
      pp(2,1)=c0
      pp(3,1)=s0
      dd(2,1)=c0
      dd(3,1)=s0
      ! onda reflejada P
      pp(2,2)=c1
      pp(3,2)=-s1
      dd(2,2)=c1
      dd(3,2)=-s1
      ! onda reflejada SV
      pp(2,3)=c2
      pp(3,3)=-s2
      dd(2,3)=-s2
      dd(3,3)=-c2
      ! amplitud onda incidente (P)
      aa(1)=(1.d0,0.d0)
      ! amplitud onda reflejada (P)
      aa(2)=(rk*rk*s20*s22-c22*c22)/(rk*rk*s20*s22+c22*c22)
      ! amplitud onda reflejada (SV)
      aa(3)=(2.d0*rk*s20*c22)/(rk*rk*s20*s22+c22*c22)
      ! numeros de onda (i*k) asociados
      zz(1)=zcp
      zz(2)=zcp
      zz(3)=zcs
    end if
    !
    ! -------
    ! Onda SV
    ! -------
    !
    if (type_of_wave.eq.fbem_harela_sv_wave) then
      non=3
      !
      ! Angulo superior al critico
      !

      ! ojo, aqui se trataba la comparacion de angulos mal creo....., hay que comparar
      ! con el original y ver....., creo que es solo valido para
      ! angulos theta 0<=theta<=90, al poner el dabs en la comparacion, ahora se admite
      ! 0<=theta<=180. Aun asi sigue habiendo algo mal cuando 180<=theta<=360, las soluciones salen
      ! salen revertidos. Es una tonteria pero hay que tenerlo en cuenta....


      if (dabs(dcos(theta)).le.rk) then
        c0=dcos(theta)
        s0=dsin(theta)
        c1=dcos(theta)
        s1=dsin(theta)
        c2=(1.d0/rk)*c0
        s2=zsqrt(1.d0-c2*c2)
      end if
      !
      ! Angulo inferior al critico
      !
      if (dabs(dcos(theta)).gt.rk) then
        c0=dcos(theta)
        s0=dsin(theta)
        c1=dcos(theta)
        s1=dsin(theta)
        c2=(1.d0/rk)*c0
        s2=-cimag*dsqrt(dreal(c2*c2)-1.d0)
      end if
      !
      s20=2.d0*s0*c0
      s22=2.d0*s2*c2
      c20=c0*c0-s0*s0
      s40=2.d0*s20*c20
      !
      ! Vectores directores p y d
      !
      ! Onda incidente SV
      pp(2,1)=c0
      pp(3,1)=s0
      dd(2,1)=s0
      dd(3,1)=-c0
      ! Onda reflejada SV
      pp(2,2)=c1
      pp(3,2)=-s1
      dd(2,2)=-s1
      dd(3,2)=-c1
      ! Onda reflejada P
      pp(2,3)=c2
      pp(3,3)=-s2
      dd(2,3)=c2
      dd(3,3)=-s2
      ! Amplitud onda incidente (SV)
      aa(1)=(1.d0,0.d0)
      ! Amplitud onda reflejada (SV)
      aa(2)=(rk*rk*s20*s22-c20*c20)/(rk*rk*s20*s22+c20*c20)
      ! Amplitud onda reflejada (P)
      aa(3)=-(rk*s40)/(rk*rk*s20*s22+c20*c20)
      ! Numeros de onda (i*k) asociados
      zz(1)=zcs
      zz(2)=zcs
      zz(3)=zcp
    end if
    !
    ! ----------------
    ! Onda de Rayleigh
    ! ----------------
    !
    if (type_of_wave.eq.fbem_harela_rayleigh_wave) then
      non=2
      !
      c0=dcos(theta)
      c1=dcos(theta)
      !
      ! Velocidad de las ondas de Rayleigh
      !
      call obtenerc(nu,rgammas)
      gammas=rgammas*rgammas
      gammap=rk*rk*gammas
      zcr=(1.d0/rgammas)*zcs
      !
      ! Vectores directores p y d
      !
      ! Onda incidente Rayleigh
      pp(2,1)=c0
      pp(3,1)=cimag*dsqrt(1.d0-gammas)
      dd(2,1)=c0
      dd(3,1)=cimag/dsqrt(1.d0-gammas)
      ! Onda reflejada Rayleigh
      pp(2,2)=c1
      pp(3,2)=cimag*dsqrt(1.d0-gammap)
      dd(2,2)=c1
      dd(3,2)=cimag*dsqrt(1.d0-gammap)
      ! Amplitud onda incidente (Rayleigh)
      aa(1)=(1.d0,0.d0)
      ! Amplitud onda reflejada (Rayleigh)
      aa(2)=-2.d0/(2.d0-gammas)
      ! no tiene sentido
      aa(3)=(0.d0,0.d0)
      ! numeros de onda (i*k) asociados
      zz(1)=zcr
      zz(2)=zcr
      zz(3)=(0.d0,0.d0)
    end if
    !
    ! --------------------------
    ! Paso a coodenadas globales
    ! --------------------------
    !
    call pgydg(pp,dd,fi,pg,dg)
    !
    ! Onda incidente (1) y reflejada (2 y 3)
    !
    if (space.eq.1) non=1
    do i=1,non
      exz(i)=cdexp(-zz(i)*(pg(1,i)*xr+pg(3,i)*zr))
      eys(i)=0.5d0*(cdexp(-zz(i)*pg(2,i)*yr)+cdexp(zz(i)*pg(2,i)*yr))
      eya(i)=0.5d0*(cdexp(-zz(i)*pg(2,i)*yr)-cdexp(zz(i)*pg(2,i)*yr))
    end do
    !
    ! -----------------------------------
    ! Desplazamientos del campo incidente
    ! -----------------------------------
    !
    ! Componente simetrica
    !
    if (nint(alfa).ge.0) then
      do i=1,non
      ! Desplazamientos
      ui(1)=ui(1)+dg(1,i)*aa(i)*(exz(i)*eys(i))
      ui(2)=ui(2)+dg(2,i)*aa(i)*(exz(i)*eya(i))
      ui(3)=ui(3)+dg(3,i)*aa(i)*(exz(i)*eys(i))
      ! Matriz de derivadas de desplazamientos
      du(1,1)=du(1,1)+dg(1,i)*aa(i)*pg(1,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(1,2)=du(1,2)+dg(1,i)*aa(i)*pg(2,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(1,3)=du(1,3)+dg(1,i)*aa(i)*pg(3,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(2,1)=du(2,1)+dg(2,i)*aa(i)*pg(1,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(2,2)=du(2,2)+dg(2,i)*aa(i)*pg(2,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(2,3)=du(2,3)+dg(2,i)*aa(i)*pg(3,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(3,1)=du(3,1)+dg(3,i)*aa(i)*pg(1,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(3,2)=du(3,2)+dg(3,i)*aa(i)*pg(2,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(3,3)=du(3,3)+dg(3,i)*aa(i)*pg(3,i)*(-1.)*zz(i)*exz(i)*eys(i)
      end do
     end if
    !
    ! Componente antimetrica
    !
    if (nint(alfa).le.0) then
      do i=1,non
      ! Desplazamientos
      ui(1)=ui(1)+dg(1,i)*aa(i)*(exz(i)*eya(i))
      ui(2)=ui(2)+dg(2,i)*aa(i)*(exz(i)*eys(i))
      ui(3)=ui(3)+dg(3,i)*aa(i)*(exz(i)*eya(i))
      ! Matriz de derivadas de desplazamientos
      du(1,1)=du(1,1)+dg(1,i)*aa(i)*pg(1,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(1,2)=du(1,2)+dg(1,i)*aa(i)*pg(2,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(1,3)=du(1,3)+dg(1,i)*aa(i)*pg(3,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(2,1)=du(2,1)+dg(2,i)*aa(i)*pg(1,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(2,2)=du(2,2)+dg(2,i)*aa(i)*pg(2,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(2,3)=du(2,3)+dg(2,i)*aa(i)*pg(3,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(3,1)=du(3,1)+dg(3,i)*aa(i)*pg(1,i)*(-1.)*zz(i)*exz(i)*eya(i)
      du(3,2)=du(3,2)+dg(3,i)*aa(i)*pg(2,i)*(-1.)*zz(i)*exz(i)*eys(i)
      du(3,3)=du(3,3)+dg(3,i)*aa(i)*pg(3,i)*(-1.)*zz(i)*exz(i)*eya(i)
      end do
    end if
    !
    ! Tensor de deformaciones
    !
    eps(1)=du(1,1)                 ! epsilon 11
    eps(2)=0.5d0*(du(1,2)+du(2,1)) ! epsilon 12
    eps(3)=0.5d0*(du(1,3)+du(3,1)) ! epsilon 13
    eps(4)=du(2,2)                 ! epsilon 22
    eps(5)=0.5d0*(du(2,3)+du(3,2)) ! epsilon 23
    eps(6)=du(3,3)                 ! epsilon 33
    !
    ! Tensor de tensiones
    !
    sigma(1)=2.d0*mu*eps(1)+lambda*(eps(1)+eps(4)+eps(6)) ! sigma 11
    sigma(2)=2.d0*mu*eps(2)                               ! sigma 12
    sigma(3)=2.d0*mu*eps(3)                               ! sigma 13
    sigma(4)=2.d0*mu*eps(4)+lambda*(eps(1)+eps(4)+eps(6)) ! sigma 22
    sigma(5)=2.d0*mu*eps(5)                               ! sigma 23
    sigma(6)=2.d0*mu*eps(6)+lambda*(eps(1)+eps(4)+eps(6)) ! sigma 33
    !
    ! --------------
    ! Vector tension
    ! --------------
    !
    ti(1)=sigma(1)*n(1)+sigma(2)*n(2)+sigma(3)*n(3)
    ti(2)=sigma(2)*n(1)+sigma(4)*n(2)+sigma(5)*n(3)
    ti(3)=sigma(3)*n(1)+sigma(5)*n(2)+sigma(6)*n(3)
    !
    !  Si el problema es 2D:
    !    x2 [y] en 2d == x3 [z] en 3d
    !    x1 [x] en 2d == x2 [y] en 3d
    !
    if (d.eq.2) then
      ui(1)=ui(2)
      ti(1)=ti(2)
      ui(2)=ui(3)
      ti(2)=ti(3)
    end if
    return
  end subroutine fbem_harela_incident_plane_wave

  !
  ! Rutina auxiliar
  !
  subroutine pgydg(pp,dd,fi,pg,dg)
    implicit none
    complex*16  pp(3,3),dd(3,3)
    complex*16  pg(3,3),dg(3,3)
    real*8      rr(3,3)          ! matriz de giro
    real*8      fi
    integer     i,j,k
    rr(1,1)=dcos(fi)
    rr(1,2)=dsin(fi)
    rr(1,3)=0.d0
    rr(2,1)=-rr(1,2)
    rr(2,2)=rr(1,1)
    rr(2,3)=0.d0
    rr(3,1)=0.d0
    rr(3,2)=0.d0
    rr(3,3)=1.d0
    do k=1,3
      do i=1,3
        do j=1,3
          dg(i,k)=dg(i,k)+rr(i,j)*dd(j,k)
          pg(i,k)=pg(i,k)+rr(i,j)*pp(j,k)
        end do
      end do
    end do
    return
  end subroutine pgydg

  !
  ! Rutina auxiliar
  !
  subroutine obtenerc(nu,gammas)
    implicit none
    real*8 nu,gammas,gammas2
    complex*16 expr2,expr3,i
    complex*16 expr21,expr22,expr31,expr32
    complex*16 den,sub2,sb1
    i=(0.0d0,1.0d0)
    ! Denominador de todas las expresiones
    if (nu.lt.0.395d0) then
      sb1=-((-0.26308206488336366d0+nu)*(0.5939211404216125d0+(-0.23691793511663636d0+nu)*nu)/((-1.d0+nu)**3))
      sub2=(448.0d0-360.0d0/(1.0d0-nu)+235.1510153071851d0*cdsqrt(sb1))**(1.0d0/3.0d0)
      den=(1.0d0-nu)*sub2
      !expr11=8.0/3.0+((-3.359789466386329+5.8193260585515848*i)+(8.399473665965822-14.548315146289621*i)*nu)/den
      !expr12=(0.13228342099734997+0.22912160616643376*i)*sub2
      !expr1=expr11-expr12
      expr21=8.0d0/3.0d0+((-3.359789466386329d0-5.8193260585515848d0*i)+(8.399473665965822d0+14.548315146289621d0*i)*nu)/den
      expr22=(0.13228342099734997d0-0.22912160616643376d0*i)*sub2
      expr2=expr21-expr22
      gammas2=dreal(expr2)
      gammas=dsqrt(gammas2)
    else if((nu.gt.0.395d0).and.(nu.lt.0.405d0)) then
      gammas=dsqrt(0.887732d0)
    else
      sb1=-((-0.26308206488336366d0+nu)*(0.5939211404216125d0+(-0.23691793511663636d0+nu)*nu)/((-1.d0+nu)**3))
      sub2=(448.0d0-360.0d0/(1.0d0-nu)+235.1510153071851d0*cdsqrt(sb1))**(1.0d0/3.0d0)
      den=(1.0d0-nu)*sub2
      expr31=8.0/3.0+(6.719578932772657d0-16.798947331931643d0*nu)/den
      expr32=0.26456684199469993d0*sub2
      expr3=expr31+expr32
      gammas2=dreal(expr3)
      gammas=dsqrt(gammas2)
    end if
    return
  end subroutine obtenerc

  !! Body wavenumbers of elastic media
  subroutine harela_calculate_body_wavenumbers(lambda,mu,rho,omega,kp,ks)
    implicit none
    ! I/O
    complex(kind=real64) :: lambda, mu, rho
    real(kind=real64)    :: omega
    complex(kind=real64) :: kp, ks
    ! Local
    integer              :: i
    complex(kind=real64) :: k(2)
    ! Calculate kp
    k(1)= omega*sqrt(rho/(lambda+2.d0*mu))
    k(2)=-omega*sqrt(rho/(lambda+2.d0*mu))
    do i=1,2
      if (real(k(i)).gt.0) then
        kp=k(i)
        exit
      end if
    end do
    ! Calculate ks
    k(1)= omega*sqrt(rho/mu)
    k(2)=-omega*sqrt(rho/mu)
    do i=1,2
      if (real(k(i)).gt.0) then
        ks=k(i)
        exit
      end if
    end do
  end subroutine harela_calculate_body_wavenumbers

  !! Rayleigh surface wavenumber of elastic media
  subroutine harela_calculate_rayleigh_wavenumber(kp,ks,relative_error,kr)
    implicit none
    ! I/O
    complex(kind=real64) :: kp
    complex(kind=real64) :: ks
    real(kind=real64)    :: relative_error
    complex(kind=real64) :: kr
    ! Local
    complex(kind=real64) :: xi
    real(kind=real64)    :: nu, bvr
    integer              :: iteration
    complex(kind=real64) :: r_1, f_1, r_2, f_2, r_n, f_n, dfdr
    logical              :: next
    !
    ! Calculate a good approximation
    !
    ! Calculate the real part of nu
    xi=ks**2/kp**2
    nu=dreal(0.5d0*(xi-2.d0)/(xi-1.d0))
    ! Use the approximation of: A.V. Pichugin, Approximation of the Rayleigh wave speed, 2008 (draft submitted to Elsevier).
    ! bvr=cr/cs=ks/kr
    bvr=256.d0/293.d0+nu*(60.d0/307.d0-nu*(4.d0/125.d0+nu*(5.d0/84.d0+4.d0/237.d0*nu)))
    kr=ks/bvr
    !
    ! Refine it by using the characteristic equation solved by the Secant Method
    !
    r_1=0.95d0*kr/ks
    f_1=harela_rayleigh_equation(kp/ks,r_1)
    r_2=1.05d0*kr/ks
    f_2=harela_rayleigh_equation(kp/ks,r_2)
    iteration=0
    next=.true.
    do while (next)
      dfdr=(f_2-f_1)/(r_2-r_1)
      r_n=r_1-f_1/dfdr
      f_n=harela_rayleigh_equation(kp/ks,r_n)
      r_1=r_2
      f_1=f_2
      r_2=r_n
      f_2=f_n
      if (abs((r_1-r_2)/r_2).lt.relative_error) then
        next=.false.
      else
        iteration=iteration+1
        if (iteration.gt.100) then
          write(error_unit,*) 'Last 2 iterations [Re(r), Im(r), Re(f(r)), Im(f(r))]:'
          write(error_unit,*) dreal(r_1), dimag(r_1), dreal(f_1), dimag(f_1)
          write(error_unit,*) dreal(r_2), dimag(r_2), dreal(f_2), dimag(f_2)
          stop 'Error: Elastic Rayleigh wavenumber calculation has reached too many iterations.'
        end if
      end if
    end do
    kr=ks*r_n
  end subroutine harela_calculate_rayleigh_wavenumber

  !! Characteristic equation of Rayleigh surface waves for elastic solids
  function harela_rayleigh_equation(p,r)
    implicit none
    ! I/O
    complex(kind=real64) :: harela_rayleigh_equation !! F(r)
    complex(kind=real64) :: p                        !! kp/ks
    complex(kind=real64) :: r                        !! kr/ks
    ! Local
    complex(kind=real64) :: p2                       ! p^2
    complex(kind=real64) :: r2                       ! r^2
    ! Auxiliary values
    r2=r**2
    p2=p**2
    ! Compute
    harela_rayleigh_equation=r2**2-(sqrt(r2-1.d0)*sqrt(r2-p2)+1.d0)*r2+0.25d0
  end function harela_rayleigh_equation

  !! Calculate all elastic constants by giving two of them
  subroutine fbem_ela_properties(pE,E,pnu,nu,plambda,lambda,pmu,mu,pK,K)
    implicit none
    ! I/O
    logical           :: pE
    real(kind=real64) :: E
    logical           :: pnu
    real(kind=real64) :: nu
    logical           :: plambda
    real(kind=real64) :: lambda
    logical           :: pmu
    real(kind=real64) :: mu
    logical           :: pK
    real(kind=real64) :: K
    ! Local
    integer :: n
    ! Check input constants
    n=0
    if (plambda) n=n+1
    if (    pmu) n=n+1
    if (    pnu) n=n+1
    if (     pE) n=n+1
    if (     pK) n=n+1
    if (n.ne.2) stop 'only 2 elastic constants are needed'
    ! si se dan mas de 2 constantes, chequear que cumplen las relaciones adecuadas
    if (pK.and.pE) then
      lambda=3.d0*K*(3.d0*K-E)/(9.d0*K-E)
      mu=3.d0*K*E/(9.d0*K-E)
      nu=(3.d0*K-E)/(6.d0*K)
    end if
    if (pK.and.plambda) then
      E=9.d0*K*(K-lambda)/(3.d0*K-lambda)
      mu=1.5d0*(K-lambda)
      nu=lambda/(3.d0*K-lambda)
    end if
    if (pK.and.pmu) then
      E=9.d0*K*mu/(3.d0*K+mu)
      lambda=K-2.d0/3.d0*mu
      nu=0.5d0*(3.d0*K-2.d0*mu)/(3.d0*K+mu)
    end if
    if (pK.and.pnu) then
      E=3.d0*K*(1.d0-2.d0*nu)
      lambda=3.d0*K*nu/(1.d0+nu)
      mu=1.5d0*K*(1.d0-2.d0*nu)/(1.d0+nu)
    end if
    if (pE.and.plambda) then
      K=(E+3.d0*lambda+sqrt(E**2+9.d0*lambda**2+2.d0*E*lambda))/6.d0
      mu=(E-3.d0*lambda+sqrt(E**2+9.d0*lambda**2+2.d0*E*lambda))/4.d0
      nu=2.d0*lambda/(E+lambda+sqrt(E**2+9.d0*lambda**2+2.d0*E*lambda))
    end if
    if (pE.and.pmu) then
      K=E*mu/(3.d0*mu-E)/3.d0
      lambda=mu*(E-2.d0*mu)/(3.d0*mu-E)
      nu=0.5d0*E/mu-1.d0
    end if
    if (pE.and.pnu) then
      K=E/(1.d0-2.d0*nu)/3.d0
      lambda=E*nu/(1.d0+nu)/(1.d0-2.d0*nu)
      mu=0.5d0*E/(1.d0+nu)
    end if
    if (plambda.and.pmu) then
      K=lambda+2.d0*mu/3.d0
      E=mu*(3.d0*lambda+2.d0*mu)/(lambda+mu)
      nu=0.5d0*lambda/(lambda+mu)
    end if
    if (plambda.and.pnu) then
      K=lambda*(1.d0+nu)/(3.d0*nu)
      E=lambda/nu*(1.d0+nu)*(1.d0-2.d0*nu)
      mu=0.5d0*lambda*(1.d0-2.d0*nu)/nu
    end if
    if (pmu.and.pnu) then
      K=2.d0*mu*(1.d0+nu)/(1.d0-2.d0*nu)/3.d0
      E=2.d0*mu*(1.d0+nu)
      lambda=2.d0*mu*nu/(1.d0-2.d0*nu)
    end if
    ! Check if constants are valid
    if (mu.le.0.d0) stop 'mu must be >0'
    if ((3.d0*lambda+2.d0*mu).le.0.d0) stop '3*lambda+2*mu must be >0'
  end subroutine fbem_ela_properties

end module fbem_harela_incident_field
