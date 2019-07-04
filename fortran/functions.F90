! Copyright 2017-2018 Tristan Salles
!
! This file is part of eSCAPE.
!
! eSCAPE is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or any later version.
!
! eSCAPE is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with eSCAPE.  If not, see <http://www.gnu.org/licenses/>.

! f2py3 --overwrite-signature -m _fortran -h functions.pyf functions.f90

#include "petsc/finclude/petsc.h"
#include "petsc/finclude/petscmat.h"

#include "petscversion.h"

#undef  CHKERRQ
#define CHKERRQ(n) if ((n) .ne. 0) return;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INTERNAL FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module meshparams

  implicit none

  integer :: nLocal
  integer, dimension(:,:), allocatable :: FVnID
  integer, dimension(:,:), allocatable :: QinID
  integer, dimension(:), allocatable :: FVnNb

  real( kind=8 ), dimension(:), allocatable :: FVarea
  real( kind=8 ), dimension(:,:), allocatable :: FVeLgt
  real( kind=8 ), dimension(:,:), allocatable :: FVvDist


end module meshparams

subroutine euclid( p1, p2, norm)
!*****************************************************************************
! Computes the Euclidean vector norm between 2 points
! along the 3 dimensions

  implicit none

  real( kind=8 ), intent(in) :: p1(3)
  real( kind=8 ), intent(in) :: p2(3)
  real( kind=8 ), intent(out) :: norm
  real( kind=8 ) :: vec(3)

  vec = p2 - p1
  norm = norm2(vec)

  return

end subroutine euclid

recursive subroutine quicksort(array, first, last, indices)
!*****************************************************************************
! quicksort routine from http://www.cgd.ucar.edu/pubsoft/TestQuicksort.html
! Reference:
! Nyhoff & Leestma, Fortran 90 for Engineers & Scientists
! (New Jersey: Prentice Hall, 1997), pp 575-577.

  real( kind=8 ), dimension(:), intent(inout) :: array
  integer, intent(in)  :: first, last
  integer, dimension(:), intent(inout) :: indices

  interface
       subroutine split(array, low, high, mid, indices)
          real( kind=8 ), dimension(:), intent(inout) :: array
          integer, intent(in) :: low, high
          integer, intent(out) :: mid
          integer, dimension(:), intent(inout) :: indices
       end subroutine split
  end interface

  integer :: mid

  if(first < last)then
    call split(array, first, last, mid, indices)
    call quicksort(array, first, mid-1, indices)
    call quicksort(array, mid+1, last,  indices)
  endif

end subroutine quicksort

subroutine split(array, low, high, mid, indices)
!*****************************************************************************
! used by quicksort  from http://www.cgd.ucar.edu/pubsoft/TestQuicksort.html
! Reference:
! Nyhoff & Leestma, Fortran 90 for Engineers & Scientists
! (New Jersey: Prentice Hall, 1997), pp 575-577.

  real( kind=8 ), dimension(:), intent(inout) :: array
  integer, intent(in) :: low, high
  integer, intent(out) :: mid
  integer, dimension(:), intent(inout) :: indices

  integer :: left, right
  real( kind=8 ) ::  pivot, swap
  integer :: ipivot, iswap

  left = low
  right = high
  pivot = array(low)
  ipivot = indices(low)

  do
    if( left >= right ) exit
    do
      if( left >= right .or. array(right) < pivot ) exit
      right = right - 1
    enddo
    do
      if(left > 12) exit
      if(array(left) > pivot) exit
      left = left + 1
    enddo

    if( left < right )then
      swap  = array(left)
      array(left)  = array(right)
      array(right) = swap
      iswap = indices(left)
      indices(left)  = indices(right)
      indices(right) = iswap
    endif
  enddo

  array(low) = array(right)
  array(right) = pivot
  mid = right
  indices(low) = indices(right)
  indices(right) = ipivot

end subroutine split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! BOUNDARY CONDITION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine meanSlope(elev, bID, gbounds, bslp, nb, b)
!*****************************************************************************
! Define average slope along boundary edges

    use meshparams
    implicit none

    integer :: b,nb
    integer, intent(in) :: bID(b)
    integer, intent(in) :: gbounds(nb)
    real( kind=8 ), intent(in) :: elev(nb)

    real( kind=8 ), intent(out) :: bslp(nb)

    integer :: k, p, n, nn, kk
    real( kind=8 ) :: smean


    bslp = 0.
    do k = 1, b
      n = bID(k)+1
      smean = 0.
      kk = 0
      do p = 1, FVnNb(n)
        nn = FVnID(n,p)+1
        if(gbounds(nn)==1 .and. FVeLgt(n,p)>0)then
          smean = smean+(elev(n)-elev(nn))/FVeLgt(n,p)
          kk = kk+1
        endif
      enddo
      if(kk>0) bslp(n) = smean/kk
    enddo

    return

end subroutine meanSlope

subroutine flatBounds(elev, erodep, bID, gbounds, be, bd, nb, b)
!*****************************************************************************
! Define flat boundary conditions

    use meshparams
    implicit none

    integer :: b,nb
    integer, intent(in) :: bID(b)
    integer, intent(in) :: gbounds(nb)
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: erodep(nb)

    real( kind=8 ), intent(out) :: be(nb)
    real( kind=8 ), intent(out) :: bd(nb)

    integer :: k, p, n, nn, kk
    real( kind=8 ) :: esum, bsum

    be = elev
    bd = erodep
    do k = 1, b
      n = bID(k)+1
      esum = 0.
      bsum = 0.
      kk = 0
      do p = 1, FVnNb(n)
        nn = FVnID(n,p)+1
        if(gbounds(nn)==1 .and. FVeLgt(n,p)>0)then
          esum = esum+elev(nn)
          bsum = bsum+erodep(nn)
          kk = kk + 1
        endif
      enddo
      if(kk>0)then
        be(n) = esum/kk
        bd(n) = bsum/kk
      endif
    enddo

    return

end subroutine flatBounds

subroutine slpBounds(elev, erodep, bID, gbounds, be, bd, nb, b)
!*****************************************************************************
! Define slope boundary conditions

    use meshparams
    implicit none

    integer :: b,nb
    integer, intent(in) :: bID(b)
    integer, intent(in) :: gbounds(nb)
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: erodep(nb)

    real( kind=8 ), intent(out) :: be(nb)
    real( kind=8 ), intent(out) :: bd(nb)

    integer :: k, p, n, nn
    real( kind=8 ) :: emin, bmin

    be = elev
    bd = erodep
    do k = 1, b
      n = bID(k)+1
      emin = 1.e8
      bmin = 1.e8
      do p = 1, FVnNb(n)
        nn = FVnID(n,p)+1
        if(gbounds(nn)==1.and.FVeLgt(n,p)>0.)then
          emin = min(emin,elev(nn))
          bmin = min(bmin,erodep(nn))
        endif
      enddo
      if(emin<1.e8) be(n) = emin-1.e-12
      if(bmin<1.e8) bd(n) = bmin-1.e-12
    enddo

    return

end subroutine slpBounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HILLSLOPE PROCESSES FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine setHillslopeCoeff(nb, Kd, dcoeff)
!*****************************************************************************
! Define hillslope coefficients

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: Kd

    real( kind=8 ), intent(out) :: dcoeff(nb,13)

    integer :: k, p
    real( kind=8 ) :: s1, c, v

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        c = Kd/FVarea(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            s1 = s1 + v
            dcoeff(k,p+1) = -v
          endif
        enddo
        dcoeff(k,1) = 1.0 + s1
      endif
    enddo

    return

end subroutine setHillslopeCoeff

subroutine setDiffusionCoeff(Kd, limit, elev, elev0, dh, dcoeff, nb)
!*****************************************************************************
! Define freshly deposited sediments diffusion implicit matrix coefficients

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: Kd
    real( kind=8 ), intent(in) :: limit
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: elev0(nb)
    real( kind=8 ), intent(in) :: dh(nb)

    real( kind=8 ), intent(out) :: dcoeff(nb,13)

    integer :: k, p, n
    real( kind=8 ) :: s1, c, v, limiter

    dcoeff = 0.
    do k = 1, nb
      s1 = 0.
      if(FVarea(k)>0)then
        c = Kd/FVarea(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0.)then
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            n = FVnID(k,p)+1
            limiter = 0.
            if(elev(n)>elev(k) .and. elev0(n)<0.)then
              limiter = dh(n)/(dh(n)+limit)
            elseif(elev(n)<elev(k) .and. elev0(k)<0.)then ! elev(k)>elev0(k))then
              limiter = dh(k)/(dh(k)+limit)
            endif
            s1 = s1 + v*limiter
            dcoeff(k,p+1) = -v*limiter
          endif
        enddo
        dcoeff(k,1) = 1.0 + s1
      endif
    enddo

    return

end subroutine setDiffusionCoeff

subroutine explicitDiff(Kd, limit, elev, elev0, dh, newz, nb)
!*****************************************************************************
! Define freshly deposited sediments diffusion explicitly

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: Kd
    real( kind=8 ), intent(in) :: limit
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: elev0(nb)
    real( kind=8 ), intent(in) :: dh(nb)

    real( kind=8 ), intent(out) :: newz(nb)

    integer :: k, p, n, i, pp
    real( kind=8 ) :: c, v, Qsout, zm, dhdt_dt
    real( kind=8 ) :: SumQs(nb), sf(nb), Qsin(nb,12)


    newz = elev
    sf = 0.
    SumQs = 0.
    Qsin = 0.

    do k = 1, nb
      if(FVarea(k)>0 .and. dh(k)>0. .and. elev0(k)<=0.)then
        c = Kd/FVarea(k)
        i = 0
        zm = elev(k)
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          v = c*FVvDist(k,p)/FVeLgt(k,p)
          if(elev(k)>elev(n))then
            Qsout = v*(elev(k)-elev(n))
            pp = QinID(k,p)
            Qsin(n,pp) = Qsin(n,pp) + Qsout
            SumQs(k) = SumQs(k) + Qsout
            zm = min(zm,elev(n))
          endif
        enddo
        if(SumQs(k)>0.)then
          sf(k) = min(dh(k),limit*(elev(k)-zm))
          sf(k) = sf(k)/sumQs(k)
          sf(k) = min(sf(k),1.)
        endif
      endif
    enddo

    do k = 1, nb
      dhdt_dt = 0.
      if(FVarea(k)>0 .and. elev0(k)<=0.)then
        dhdt_dt = -sf(k)*SumQs(k)
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          dhdt_dt = dhdt_dt + sf(n)*Qsin(k,p)
        enddo
        newz(k) = newz(k) + dhdt_dt
      endif
    enddo

    return

end subroutine explicitDiff

subroutine explicitDiffOld(Kd, limit, elev, elev0, dh, newz, nb)
!*****************************************************************************
! Define freshly deposited sediments diffusion explicitly

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: Kd
    real( kind=8 ), intent(in) :: limit
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: elev0(nb)
    real( kind=8 ), intent(in) :: dh(nb)

    real( kind=8 ), intent(out) :: newz(nb)

    integer :: k, p, n
    real( kind=8 ) :: vals !, maxe
    real( kind=8 ) :: s1, c, v, limiter

    newz = elev
    do k = 1, nb
      s1 = 0.
      vals = 0.
      if(FVarea(k)>0)then
        c = Kd/FVarea(k)
        do p = 1, FVnNb(k)
          ! maxe = 0.
          if(FVvDist(k,p)>0.)then
            v = c*FVvDist(k,p)/FVeLgt(k,p)
            n = FVnID(k,p)+1
            limiter = 0.
            if(elev(n)>elev(k) .and. elev0(n)<0.)then
              limiter = dh(n)/(dh(n)+limit)
              ! maxe = 0.1*dh(n)
            elseif(elev(n)<elev(k) .and. elev0(k)<0.)then
              limiter = dh(k)/(dh(k)+limit)
              ! maxe = 0.1*dh(k)
            endif
            s1 = s1 + v*limiter
            vals = vals + v*limiter*elev(n)
          endif
        enddo
        newz(k) = elev(k)*(1.0 - s1) + vals
      endif
    enddo

    return

end subroutine explicitDiffOld

subroutine distributeHeight( inIDs, sl, elev, elev0, sed, nelev, nsed, nb )
!*****************************************************************************
! Update elevation based on incoming sedimentary volume

  use meshparams
  implicit none

  integer :: nb

  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: sl
  real( kind=8 ), intent(in) :: elev(nb)
  real( kind=8 ), intent(in) :: elev0(nb)
  real( kind=8 ), intent(in) :: sed(nb)

  real( kind=8 ), intent(out) :: nelev(nb)
  real( kind=8 ), intent(out) :: nsed(nb)

  integer :: k,p,n,nbr
  real( kind=8 ) :: hmax,dh,dv

  nelev = elev
  nsed = sed
  do k = 1, nb
    if(inIDs(k)>0)then
      ! Deposition on the local node
      if(sed(k)>0)then
        nbr = 0
        hmax = 1.e8
        if(elev0(k)<sl)  hmax = elev(k)
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          if(elev0(k)<sl)then
            if(elev(n)<sl .and. elev(n)>elev(k))then
              hmax = max(elev(n)-1.e-3,hmax)
            endif
            if(elev(n)<elev(k)) nbr = nbr+1
          else
            if(elev(n)>elev(k))  hmax = min(hmax,elev(n)-1.e-3)
            if(elev(n)<elev(k)) nbr = nbr+1
          endif
        enddo
        if(elev0(k)<sl .and. hmax==elev(k) .and. nbr == 0) hmax = hmax + 1.e-4
        if(elev0(k)>=sl .and. hmax==elev(k) .and. nbr == 0) hmax = hmax + 1.e-4
        if(elev0(k)<sl) hmax = min(hmax,0.1*(sl-elev(k))+elev(k))
        if(elev(k)<hmax)then
          dh = (hmax-elev(k))
          dv = dh*FVarea(k)
          if(dv>sed(k))then
            nelev(k) = elev(k)+sed(k)/FVarea(k)
            nsed(k) = 0.
          else
            nelev(k) = elev(k)+dh
            nsed(k) = sed(k) - dv
          endif
        endif
      endif
    endif
  enddo

  return

end subroutine distributeHeight

subroutine distributeVolume( inIDs, sl, elev, elev0, sed, nsed, nb )
!*****************************************************************************
! Compute remaining volume to distribute

  use meshparams
  implicit none

  integer :: nb

  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: sl
  real( kind=8 ), intent(in) :: elev(nb)
  real( kind=8 ), intent(in) :: elev0(nb)
  real( kind=8 ), intent(in) :: sed(nb)

  real( kind=8 ), intent(out) :: nsed(nb)

  integer :: k,p,n,nbr,num
  real( kind=8 ) :: prop

  nsed = 0.
  do k = 1, nb
    if(inIDs(k)>=0)then
      if(sed(k)>0)then
        nbr = 0
        num = 0
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          if(elev0(k)<sl .and. elev(n) < sl) num = num+1
          if(elev(n)<elev(k)) nbr = nbr+1
        enddo
        if(nbr>0)then
          prop = sed(k)/nbr
          do p = 1, FVnNb(k)
            n = FVnID(k,p)+1
            if(elev(n)<elev(k)) nsed(n) = nsed(n) + prop
          enddo
        else
          if(elev0(k)<sl .and. num>0)then
            prop = sed(k)/(num+1)
          else
            nbr = FVnNb(k)+1
            prop = sed(k)/nbr
          endif
          nsed(k) = nsed(k) + prop
          do p = 1, FVnNb(k)
            n = FVnID(k,p)+1
            if(num>0)then
              if(elev(n)<sl) nsed(n) = nsed(n) + prop
            else
              nsed(n) = nsed(n) + prop
            endif
          enddo
        endif
      endif
    endif
  enddo

  return

end subroutine distributeVolume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SEDIMENT DIFFUSION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine diffusionDT(dm, hLocal, hL0, bounds, iters, itflx, inIDs, sFlux, sKd, &
                       oKd, sl, ierr, nb)
!*****************************************************************************
! Internal loop for marine and aerial diffusion on currently deposited sediments.

  use petsc
  use petscmat
  use meshparams
  implicit none

  integer :: nb
  DM, intent(in) :: dm
  integer, intent(in) :: iters
  integer, intent(in) :: itflx

  real( kind=8 ), intent(in) :: sFlux(nb)
  real( kind=8 ), intent(in) :: sKd(nb,12)
  real( kind=8 ), intent(in) :: oKd(nb,12)
  real( kind=8 ), intent(in) :: sl
  integer, intent(in) :: inIDs(nb)
  integer, intent(in) :: bounds(nb)

  Vec, intent(inout) :: hLocal
  Vec, intent(in) :: hL0
  PetscErrorCode,intent(out) :: ierr

  PetscScalar, pointer :: hArr0(:)
  PetscScalar, pointer :: hArr(:)
  PetscScalar, pointer :: Cero(:)
  real( kind=8 ) :: dh(nb)
  Vec :: vLoc
  Vec :: vGlob

  integer :: tstep
  integer :: k, n, p
  real( kind=8 ) :: kd, val0, val1(nb,12)

  ! Define vectors
  call DMCreateGlobalVector(dm,vGlob,ierr);CHKERRQ(ierr)
  call VecDuplicate(hLocal,vLoc,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(hL0,hArr0,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(hLocal,hArr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(vLoc,Cero,ierr)

  do tstep = 1, iters

    dh = 0.
    val1 = 0.

    ! Compute maximum erosion thickness for diffusion to ensure stability
    do k = 1, nb
      Cero(k) = 1.0
      val0 = 0.
      if(inIDs(k)>=0)then
        do p = 1, FVnNb(k)
          n = FVnID(k,p)+1
          kd = sKd(k,p)
          if(n>0 .and. FVeLgt(k,p)>0.)then
            if(hArr(k)<sl .and. hArr(n)<sl)then
              kd = oKd(k,p)
            elseif(hArr(k)>=sl .and. hArr(n)>=sl)then
              kd = sKd(k,p)
            elseif(hArr(k)<sl .and. hArr(n)>=sl)then
              kd = sKd(k,p)
            elseif(hArr(k)>=sl .and. hArr(n)<sl)then
              kd = sKd(k,p)
            else
              if(hArr(k)>=sl) kd = sKd(k,p)
              if(hArr(k)<sl) kd = oKd(k,p)
            endif
            if(hArr(k)>hArr(n))then
              val0 = val0+kd*(hArr(n)-hArr(k))
            endif
            val1(k,p) = kd*(hArr(n)-hArr(k))
          endif
        enddo
      endif
      if(val0<0 .and. hArr(k)>hArr0(k))then
        if(val0<hArr0(k)-hArr(k)) Cero(k) = (hArr0(k)-hArr(k))/val0
      elseif(hArr(k)==hArr0(k))then
        Cero(k) = 0.
      endif
    enddo
    call VecRestoreArrayF90(vLoc,Cero,ierr)

    ! Update ghosts
    call DMLocalToGlobalBegin(dm,vLoc,INSERT_VALUES,vGlob,ierr)
    call DMLocalToGlobalEnd(dm,vLoc,INSERT_VALUES,vGlob,ierr)
    call DMGlobalToLocalBegin(dm,vGlob,INSERT_VALUES,vLoc,ierr)
    call DMGlobalToLocalEnd(dm,vGlob,INSERT_VALUES,vLoc,ierr)

    call VecGetArrayF90(vLoc,Cero,ierr)

    ! Compute elevation change due to diffusion
    do k = 1, nb
      if(inIDs(k)>0 .and. bounds(k)==0)then
        do p = 1, FVnNb(k)
          if(val1(k,p)>0.)then
            n = FVnID(k,p)+1
            dh(k) = dh(k) + val1(k,p)*Cero(n)
          elseif(val1(k,p)<0.)then
            dh(k) = dh(k) + val1(k,p)*Cero(k)
          endif
        enddo
      endif
    enddo
    hArr = hArr + dh
    if(tstep < itflx) hArr = hArr + sFlux

    call VecRestoreArrayF90(hLocal,hArr,ierr)
    call DMLocalToGlobalBegin(dm,hLocal,INSERT_VALUES,vGlob,ierr)
    call DMLocalToGlobalEnd(dm,hLocal,INSERT_VALUES,vGlob,ierr)
    call DMGlobalToLocalBegin(dm,vGlob,INSERT_VALUES,hLocal,ierr)
    call DMGlobalToLocalEnd(dm,vGlob,INSERT_VALUES,hLocal,ierr)

    ! Update upper layer fraction
    call VecGetArrayF90(hLocal,hArr,ierr)
  enddo

  call VecRestoreArrayF90(vLoc,Cero,ierr)
  call VecRestoreArrayF90(hLocal,hArr,ierr)
  call VecRestoreArrayF90(hL0,hArr0,ierr)

  call VecDestroy(vGlob,ierr)
  call VecDestroy(vLoc,ierr)

  return

end subroutine diffusionDT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FLOW DIRECTION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MFDreceivers( nRcv, inIDs, elev, rcv, slope, dist, wgt, nb)
!*****************************************************************************
! Compute receiver characteristics based on multiple flow direction algorithm

  use meshparams
  implicit none

  interface
    recursive subroutine quicksort(array, first, last, indices)
      real( kind=8 ), dimension(:), intent(inout) :: array
      integer, intent(in)  :: first, last
      integer, dimension(:), intent(inout) :: indices
    end subroutine quicksort
  end interface

  integer :: nb

  integer, intent(in) :: nRcv
  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: elev(nb)

  integer, intent(out) :: rcv(nb,nRcv)
  real( kind=8 ), intent(out) :: slope(nb,nRcv)
  real( kind=8 ), intent(out) :: dist(nb,nRcv)
  real( kind=8 ), intent(out) :: wgt(nb,nRcv)

  integer :: k, n, p, kk
  real( kind=8 ) :: slp(12),dst(12),val
  integer :: id(12)

  rcv = -1
  slope = 0.
  dist = 0.
  wgt = 0.

  do k = 1, nb
    if(inIDs(k)>0)then
      slp = 0.
      id = 0
      val = 0.
      kk = 0
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        if(n>0 .and. FVeLgt(k,p)>0.)then
          val = (elev(k) - elev(n))/FVeLgt(k,p)
          if(val>0.)then
            kk = kk + 1
            slp(kk) = val
            id(kk) = n-1
            dst(kk) = FVeLgt(k,p)
          endif
        endif
      enddo

      if(kk == 0)then
        rcv(k,1:nRcv) = k-1
      elseif(kk <= nRcv)then
        val = 0.
        rcv(k,1:nRcv) = k-1
        do p = 1, kk
          rcv(k,p) = id(p)
          slope(k,p) = slp(p)
          dist(k,p) = dst(p)
          val = val + slp(p)
        enddo
        do p = 1, nRcv
          wgt(k,p) = slp(p)/val
        enddo
      else
        rcv(k,1:nRcv) = k-1
        call quicksort(slp,1,kk,id)
        n = 0
        val = 0.
        do p = kk,kk-nRcv+1,-1
          n = n + 1
          rcv(k,n) = id(p)
          slope(k,n) = slp(p)
          dist(k,n) = dst(p)
          val = val + slp(p)
        enddo
        do p = 1, nRcv
          wgt(k,p) = slope(k,p)/val
        enddo
      endif
    endif
  enddo

  return

end subroutine MFDreceivers

subroutine minHeight( inIDs, elev, hmin, nb)
!*****************************************************************************
! Compute minimum height of donor nodes

  use meshparams
  implicit none

  integer :: nb

  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: elev(nb)

  real( kind=8 ), intent(out) :: hmin(nb)

  integer :: k, n, p

  hmin = 0.
  do k = 1, nb
    if(inIDs(k)>=0)then
      hmin(k) = 1.e8
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        if(n>0 .and. FVeLgt(k,p)>0.)then
          if(elev(n)>elev(k) .and. elev(n)-elev(k)<hmin(k))then
            hmin(k) = elev(n)-elev(k)
          endif
        endif
      enddo
      if(hmin(k) > 9.9e7) hmin(k) = 0.
    endif
  enddo

  return

end subroutine minHeight

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PIT FILLING FUNCTION !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pitData( rank, pnode, inids, natural, tmp3, nb, nb1, nb2 )
!*****************************************************************************
! Compute pit information: processor rank and local point ID

  implicit none

  integer :: nb, nb1, nb2

  integer, intent(in) :: rank
  integer, intent(in) :: pnode(nb)
  integer, intent(in) :: inids(nb1)
  integer, intent(in) :: natural(nb2)

  integer, intent(out) :: tmp3(nb*2)

  integer :: k, p

  tmp3 = -1

  do k = 1, nb
    loop: do p = 1, nb2
      if( natural(p) == pnode(k) )then
        if( inids(p) == 1 )then
          tmp3(k) = rank
          tmp3(k+nb) = p-1
          exit loop
        endif
      endif
    enddo loop
  enddo

  return

end subroutine pitData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MESH DECLARATION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine defineTIN( coords, cells_nodes, cells_edges, edges_nodes, area, circumcenter, &
                                    ngbNb, ngbID, edgeLgt, voroDist, n, nb, m  )
!*****************************************************************************
! Compute for a specific triangulation the characteristics of each node and
! associated voronoi for finite volume discretizations

  use meshparams
  implicit none

  integer :: m, n, nb
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: cells_edges(n,3)
  integer, intent(in) :: edges_nodes(m, 2)

  real( kind=8 ), intent(in) :: coords(nb,3)
  real( kind=8 ), intent(in) :: area(nb)
  real( kind=8 ), intent(in) :: circumcenter(3,n)

  integer, intent(out) :: ngbID(nb, 12)
  integer, intent(out) :: ngbNb(nb)

  real( kind=8 ), intent(out) :: edgeLgt(nb,12)
  real( kind=8 ), intent(out) :: voroDist(nb,12)

  integer :: i, n1, n2, k, l, p, eid, cid, e, id
  integer :: nid(2), nc(3), edge(nb, 12)
  integer :: edgeNb(3), edges(3,2), cell_ids(nb, 12)

  real( kind=8 ) :: coords0(3), coordsID(3)
  real( kind=8 ) :: midpoint(3), dist

  logical :: inside

  cell_ids(:,:) = -1
  edge(:,:) = -1
  ngbNb(:) = 0
  ngbID(:,:) = -1
  edgeLgt(:,:) = 0.

  ! Find all cells surrounding a given vertice
  do i = 1, n
    nc = cells_nodes(i,1:3)+1
    do p = 1, 3
      inside = .False.
      lp: do k = 1, 12
        if( cell_ids(nc(p),k) == i-1 )then
          exit lp
        elseif( cell_ids(nc(p),k) == -1 )then
          inside = .True.
          exit lp
        endif
      enddo lp
      if( inside )then
        cell_ids(nc(p),k)  = i-1
      endif
    enddo
  enddo

  ! Find all edges connected to a given vertice
  do i = 1, m
    n1 = edges_nodes(i,1)+1
    n2 = edges_nodes(i,2)+1
    inside = .False.
    lp0: do k = 1, 12
      if(edge(n1,k) == i-1)then
        exit lp0
      elseif(edge(n1,k) == -1)then
        inside = .True.
        exit lp0
      endif
    enddo lp0
    if( inside )then
      edge(n1,k)  = i-1
      ngbNb(n1) = ngbNb(n1) + 1
    endif
    inside = .False.
    lp1: do k = 1, 12
      if(edge(n2,k) == i-1)then
        exit lp1
      elseif(edge(n2,k) == -1)then
        inside = .True.
        exit lp1
      endif
    enddo lp1
    if( inside )then
      edge(n2,k)  = i-1
      ngbNb(n2) = ngbNb(n2) + 1
    endif
  enddo

  do k = 1, nb

    ! Get triangulation edge lengths
    coords0 = coords(k,1:3)
    l = 0
    do eid = 1, ngbNb(k)
      nid = edges_nodes(edge(k,eid)+1,1:2)
      if( nid(1) == k-1)then
        l = l + 1
        ngbID(k,l) = nid(2)
        coordsID = coords(nid(2)+1,1:3)
      else
        l = l + 1
        ngbID(k,l) = nid(1)
        coordsID = coords(nid(1)+1,1:3)
      endif
      call euclid( coords0, coordsID, edgeLgt(k,l) )
    enddo

    ! Get voronoi edge lengths
    lp2: do cid = 1, 12
      if( cell_ids(k,cid) == -1 ) exit lp2
      edgeNb(1:3) = cells_edges( cell_ids(k,cid)+1,1:3 )
      do e = 1, 3
        edges(e,1:2) = edges_nodes(edgeNb(e)+1,1:2)
        if( k-1 == edges(e,1) .or. k-1 == edges(e,2))then
          midpoint(1:3) = 0.5 * (coords(edges(e,1)+1,1:3)+coords(edges(e,2)+1,1:3))
          id = -1
          if( edges(e,1) == k-1 )then
            lp3: do i = 1, ngbNb(k)
              if(ngbID(k,i) == edges(e,2))then
                id = i
                exit lp3
              endif
            enddo lp3
          else
            lp4: do i = 1, ngbNb(k)
              if(ngbID(k,i) == edges(e,1))then
                id = i
                exit lp4
              endif
            enddo lp4
          endif
          call euclid( midpoint(1:3), circumcenter(1:3,cell_ids(k,cid)+1),  dist)
          voroDist(k,id) = voroDist(k,id) + dist
        endif
      enddo
    enddo lp2

  enddo

  ! Store mesh parameters
  nLocal = nb
  if(allocated(FVarea)) deallocate(FVarea)
  if(allocated(FVnID)) deallocate(FVnID)
  if(allocated(FVnNb)) deallocate(FVnNb)
  if(allocated(FVeLgt)) deallocate(FVeLgt)
  if(allocated(FVvDist)) deallocate(FVvDist)
  if(allocated(QinID)) deallocate(QinID)

  allocate(FVarea(nLocal))
  allocate(FVnNb(nLocal))
  allocate(FVnID(nLocal,12))
  allocate(FVeLgt(nLocal,12))
  allocate(FVvDist(nLocal,12))
  allocate(QinID(nLocal,12))

  FVarea = area
  FVnNb = ngbNb
  FVnID = ngbID
  FVeLgt = edgeLgt
  FVvDist = voroDist

  QinID = -1
  do k = 1, nLocal
    do p = 1, FVnNb(k)
      l = FVnID(k,p)+1
      lpp: do e = 1, FVnNb(l)
        if(FVnID(l,e)+1 == k)then
          QinID(k,p) = e
          exit lpp
        endif
      enddo lpp
    enddo
  enddo

end subroutine defineTIN


subroutine defineGTIN( nb, cells_nodes, edges_nodes, ngbNb, ngbID, n, m)
!*****************************************************************************
! Compute for global triangulation the characteristics of each node

  use meshparams
  implicit none

  integer :: m, n
  integer, intent(in) :: nb
  integer, intent(in) :: cells_nodes(n, 3)
  integer, intent(in) :: edges_nodes(m, 2)

  integer, intent(out) :: ngbID(nb, 12)
  integer, intent(out) :: ngbNb(nb)

  integer :: i, n1, n2, k, l, eid, p
  integer :: nid(2), nc(3), edge(nb, 12)
  integer :: cell_ids(nb, 12)

  logical :: inside

  cell_ids(:,:) = -1
  edge(:,:) = -1
  ngbNb(:) = 0
  ngbID(:,:) = -1

  ! Find all cells surrounding a given vertice
  do i = 1, n
    nc = cells_nodes(i,1:3)+1
    do p = 1, 3
      inside = .False.
      lp: do k = 1, 12
        if( cell_ids(nc(p),k) == i-1 )then
          exit lp
        elseif( cell_ids(nc(p),k) == -1 )then
          inside = .True.
          exit lp
        endif
      enddo lp
      if( inside )then
        cell_ids(nc(p),k)  = i-1
      endif
    enddo
  enddo

  ! Find all edges connected to a given vertice
  do i = 1, m
    n1 = edges_nodes(i,1)+1
    n2 = edges_nodes(i,2)+1
    inside = .False.
    lp0: do k = 1, 12
      if(edge(n1,k) == i-1)then
        exit lp0
      elseif(edge(n1,k) == -1)then
        inside = .True.
        exit lp0
      endif
    enddo lp0
    if( inside )then
      edge(n1,k)  = i-1
      ngbNb(n1) = ngbNb(n1) + 1
    endif
    inside = .False.
    lp1: do k = 1, 12
      if(edge(n2,k) == i-1)then
        exit lp1
      elseif(edge(n2,k) == -1)then
        inside = .True.
        exit lp1
      endif
    enddo lp1
    if( inside )then
      edge(n2,k)  = i-1
      ngbNb(n2) = ngbNb(n2) + 1
    endif
  enddo

  do k = 1, nb
    ! Get triangulation edge lengths
    l = 0
    do eid = 1, ngbNb(k)
      nid = edges_nodes(edge(k,eid)+1,1:2)
      if( nid(1) == k-1)then
        l = l + 1
        ngbID(k,l) = nid(2)
      else
        l = l + 1
        ngbID(k,l) = nid(1)
      endif
    enddo

  enddo

end subroutine defineGTIN
