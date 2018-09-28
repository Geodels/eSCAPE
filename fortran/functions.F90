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

! f2py --overwrite-signature -m _fortran -h functions.pyf functions.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INTERNAL FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module meshparams

  implicit none

  integer :: nLocal
  integer, dimension(:,:), allocatable :: FVnID
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PIT FILLING RELATED FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine pitVolume(depLocal, pitID, pitNb, pitVol, notpitVol, m)
!*****************************************************************************
! Get local volume of sediment deposited on each pit for unstructured grids

  implicit none

  integer :: m
  integer, intent(in) :: pitNb
  real( kind=8 ), intent(in) :: depLocal(m)
  integer, intent(in) :: pitID(m)

  real( kind=8 ), intent(out) :: notpitVol(m)
  real( kind=8 ), intent(out) :: pitVol(pitNb)

  integer :: k, c

  pitVol = 0.
  notpitVol = 0.

  do k = 1, m
    if(depLocal(k)>0.)then
      c = pitID(k)
      if(c>0)then
        pitVol(c) = pitVol(c)+depLocal(k)
      else
        notpitVol(k) = depLocal(k)
      endif
    endif
  enddo

  return

end subroutine pitVolume

subroutine pitHeight(elev, fillZ, pitID, pitVol, pitsedVol, newZ, remain, totnodes, m, nn, nb)
!*****************************************************************************
! Update elevation each pit

  implicit none

  integer :: m, nn, nb
  real( kind=8 ), intent(in) :: elev(m)
  real( kind=8 ), intent(in) :: fillZ(m)
  integer, intent(in) :: pitID(m)
  real( kind=8 ), intent(in) :: pitVol(nn)
  real( kind=8 ), intent(in) :: pitsedVol(nb)

  real( kind=8 ), intent(out) :: newZ(m)
  real( kind=8 ), intent(out) :: remain(nb)
  integer, intent(out) :: totnodes(nb)

  integer :: k, c
  real( kind=8 ) :: frac

  newZ = elev
  remain = 0.
  totnodes = 0

  ! Update elevation
  do k = 1, m
    c = pitID(k)
    if(c>0)then
      if(pitsedVol(c)>=pitVol(c))then
        newZ(k) = fillZ(k)
        remain(c) = pitsedVol(c)-pitVol(c)
        totnodes(c) = totnodes(c) + 1
      else
        if(pitVol(c)>0.)then
          frac = pitsedVol(c)/pitVol(c)
          newZ(k) = frac*(fillZ(k)-elev(k))+elev(k)
        endif
        remain(c) = 0.
      endif
    endif
  enddo

  return

end subroutine pitHeight

subroutine addExcess(excess, pitID, addVol, m, n)
!*****************************************************************************
! Add excess sediment volume on pit nodes

  implicit none

  integer :: m, n
  integer, intent(in) :: pitID(m)
  real( kind=8 ), intent(in) :: excess(n)

  real( kind=8 ), intent(out) :: addVol(m)

  integer :: k, c

  addVol = 0.

  ! Update elevation
  do k = 1, m
    c = pitID(k)
    if(c>0)then
      if(excess(c)>0.)then
        addVol(k) = excess(c)
      endif
    endif
  enddo

  return

end subroutine addExcess

subroutine spillPoints(nb,remainSed, pitData, spillPts, n, m)
!*****************************************************************************
! Update spill over point sediment load

  implicit none

  integer :: m, n, nb
  real( kind=8 ), intent(in) :: remainSed(n)
  real( kind=8 ), intent(in) :: pitData(m,5)

  real( kind=8 ), intent(out) :: spillPts(nb,3)

  integer :: k, p, c
  spillPts = 0.

  ! Update elevation
  c = 1
  do k = 1, m
    p = int(pitData(k,3))
    if(p>0)then
      if(remainSed(p)>0.)then
        spillPts(c,1) = pitData(k,1)
        spillPts(c,2) = pitData(k,4)
        spillPts(c,3) = remainSed(p)
        c = c+1
      endif
    endif
  enddo

  return

end subroutine spillPoints

subroutine fillDepression(dem, fillp, pitID, watershed, graph, nv, elev, depressionID, m, nb)
!*****************************************************************************
! Fill depressions from priority-flood calculation

  implicit none

  integer :: m, nb
  integer, intent(in) :: nv
  integer, intent(in) :: watershed(m)
  real( kind=8 ), intent(in) :: dem(m)
  real( kind=8 ), intent(in) :: fillp(m)
  integer, intent(in) :: pitID(m)
  real( kind=8 ), intent(in) :: graph(nb)

  real( kind=8 ), intent(out) :: elev(m)
  integer, intent(out) :: depressionID(m)

  integer :: k, pitnb, nn
  integer :: pitVal(nv)
  real( kind=8 ) :: fillwatershed

  pitnb = 0
  depressionID = -1
  pitVal = -1

  do k = 1, m
    nn = watershed(k)+1
    fillwatershed = graph(nn)
    if(dem(k) < fillp(k) .and. fillp(k) >= fillwatershed)then
      elev(k) = fillp(k)
      if(pitID(k)>0)then
        if(pitVal(pitID(k))<0)then
          pitnb = pitnb+1
          pitVal(pitID(k)) = pitnb
        endif
        depressionID(k) = pitVal(pitID(k))
      endif
    elseif(dem(k) < fillwatershed)then
      elev(k) = fillwatershed
      if(pitVal(nn)<0)then
        pitnb = pitnb+1
        pitVal(nn) = pitnb
      endif
      depressionID(k) = pitVal(nn)
    else
      elev(k) = dem(k)
    endif
  enddo

  return

end subroutine fillDepression

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
        if(gbounds(nn)==0 .and. FVeLgt(n,p)>0)then
          smean = smean+(elev(n)-elev(nn))/FVeLgt(n,p)
          kk = kk+1
        endif
      enddo
      if(kk>0) bslp(n) = smean/kk
    enddo

    return

end subroutine meanSlope

subroutine flatBounds(elev, bID, gbounds, be, nb, b)
!*****************************************************************************
! Define flat boundary conditions

    use meshparams
    implicit none

    integer :: b,nb
    integer, intent(in) :: bID(b)
    integer, intent(in) :: gbounds(nb)
    real( kind=8 ), intent(in) :: elev(nb)

    real( kind=8 ), intent(out) :: be(nb)

    integer :: k, p, n, nn, kk
    real( kind=8 ) :: esum

    be = elev
    do k = 1, b
      n = bID(k)+1
      esum = 0.
      kk = 0
      do p = 1, FVnNb(n)
        nn = FVnID(n,p)+1
        if(gbounds(nn)==0)then
          esum = esum+elev(nn)
          kk = kk + 1
        endif
      enddo
      if(kk>0) be(n) = esum/kk
    enddo

    return

end subroutine flatBounds

subroutine slpBounds(elev, bslp, bID, gbounds, be, nb, b)
!*****************************************************************************
! Define slope boundary conditions

    use meshparams
    implicit none

    integer :: b,nb
    integer, intent(in) :: bID(b)
    integer, intent(in) :: gbounds(nb)
    real( kind=8 ), intent(in) :: elev(nb)
    real( kind=8 ), intent(in) :: bslp(nb)

    real( kind=8 ), intent(out) :: be(nb)

    integer :: k, p, n, nn, kk
    real( kind=8 ) :: v1,v2

    be = elev
    do k = 1, b
      n = bID(k)+1
      kk = 0
      v1 = 0.
      v2 = 0.
      do p = 1, FVnNb(n)
        nn = FVnID(n,p)+1
        if(gbounds(nn)==0.and.FVeLgt(n,p)>0.)then
          v1 = v1+elev(nn)/FVeLgt(n,p)
          v2 = v2+1./FVeLgt(n,p)
          kk = kk + 1
        endif
      enddo
      if(kk>0 .and. v2.ne.0.)then
        be(n) = (kk*bslp(n)+v1)/v2
      endif
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SEDIMENT DIFFUSION FUNCTIONS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initDiffCoeff(nb, dt, Kds, Kdm, sC, sM, mindt)
!*****************************************************************************
! Initialise diffusion coefficients

    use meshparams
    implicit none

    integer :: nb

    real( kind=8 ), intent(in) :: dt
    real( kind=8 ), intent(in) :: Kds
    real( kind=8 ), intent(in) :: Kdm

    real( kind=8 ), intent(out) :: sC(nb,12)
    real( kind=8 ), intent(out) :: sM(nb,12)
    real( kind=8 ), intent(out) :: mindt

    integer :: k, p
    real( kind=8 ) :: c1,c2,Kdmax

    sC = 0.
    sM = 0.

    mindt = dt
    Kdmax = Kds
    if(Kds<Kdm) Kdmax = Kdm

    do k = 1, nb
      if(FVarea(k)>0.)then
        c1 = Kds/FVarea(k)
        c2 = Kdm/FVarea(k)
        do p = 1, FVnNb(k)
          if(FVvDist(k,p)>0. .and. FVeLgt(k,p)>0.)then
            mindt = min(mindt,FVeLgt(k,p)*FVeLgt(k,p)*0.25/Kdmax)
            sC(k,p) = c1*FVvDist(k,p)/FVeLgt(k,p)
            sM(k,p) = c2*FVvDist(k,p)/FVeLgt(k,p)
          endif
        enddo
      endif
    enddo

    mindt = max(1.,mindt)
    sC = sC
    sM= sM

    return

end subroutine initDiffCoeff

subroutine getMaxEro( sealvl, inIDs, elev, elev0, clDi, cmDi, Cero, nb )
!*****************************************************************************
! Compute maximum erosion thickness for diffusion to ensure stability

  use meshparams
  implicit none

  integer :: nb
  real( kind=8 ), intent(in) :: sealvl
  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: cmDi(nb,12)
  real( kind=8 ), intent(in) :: clDi(nb,12)
  real( kind=8 ), intent(in) :: elev(nb)
  real( kind=8 ), intent(in) :: elev0(nb)

  real( kind=8 ), intent(out) :: Cero(nb)

  integer :: k, n, p
  real( kind=8 ) :: kd, val0

  Cero = 1.

  do k = 1, nb
    val0 = 0.
    if(inIDs(k)>0)then
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        kd = clDi(k,p)
        if(n>0 .and. FVeLgt(k,p)>0.)then
          if(elev(k)<sealvl .and. elev(n)<sealvl)then
            kd = cmDi(k,p)
          elseif(elev(k)>=sealvl .and. elev(n)>=sealvl)then
            kd = clDi(k,p)
          elseif(elev(k)<sealvl .and. elev(n)>=sealvl)then
            kd = clDi(k,p)
          elseif(elev(k)>=sealvl .and. elev(n)<sealvl)then
            kd = clDi(k,p)
          else
            if(elev(k)>=sealvl) kd = clDi(k,p)
            if(elev(k)<sealvl) kd = cmDi(k,p)
          endif
          if(elev(k)>elev(n))then
            val0 = val0+kd*(elev(n)-elev(k))
          endif
        endif
      enddo
    endif
    if(val0<0 .and. elev(k)>elev0(k))then
      if(val0<elev0(k)-elev(k))then
        Cero(k) = (elev0(k)-elev(k))/val0
      endif
    elseif(elev(k)==elev0(k))then
      Cero(k) = 0.
    endif
  enddo

end subroutine getMaxEro

subroutine getDiffElev( sealvl, inIDs, elev, Cero, clDi, cmDi, dh, nb )
!*****************************************************************************
! Compute elevation change due to diffusion

  use meshparams
  implicit none

  integer :: nb
  real( kind=8 ), intent(in) :: sealvl
  integer, intent(in) :: inIDs(nb)
  real( kind=8 ), intent(in) :: cmDi(nb,12)
  real( kind=8 ), intent(in) :: clDi(nb,12)
  real( kind=8 ), intent(in) :: elev(nb)
  real( kind=8 ), intent(in) :: Cero(nb)

  real( kind=8 ), intent(out) :: dh(nb)

  integer :: k, n, p
  real( kind=8 ) :: kd, val, val0

  dh = 0.

  do k = 1, nb
    val0 = 0.
    if(inIDs(k)>0)then
      do p = 1, FVnNb(k)
        n = FVnID(k,p)+1
        kd = clDi(k,p)
        if(n>0 .and. FVeLgt(k,p)>0.)then
          if(elev(k)<sealvl .and. elev(n)<sealvl)then
            kd = cmDi(k,p)
          elseif(elev(k)>=sealvl .and. elev(n)>=sealvl)then
            kd = clDi(k,p)
          elseif(elev(k)<sealvl .and. elev(n)>=sealvl)then
            kd = clDi(k,p)
          elseif(elev(k)>=sealvl .and. elev(n)<sealvl)then
            kd = clDi(k,p)
          else
            if(elev(k)>=sealvl) kd = clDi(k,p)
            if(elev(k)<sealvl) kd = cmDi(k,p)
          endif
          val = elev(n)-elev(k)
          if(val>0.)then
            val0 = val0+kd*val*Cero(n)
          elseif(val<0.)then
            val0 = val0+kd*val*Cero(k)
          endif
        endif
      enddo
      dh(k) = val0
    endif
  enddo

  return

end subroutine getDiffElev

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

  integer :: k, n, p,kk
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
  real( kind=8 ), intent(in) :: circumcenter(n,3)

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
          call euclid( midpoint(1:3), circumcenter(cell_ids(k,cid)+1,1:3),  dist)
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

  allocate(FVarea(nLocal))
  allocate(FVnNb(nLocal))
  allocate(FVnID(nLocal,12))
  allocate(FVeLgt(nLocal,12))
  allocate(FVvDist(nLocal,12))

  FVarea = area
  FVnNb = ngbNb
  FVnID = ngbID
  FVeLgt = edgeLgt
  FVvDist = voroDist

end subroutine defineTIN
