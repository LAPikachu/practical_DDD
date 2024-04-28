MODULE VARGLOB

IMPLICIT NONE

integer,parameter   :: DPI = selected_int_kind(9)
integer,parameter   :: DP = selected_real_kind(p=14)
integer,parameter   :: nsegmax = 100000
integer, parameter  :: TabDim=1000000                     ! Max number of segments accounted for in paraview

integer(kind=DPI),dimension(2,nsegmax) :: tabvois
integer(kind=DPI),dimension(5,TabDim) :: iseg
integer(kind=DPI),dimension(TabDim) :: junctt
integer(kind=DPI)   :: slipsys
real(kind=DP),dimension (3,nsegmax)       :: Vcourb

END MODULE VARGLOB