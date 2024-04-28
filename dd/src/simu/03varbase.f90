
!===================================================================================================
!========================    DEBUT    MODULE  "VARBASE"   ==========================================
!===================================================================================================

include '00copyright.f90'

!> \brief This module does not have a brief description yet
!! \todo  Change some comments into doxygen comments
!! \todo  Translate some comments into English
module VARBASE

!> \ingroup Constantes
use CONSTANTES

implicit none

!#########################################################################
!# Variables                                                             #
!#########################################################################

integer(kind=4),allocatable,save ::BVECLINCP(:,:) ! vecteurs ligne en precision simple pour Graph.c

integer(kind=DPI),allocatable,save ::  &

BVECLIN(:,:)    ,&!*** vecteurs ligne de la simulation
BVECDEP(:,:)    ,&!*** vecteurs deplacement associes
BVECNOR(:,:)    ,&!*** normale au sys de glissement
courbuO(:,:)    ,&!*** type des rotules pour
courbuE(:,:)    ,&!***   decrire la bonne courbure
ASSOC(:,:)      ,&!*** index des vecteurs associes (base 2D)
INVSI(:)        ,&!*** renvoie le vecteur de signe oppose
DEASSOC(:)      ,&!*** operation inverse : connaissant i et j donne l'indice reduit de j
TYSEG(:)        ,&!*** type de segment
SISEG(:)        ,&!*** signe du vecteur ligne du segment
SEGdev(:,:)     ,&!*** devseg(i,j) = l'indice du vecteur associe a I dans le system d'ordre J
SEGdevIJ(:,:)   ,&!*** devseg(i,j) = l'indice du vecteur associe a I dans le system du seg J
CONEC(:,:)      ,&!*** segments connectables
DECONEC(:,:)    ,&!*** facdep(i,deconec(i,conec(i,1)),...
SC(:,:)         ,&!*** Pour changer la config au pt d'ancrage SimplConnec
CoefSC(:,:)     ,&!*** le coefficient de remplacementde SC
AXEJONCi(:,:)   ,&!*** axe de jonction entre deux sys de gliss
AXEJONCj(:,:)   ,&!*** axe de jonction entre deux sys de gliss
SYSOBST(:,:)    ,&!*** quel type d'obtacle ???
SYSCONNEC(:,:)  ,&!*** pour connecter les segments entre eux
EtatSYS(:,:)    ,&!*** relation entre les plans de glissment de I et J
NBROT(:,:,:)    ,&!*** Combien de rotules et la ou les quelles
GDROT(:,:)      ,&!*** vis intermediaire eventuelle
INVCONEC(:,:)   ,&!*** permutation de type de seg conect...
SYSEG(:)        ,&!*** systeme de glissement d'un segment
LOMA(:)         ,&!*** Longueur de discretisation tabulee des segments
JSLOMA(:)       ,&!*** Longueur de discretisation tabulee des jonctions
numero_loi(:)   ,&!*** renvoi le numero de la loi de vitesse auquel le vecteur linge est associe
MODDEPLA(:,:,:) ,&!*** le modulo de deplacement /generalisation des mixt/mixt pour les CFC
Inn(:)          ,&!*** indice non nul du vecteur de la base
SysDev(:,:)       !*** table of indexes of systems of cross slip


real(kind=DP),allocatable,save ::          &

FAC1DEP(:,:)       ,& !*** facteurs intervenant dans dynam
FAC2DEP(:,:,:)     ,& !*** facteurs intervenant dans dynam
VECNORLIN(:,:)     ,& !*** Base de vecteurs ligne normalises
VECNORDEP(:,:)     ,& !***                  deplacement
VECNORNOR(:,:)     ,& !***                  normal
NORMLIN(:)         ,& !*** Normes des vecteurs ligne de la base
NORMDEP(:)         ,& !***                  deplacement
partvis(:)         ,& !*** cosinus entre VL et la direction vis du system
projvis(:)         ,& !*** projection of segments direction on Burgers vector direction (not normalized)
projcoin(:)           !*** projection of segments direction on edge direction (not normalized)

end module VARBASE

!===================================================================================================
!========================    FIN    MODULE     =====================================================
!===================================================================================================
