/*
!==================================================================================================
! mM (for microMegas) is an open source program of DD (Dislocation Dynamics) simulation
! originally developed at the 'Laboratoire d'Etude des Microstructure' CNRS-ONERA, FRANCE
! Copyright (C) 1993, 1996, 2000, 2002, 2004   Benoit Devincre, Ladislas Kubin, Marc Condat,
! Christophe Lemarchand, Ronan Madec, Sebastien Groh, Ghiath Monnet, Julien Durinck,
! Phillipe Carrez, Christophe de Sansal, Mohamed-Gazi Tagorti.
!
! This file is part of mM.
!
! mM is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! any later version.
!
! mM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with MM; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
!
! For more information see http://zig.onera.fr/mm_home_page/
!===================================================================================================
*/

#include <stdio.h>

/**********************************************************************/

void aide()

/***********************************************************************
*******
*******	Ce programme est une aide en ligne pour un rappel des touches
******* qui permettent une action a partir de la fenetre de
******* visualisation.
*******
***********************************************************************/

{
fprintf (stderr,"\n   KEY   | EFFECT \
               \n  --------------- \
               \n    a    | Help \
               \n    b    | Bulk:crystal state restoring \
               \n    c    | Color on/off \
               \n    d    | Switch off/on the crystal cube \
               \n    e    | Bulk:crystal state restoring and visualising \
               \n    f    | Scale factor \
               \n    g    | Visualising delay (1-1000)\
               \n    h    | Horizontal rotation (in degrees)\
               \n    i    | Switch every/one segment(s)\
               \n    j    | Thin foil with Diophantien tilming \
               \n    k    | Film and .epsf captures\
               \n    l    | Thin specimen \
               \n    m    | Miller's indices entry (three following integers)\
               \n    N    | Print or Not Seg Number \
               \n    o    | Segments : listing \
               \n    p    | View point (three following numbers)\
               \n    q    | Quit the simulation      \
               \n    r    | Automatic mode (revolver) \
               \n    s    | One step \
               \n    t    | Stepping mode \
               \n    u    | Zoom for forest density when DIOPHANTIEN TILING is used\
               \n    v    | Vertical rotation (in degrees) \
               \n    x    | Info on the segment clicked by the mouse \
               \n    y    | Link with the simulation code -> Datas & Stats \
               \n    z    | Zoom factor \n\n\n");
}
/**********************************************************************/
