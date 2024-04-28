/*
!==================================================================================================
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
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>

#define STRING "microMegas"
#define BORDER 1
#define FONT "fixed"
#define mc_cursor_width 17      /* curseur special pour la fenetre */
#define mc_cursor_height 17     /* dessine pour la circonstance */
#define mc_cursor_x_hot 0       /* position du point actif */
#define mc_cursor_y_hot 0

static char mc_cursor_bits[] = {
0x03, 0x00, 0x00, 0x0f, 0x00, 0x00, 0x3e, 0x00, 0x00, 0xfe, 0x00, 0x00,
0xfc, 0x03, 0x00, 0xfc, 0x0f, 0x00, 0xf8, 0x00, 0x00, 0xf8, 0x00, 0x00,
0xb0, 0x01, 0x00, 0x30, 0x03, 0x00, 0x20, 0x06, 0x00, 0x20, 0x0c, 0x00,
0x00, 0x18, 0x00, 0x00, 0x30, 0x00, 0x00, 0x60, 0x00, 0x00, 0xc0, 0x00,
0x00, 0x80, 0x01};

Display   *dpy;  /* numero de display */
Window    win;    /* numero de fenetre */
GC        gc;
Pixmap    source,mask;
XEvent    event;
XSizeHints  xsh;
XColor    m;

unsigned int	taille_fenetre = 400; // taille de la fenetre
unsigned long	fg,bg;		/* teintes foreground et background :ici blanc et noir */


// The line setting
unsigned int line_width = 2;    /* 0 would be fast line of width 1 */
int line_style = LineSolid;     /* If LineOnOffDash or LineDoubleDash,
* must set dashes */
int cap_style = CapRound;       /* else CapNotLast, CapButt, or
* CapProjecting */
int join_style = JoinRound;     /* else JoinMiter or JoinBevel */


XWMHints xwmh = {		/* definition des Hints de la fenetre */
(InputHint|StateHint),
False,
ZoomState,0,0,0,0,0,0,};


/**********************************************************************/

void window_init(int argc, char **argv)

/***********************************************************************
*******
******* Ce programme realise une initialisation de fenetre simplifiee
******* permettant de tracer des traits colores a l'aide d'une palette
******* qui y est definie. Le parametre d'entree est taille_fenetre qui
******* est explicite et ouvre une fenetre carree. Le curseur de la
******* fenetre a ete specialement redessine.
*******
***********************************************************************/

{

unsigned long 	bd,bw;      /* definition des parametres entourage de fenetre */
XGCValues   gcv;            /* parametres de contexte graphique */
XSetWindowAttributes xswa;  /* attributs de la fenetre */
Cursor    the_cursor;       /* pointeur de definition d'un curseur */
/*Visual	*visual;*/

getcool(&modecool,&fondcool,&boitecool,trcool,trR,trG,trB,tron,&cooljonc,&cooldev,&nbube);

if ((dpy=XOpenDisplay(NULL)) == NULL){ /* peut-on ouvrir la session ? */
fprintf(stderr,"%s:can't open %s\n",argv[0],XDisplayName(NULL));
exit(1); }

bd = WhitePixel(dpy,DefaultScreen(dpy));/* recherche des couleurs par
                                          defaut */
bg = fondcool;// BlackPixel(dpy,DefaultScreen(dpy));
fg = boitecool;// WhitePixel(dpy,DefaultScreen(dpy));
bw = fondcool; //1

/* annonce a l'ecran sur la facon d'operer */
      fprintf(stderr,"\n\n");
      fprintf(stderr,"Click inside the graphical xterm to access commands mode\n");
      fprintf(stderr,"In commands mode type --> a <-- to see the list of commands \n\n");

xsh.flags = (PPosition );
xsh.height = taille_fenetre;
xsh.width  = taille_fenetre;
xsh.x = 0 ;
xsh.y = 0 ;
win = XCreateSimpleWindow(dpy,DefaultRootWindow(dpy),
                  xsh.x,xsh.y,xsh.width,xsh.height,
                  bw,bd,bg);                      /* creation d'une
                                                   fenetre virtuelle */
XSetStandardProperties(dpy,win,STRING,STRING,None,argv,argc,&xsh);
XSetWMHints(dpy,win,&xwmh);

/*visual = XDefaultVisual(dpy,DefaultScreen(dpy));
xswa.colormap = XCreateColormap(dpy,win,visual,AllocAll);*/
xswa.colormap = DefaultColormap(dpy,DefaultScreen(dpy));/* chargement
                                          d'une palette de couleurs */

xswa.bit_gravity = CenterGravity; /* gravite pour les rafraichissements
                                    d'ecran */
source = XCreateBitmapFromData(dpy,win,mc_cursor_bits,mc_cursor_width,mc_cursor_height);
mask = XCreateBitmapFromData(dpy,win,mc_cursor_bits,mc_cursor_width,mc_cursor_height);
                              /* creation d'un cursuer special */
gcv.foreground = fg;
gcv.background = bg;
gc = XCreateGC(dpy,win,(GCForeground | GCBackground ),&gcv); /* creation d'un contexte graphique */
//gcv.line_width = 2; /* RM 2014 Manifestement il ne trace pas les segment tout petit même pas sous forme de point avec cette option */
//gc = XCreateGC(dpy,win,(GCForeground | GCBackground | GCLineWidth),&gcv);
XSetLineAttributes(dpy, gc, line_width, line_style, cap_style, join_style);


XParseColor(dpy,xswa.colormap,"magenta",&m);
XAllocColor(dpy,xswa.colormap,&m); /* POUR LE CURSEUR */

/*m.pixel = trcool[1];
m.red = trR[1]*65535;
m.green = trG[1]*65535;
m.blue = trB[1]*65535;
m.flags = DoRed&&DoGreen&&DoBlue;
XStoreColor(dpy,xswa.colormap,&m); PB ON A PAS LES DROIT POUR MODIFIER LES PALETTES PEUT ETRE CAR TRUE COLOR*/

the_cursor = XCreatePixmapCursor(dpy,source,mask,&m,&m,mc_cursor_x_hot,mc_cursor_y_hot);
xswa.cursor = the_cursor; 	/* mise en place du curseur */
/* definition des attributs de la fenetre: mise en place des masques */
XChangeWindowAttributes(dpy,win,(CWColormap | CWBitGravity | CWCursor),&xswa);
XSelectInput(dpy,win,ExposureMask | ButtonPressMask | KeyPressMask | StructureNotifyMask );
               /* selection des evenements a prendre en compte */

XMapWindow(dpy,win); 		/* creation d'une fenetre reelle */
}
/**********************************************************************/
