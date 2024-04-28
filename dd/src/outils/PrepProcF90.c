

/* Filte un fichier *.F ou *.INC
	 et écrit le résultat dans fichier *.f ou *.inc */

#include <stdio.h>
#include <string.h>

main(int argc, char **argv, char **envp)
{
   FILE *fichier,*fichier2;
   char *code,*nom,*nom2,*bufnom;
   char pclef[11],mclef[11],bufclef[11];
   char cas;
   int i,j,x,y,mode;
   
argv++;
code=*argv;
argv++;
nom=*argv;  
argv++;
nom2=*argv;  
fichier = fopen(nom,"rt");

if (fichier == NULL)
{
printf("\n\n\n ATTENTION !!! %s pas trouve...\n\n\n",nom);
exit(1);
}
else
{

/*bufnom=nom;
for(;*nom!='.';nom++);
nom++;
for(;*nom!='\0';nom++)*nom=tolower(*nom);
nom=bufnom;*/
fichier2 = fopen(nom2,"w");

   strcpy(pclef,"!+++ ");
   strcat(pclef,code);
   strcat(pclef," +++");
   strcpy(mclef,"!--- ");
   strcat(mclef,code);
   strcat(mclef," ---");
   strcpy(bufclef,"01234567890");
   mode = 2;

   /*
      mode = 0 on verifie qu'il n'y a pas de ! en premiere colonne
           = 1 on verifie qu'il y a des ! en premiere colonne
           = 2 on recopie a l'identique
   */
   
   do
    {
    if(strncmp(bufclef,pclef,11) == 0) switch(mode)
     {
     case 2 : {mode = 0;break;}
     case 0 : {mode = 2;break;}
     }
    
    if(strncmp(bufclef,mclef,11) == 0) switch(mode)
     {
     case 2 : {mode = 1;break;}
     case 1 : {mode = 2;break;}
     }
    x = 0;
    do
     { 
     cas = getc(fichier);
     if(x<11){bufclef[x] = cas;x++;}; 
     if(!feof(fichier)) switch(mode) 
      {
      case 0 : {
	       if(cas !='!' || x>1)
		{
		fputc(cas,fichier2);
		}
	       else
		{
		do
		{
		    cas = getc(fichier);
		    bufclef[x] = cas;
		    x++;    
		}
		while((x<11  && cas != '\n') && !feof(fichier));
		if(x == 11 && strncmp(bufclef,pclef,11) == 0){y=0;}else{y=1;}
		do
		{
		    fputc(bufclef[y],fichier2);
		    y++;    
		}
		while(y<x);
	        }
		break;
	       }
      case 1 : {
	       if(cas !='!' && x==1)
                {fputc('!',fichier2);fputc(cas,fichier2);}
	       else
	        {fputc(cas,fichier2);} 
               break;
	       }
      case 2 : {fputc(cas,fichier2);break;}
      }
     }
    while(cas != '\n' && !feof(fichier));
    }
   while (!feof(fichier));
   fclose(fichier);
   fclose(fichier2);
   printf("%s %s \n",code,nom);
}	 
}
