// PrepIncF90

#include <stdio.h>
#include <string.h>

main(int argc, char **argv, char **envp)
{
   char *nom;
   
argv++;
nom=*argv; 

printf("<%s]...\n",nom);

preinc(nom);
}
 
preinc(char *nom)
{

   FILE *fichier,*fichier2,*fichier3;
   char *bufnom;
   char *nom2,*bufnom2,buf[280];   
   char pclef[8],claclef[255],bufclef[255],bufx[255];
   char cas,cas2,mode[2]=" ";
   int i,j,x,y,z;


fichier = fopen(nom,"rt");
if (fichier == NULL)
{
printf("\n\n\n ATTENTION !!! %s pas trouve...\n\n\n",nom);
exit(1);
}
else
{
bufnom=nom;
for(;*nom!='.';nom++);
nom++;
for(;*nom!='\0';nom++)*nom=toupper(*nom);
nom=bufnom;
fichier2 = fopen(nom,"w");
strcpy(claclef  ,"0123456");
strcpy(pclef  ,"INCLUDE");
strcpy(bufclef,"0123456");
      
   do
    {
    if(strncmp(claclef,pclef,7) == 0)
     {
     nom2=&bufclef[7];
     for(;*nom2!=34 && *nom2!=39;nom2++);   //On relance sur le fils
     nom2++;
     bufnom2=nom2;
     for(;*nom2!=34 && *nom2!=39;nom2++);
     *nom2='\0'; 
printf("<(%c)%s]",mode[0],bufnom2);

if (strcmp(nom,bufnom2) == 0)
{
printf("\n\n\n ATTENTION !!! la procedure est recursive...\n\n\n");
exit(1);
}
else
{
     preinc(bufnom2);
}		 
     fichier3 = fopen(bufnom2,"rt");  
     strcpy(buf,"!<(\0");strcat(buf,mode);strcat(buf,")\0");strcat(buf,bufnom2);strcat(buf,"]\n"); 
     for(y=0;buf[y]!='\n';y++)fputc(buf[y],fichier2);
     fputc('\n',fichier2);
     if(mode[0]==33)fputc(mode[0],fichier2);
     do{
      cas2 = getc(fichier3);
      if(!feof(fichier3)){
       fputc(cas2,fichier2);
       if(cas2=='\n' && mode[0]=='!')fputc(mode[0],fichier2);
      }
     }while(!feof(fichier3));
     fclose(fichier3);
     fputc('\n',fichier2);
     strcpy(buf,"![\0");strcat(buf,bufnom2);strcat(buf,">\n"); 
     for(y=0;buf[y]!='\n';y++)fputc(buf[y],fichier2);
     fputc('\n',fichier2);
     }
    x = 0;
    z = 0;
    mode[0] = 32;
    do
     { 
     cas = getc(fichier);
     if(!feof(fichier)) 
      {
       if(x<255){
       
        bufx[x]=cas;
	x++;
	bufclef[z] = cas;claclef[z] = toupper(cas);
	if(cas!=33 && cas!=32 || z!=0)z++;
	if(cas==33 && z<=7) mode[0]=33;
	
	};
	
       if(z>=7){
        if(strncmp(claclef,pclef,7) != 0){
	 if(z==7){
	  for(y=0;y<(x-z);y++)fputc(bufx[y],fichier2);
	  for(y=0;y<z;y++)fputc(bufclef[y],fichier2);
	 }
	 if(z>7)fputc(cas,fichier2);
	}
       };
      }
     if(cas=='\n' && z<7)
      {
      for(y=0;y<x;y++)fputc(bufx[y],fichier2);
      } 
     }
    while(cas != '\n' && !feof(fichier));
    }
   while (!feof(fichier));
   fclose(fichier);
   fclose(fichier2);
printf("[%s>\n",nom);   
}
}
