# include <stdlib.h>
# include <stdio.h>
# include <string.h>

int main(int argc, char *argv[]){
  
  if(argc!=2){
    printf("USAGE: %s <inputfile>\n", argv[0]);
    return 1;
  }
  char param[128];
  FILE *fp;

  system("clear");  
  fp = fopen(argv[1], "r");
  while(feof(fp)==0){
    if(fgets(param,100,fp),(param[0]=='$'))
      printf("%s", param);
  }
  fclose(fp);
  system("sleep 5s");
  
  return 0;
}
