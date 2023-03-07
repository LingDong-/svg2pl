#define S2P_SETDIM(w,h) {printf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%f\" height=\"%f\"><path stroke=\"black\" fill=\"none\" d=\"\n",(w),(h));}

#include "../svg2pl.c"

int main(int argc, char** argv){
  if (argc > 1){
    s2p_parse_from_file(argv[1]);
  }else{
    FILE* fd;
    char* buf;
    size_t sz;
    fd = open_memstream(&buf,&sz);
    int c;
    while ((c=fgetc(stdin)) != EOF){
      fputc(c, fd);
    }
    fflush(fd);
    fclose(fd);

    fd = fmemopen(buf,sz,"r");
    s2p_parse(fd);

    fclose(fd);
  }
  printf("\" /></svg>");
  
}