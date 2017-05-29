#include <ctype.h>

// reverse the given null-terminated string in place
void inplace_reverse(char * str){
  if(str){
    char * end = str + strlen(str) - 1;

    // swap the values in the two given variables
    // XXX: fails when a and b refer to same memory location
    #define XOR_SWAP(a,b) do\
        {\
        a ^= b;\
        b ^= a;\
        a ^= b;\
        } while(0)

      // walk inwards from both ends of the string,
      // swapping until we get to the middle
      while(str < end){
        XOR_SWAP(*str, *end);
        str++;
        end--;
      }

      #undef XOR_SWAP
    }
}


int get_qsnapshot(double* n, double* wn, double* r, double* wmu, Node* node1, Node* node2){
  int        Index = 0;

  char*              p;
  
  char     header[200];
  long  node_counts[2];
  
  // format:
  // ## Walltime: %.2lf seconds; Progress: %.6lf \%; number of nodes in tree: %ld; current nodes: %ld \t %ld \n", node1->label, node2->label);
  // ## number of pairs; sum over pairs of weights; sum over pairs of weighted r, sum over pairs of weighted mu. \n"); 
  
  sprintf(filepath, "%s/Qmultipoles/%s_raw.dat", outputdir, surveyType);  // Save raw counts.

  inputfile = fopen(filepath, "r");

  fgets(header, 200, inputfile);
  
  printf("\n\nHeader retrieved: %s", header);

  // cut string at current nodes. 
  p = strstr(header, "current nodes:"); // will have the address of ":"
  
  while(*p){
    if(isdigit(*p)){
      long val     = strtol(p, &p, 10);
     
      node_counts[Index] = val;
    }

    else{
      p++;
    }
  }

  for(j=0; j<2; j++)  printf("\n\nNode counts: %ld", node_counts[j]);
  
  // get second header line to move file pointer. 
  fgets(header, 200, inputfile); 

  printf("\n\nInputing Q-multipoles snapshot.");
  
  for(j=0; j<NBINS; j++){
    fscanf(inputfile, "%.4le \t %.4le \t %.4le \t %.4le \t %*.4le \n", &n[j], &wn[j], &r[j], &wmu[j]);

    printf("%.4le \t %.4le \t %.4le \t %.4le \n", n[j], wn[j], r[j], wmu[j]);
  }
  
  fclose(inputfile);
  
  return 0;
}
