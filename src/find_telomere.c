/**
 * @file find_telomere.c
 * @author Kavindu Jayasooriya (k.jayasooriya@unsw.edu.au)
 *
Adapted from the `find_telomere.c` in the VGP assembly repository:
https://github.com/VGP/vgp-assembly/blob/master/pipeline/telomere/find_telomere.c

VGP Assembly project is licensed under BSD license
https://github.com/VGP/vgp-assembly/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "error.h"

char* rc(const char* str) {
   size_t len = strlen(str);
   char* DNAseq = (char*)malloc(len + 1);
   if (!DNAseq) return NULL;

   // reverse and complement
   for (size_t i = 0; i < len; ++i) {
      switch (str[len - 1 - i]) {
         case 'A': DNAseq[i] = 'T'; break;
         case 'C': DNAseq[i] = 'G'; break;
         case 'G': DNAseq[i] = 'C'; break;
         case 'T': DNAseq[i] = 'A'; break;
         default: DNAseq[i] = str[len - 1 - i]; break; // handle unexpected characters
      }
   }
   DNAseq[len] = '\0';
   return DNAseq;

}

void find(const char* query, const char* name, const char* target) {
   size_t pos = 0;
   size_t len_query = strlen(query);
   size_t len_target = strlen(target);

   while ((pos = strstr(query + pos, target) - query) < len_query) {
      size_t len = 0;
      printf("%s\t%zu\t0\t%zu", name, len_query, pos);
      while (strncmp(query + pos, target, len_target) == 0) {
         len += len_target;
         pos += len_target;
      }
      printf("\t%zu\t%zu\n", pos, len);
      pos++;
   }

   // reverse search
   char* rev = rc(target);
   pos = 0;
   while ((pos = strstr(query + pos, rev) - query) < len_query) {
      size_t len = 0;
      printf("%s\t%zu\t1\t%zu", name, len_query, pos);
      while (strncmp(query + pos, rev, len_target) == 0) {
         len += len_target;
         pos += len_target;
      }
      printf("\t%zu\t%zu\n", pos, len);
      pos++;
   }
   free(rev);
}

int find_telomere_main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Error: invalid number of parameters\n");
        fprintf(stderr, "Usage: find <input fasta> [optional sequence to search for, default is vertebrate TTAGGG]\n");
        exit(EXIT_FAILURE);
    }

    FILE* infile = fopen(argv[1], "r");
    F_CHK(infile,argv[1]);

    char* str = NULL;
    size_t str_size = 0;
    char line[2048];
    char header[2048] = "";

    while (fgets(line, sizeof(line), infile)) {
        if (strchr(line, '>') == NULL) {
            size_t line_len = strlen(line);
            if (line[line_len - 1] == '\n') {
            line[--line_len] = '\0'; // remove newline character
            }
            str = (char*)realloc(str, str_size + line_len + 1);
            MALLOC_CHK(str);
            strcpy(str + str_size, line);
            str_size += line_len;
        } else {
            if (str_size > 0) {
                find(str, header, (argc >= 3 ? argv[2] : "TTAGGG"));
            }
            free(str);
            str = NULL;
            str_size = 0;
            strncpy(header, line, sizeof(header) - 1);
            header[sizeof(header) - 1] = '\0';
            header[strcspn(header, "\n")] = '\0';
        }
    }
    if (str_size > 0) {
        find(str, header, (argc >= 3 ? argv[2] : "TTAGGG"));
        free(str);
    }
    fclose(infile);

    return EXIT_SUCCESS;
}
