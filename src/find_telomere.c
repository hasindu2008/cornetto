/**
 * @file find_telomere.c
 * @author Kavindu Jayasooriya (k.jayasooriya@unsw.edu.au)
 *
rc and find functions are adapted from the `find_telomere.c` in the VGP assembly repository:
https://github.com/VGP/vgp-assembly/blob/master/pipeline/telomere/find_telomere.c

VGP Assembly project is licensed under BSD license
https://github.com/VGP/vgp-assembly/blob/master/LICENSE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <zlib.h>
#include <assert.h>
#include "error.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

//todo lowercase?
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

static void disambiguate(char* str) {
   while (*str) {
      *str = toupper(*str);
      str++;
   }
}

int find_telomere_main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Error: invalid number of parameters\n");
        fprintf(stderr, "Usage: find <input fasta> [optional sequence to search for, default is vertebrate TTAGGG]\n");
        exit(EXIT_FAILURE);
    }

   const char* fasta = argv[1];
   const char* query = (argc >= 3 ? argv[2] : "TTAGGG");

   gzFile fp;
   kseq_t *seq;
   int l;
   fp = gzopen(fasta, "r");
   F_CHK(fp,fasta);
   seq = kseq_init(fp);
   MALLOC_CHK(seq);

   while ((l = kseq_read(seq)) >= 0) {
      assert(l==(int)strlen(seq->seq.s));
      disambiguate(seq->seq.s);
      find(seq->seq.s, seq->name.s, query);
   }

   kseq_destroy(seq);
   gzclose(fp);

   return EXIT_SUCCESS;
}
