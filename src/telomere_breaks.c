#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "error.h"
#include "misc.h"
#include "khash.h"

#define MIN_TEL 24
#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

typedef struct {
    char* name;
    int length;
    unsigned char* bitset;
} Scaffold;

KHASH_MAP_INIT_STR(scaffold, Scaffold*)
KHASH_MAP_INIT_STR(final_scaffold, Scaffold*)
KHASH_MAP_INIT_STR(length, int)

Scaffold* create_scaffold(const char* name, int length) {
    Scaffold* scaffold = (Scaffold*)malloc(sizeof(Scaffold));
    scaffold->name = strdup(name);
    scaffold->length = length;
    scaffold->bitset = (unsigned char*)calloc((int)ceil(length / 8.0), sizeof(unsigned char));
    return scaffold;
}

void set_bit(unsigned char* bitset, int index) {
    bitset[index / 8] |= (1 << (index % 8));
}

int get_bit(unsigned char* bitset, int index) {
    return bitset[index / 8] & (1 << (index % 8));
}

void free_scaffold(Scaffold* scaffold) {
    free(scaffold->name);
    free(scaffold->bitset);
    free(scaffold);
}

int telomere_breaks_main(int argc, char* argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: telobreaks <lens_file> <sdust_file> <telomere_file>\n");
        return EXIT_FAILURE;
    }

    // Initialize scaffolds
    khash_t(scaffold) *scaffold_map = kh_init(scaffold);
    khash_t(final_scaffold) *final_scaffold_map = kh_init(final_scaffold);
    khash_t(length) *length_map = kh_init(length);

    FILE* lens_file = fopen(argv[1], "r");
    F_CHK(lens_file, argv[1]);

    char line[LINE_MAX];
    while (fgets(line, sizeof(line), lens_file)) {
        char name[LINE_MAX];
        int length;
        sscanf(line, "%s %d", name, &length);
        int ret;
        khiter_t k = kh_put(scaffold, scaffold_map, strdup(name), &ret);
        kh_value(scaffold_map, k) = create_scaffold(name, length);
        k = kh_put(final_scaffold, final_scaffold_map, strdup(name), &ret);
        kh_value(final_scaffold_map, k) = create_scaffold(name, length);
        k = kh_put(length, length_map, strdup(name), &ret);
        kh_value(length_map, k) = length;
    }
    fclose(lens_file);

    FILE* sdust_file = fopen(argv[2], "r");
    F_CHK(sdust_file, argv[2]);

    while (fgets(line, sizeof(line), sdust_file)) {
        char name[LINE_MAX];
        int start, end;
        sscanf(line, "%s %d %d", name, &start, &end);
        khiter_t k = kh_get(scaffold, scaffold_map, name);
        if (k != kh_end(scaffold_map)) {
            Scaffold* scaffold = kh_value(scaffold_map, k);
            for (int j = start; j < end; ++j) {
                set_bit(scaffold->bitset, j);
            }
        }
    }
    fclose(sdust_file);

    FILE* telomere_file = fopen(argv[3], "r");
    F_CHK(telomere_file, argv[3]);

    while (fgets(line, sizeof(line), telomere_file)) {
        char name[LINE_MAX];
        int start, end, matched_len;
        sscanf(line, "%s %*d %*d %d %d %d", name, &start, &end, &matched_len);
        if (matched_len >= MIN_TEL) {
            khiter_t k = kh_get(scaffold, scaffold_map, name);
            if (k != kh_end(scaffold_map)) {
                Scaffold* scaffold = kh_value(scaffold_map, k);
                int rStart = start - 100 < 0 ? 0 : start - 100;
                int rEnd = end + 100 > scaffold->length ? scaffold->length : end + 100;
                int all_set = 1;
                for (int j = rStart; j < rEnd; ++j) {
                    if (!get_bit(scaffold->bitset, j)) {
                        all_set = 0;
                        break;
                    }
                }
                if (all_set) {
                    int length = kh_value(length_map, kh_get(length, length_map, name));//TODO: check for return value != kh_end
                    Scaffold* final_scaffold = kh_value(final_scaffold_map, kh_get(final_scaffold, final_scaffold_map, name));
                    int rStart = start;
                    while (rStart > 0 && get_bit(scaffold->bitset, rStart - 1)) {
                        rStart--;
                    }
                    int rEnd = end;
                    while (rEnd < length && get_bit(scaffold->bitset, rEnd)) {
                        rEnd++;
                    }
                    for (int j = rStart; j < rEnd; ++j) {
                        set_bit(final_scaffold->bitset, j);
                    }
                }
            }
        }
    }
    fclose(telomere_file);

    for (khiter_t k = kh_begin(final_scaffold_map); k != kh_end(final_scaffold_map); ++k) {
        if (kh_exist(final_scaffold_map, k)) {
            Scaffold* scaffold = kh_value(final_scaffold_map, k);
            for (int i = 0; i < scaffold->length; ++i) {
                if (get_bit(scaffold->bitset, i)) {
                    int end = i;
                    while (end < scaffold->length && get_bit(scaffold->bitset, end)) {
                        ++end;
                    }
                    i=i-1 < 0 ? 0 : i-1;
                    printf("Found telomere positions %d to %d is a telomere in %s of length %d\n", i, end - 1, scaffold->name, scaffold->length);
                    i = end;
                }
            }
        }
    }

    for (khiter_t k = kh_begin(scaffold_map); k != kh_end(scaffold_map); ++k) {
        if (kh_exist(scaffold_map, k)) {
            free_scaffold(kh_value(scaffold_map, k));
        }
    }
    kh_destroy(scaffold, scaffold_map);

    for (khiter_t k = kh_begin(final_scaffold_map); k != kh_end(final_scaffold_map); ++k) {
        if (kh_exist(final_scaffold_map, k)) {
            free_scaffold(kh_value(final_scaffold_map, k));
        }
    }
    kh_destroy(final_scaffold, final_scaffold_map);

    for (khiter_t k = kh_begin(length_map); k != kh_end(length_map); ++k) {
        if (kh_exist(length_map, k)) {
            free((char*)kh_key(length_map, k));
        }
    }
    kh_destroy(length, length_map);

    return EXIT_SUCCESS;
}
