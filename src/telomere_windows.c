// This code is adapted from the VGP Assembly repository
// Original source: https://github.com/VGP/vgp-assembly/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>

#define WINDOW_SIZE 1000
#define MIN_OFFSET 0

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

static double THRESHOLD = 0.4;

void print_usage_telomere_windows() {
    fprintf(stderr, "Usage: cornetto telomere --windows <in> <identity> <threshold>\n");
    fprintf(stderr, "This program sizes a fasta or fastq file. Multiple fasta files can be supplied by using a comma-separated list.\n");
    fprintf(stderr, "Example usage: cornetto telomere --windows fasta1.fasta,fasta2.fasta\n");
}

void process_scaffold(const char* name, uint8_t* b, int length) {
    if (b == NULL) { return; }

    for (int i = MIN_OFFSET; i <= length; i += WINDOW_SIZE / 5) {
        int car = 0;
        for (int j = i; j < i + WINDOW_SIZE && j < length; ++j) {
            if (b[j]) car++;
        }
        int den = (i + WINDOW_SIZE < length) ? WINDOW_SIZE : length - i;
        if ((double)car / den >= THRESHOLD)
            printf("Window\t%s\t%d\t%d\t%d\t%.3g\n", name, length, i, i + den, (double)car / den);

        if (i + WINDOW_SIZE >= length)
            break;
    }
}

int telomere_windows_main(int argc, char* argv[]) {
    if (argc < 3) { print_usage_telomere_windows(); return EXIT_FAILURE; }

    if (argc == 4) {
        THRESHOLD = atof(argv[3]);
    }

    char* input_file = argv[1];
    double identity = atof(argv[2]) / 100;
    THRESHOLD = THRESHOLD * pow(identity, 6); // NOTE:assumes the pattern is 6 bases long
    fprintf(stderr, "Given error rate of %.6f running with adjusted threshold of %.6f due to survival prob %.6f\n", identity, THRESHOLD, pow(identity, 6));

    FILE* fp = fopen(input_file, "r");
    if (!fp) {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    char line[LINE_MAX];
    uint8_t* scaffold = NULL;
    char name[LINE_MAX] = {0};
    int length = 0;

    while (fgets(line, sizeof(line), fp)) {
        char split[6][LINE_MAX];
        sscanf(line, "%s %s %s %s %s %s", split[0], split[1], split[2], split[3], split[4], split[5]);//TODO: buffer overflow

        if (scaffold == NULL || strcmp(split[0], name) != 0) {
            process_scaffold(name, scaffold, length);
            length = atoi(split[1]);
            scaffold = (uint8_t*)calloc(length, sizeof(uint8_t));
            strcpy(name, split[0]);
        }
        int start = atoi(split[3]);
        int end = atoi(split[4]);
        for (int i = start; i < end; ++i) {
            scaffold[i] = 1;
        }
    }
    process_scaffold(name, scaffold, length);
    free(scaffold);
    fclose(fp);

    return EXIT_SUCCESS;
}
