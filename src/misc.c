/* @file misc.c
**
** miscellaneous definitions and inline functions
** @@
******************************************************************************/

#include <sys/resource.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "error.h"

/*

realtime, cputime, peakrss and mm_parse_num
are adapted from https://github.com/lh3/minimap2/blob/master/misc.c

The MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
              2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

double realtime(void) {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

double cputime(void) {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

int64_t mm_parse_num(const char* str)
{
    double x;
    char* p;
    x = strtod(str, &p);
    if (*p == 'G' || *p == 'g')
        x *= 1e9;
    else if (*p == 'M' || *p == 'm')
        x *= 1e6;
    else if (*p == 'K' || *p == 'k')
        x *= 1e3;
    return (int64_t)(x + .499);
}

//parse yes or no arguments
void yes_or_no(uint64_t* flag_a, uint64_t flag, const char* opt_name, const char* arg, int yes_to_set)
{
    if (yes_to_set) {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            *flag_a |= flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            *flag_a &= ~flag;
        } else {
            fprintf(stderr, "option '--%s' only accepts 'yes' or 'no'.\n", opt_name);
        }
    } else {
        if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) {
            *flag_a &= ~flag;
        } else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) {
            *flag_a |= flag;
        } else {
            fprintf(stderr, "option '--%s' only accepts 'yes' or 'no'.\n", opt_name);
        }
    }
}

/* alpha-numeric sort is adapted from https://github.com/samtools/samtools/blob/c83652d0d97544370b5f50abef9141df7a113d03/bam_sort.c

The MIT License

    Copyright (C) 2008-2025 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>
    Author: Martin Pollard <mp15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */



#define is_digit(c) ((c)<='9' && (c)>='0')
int strnum_cmp(const char *_a, const char *_b){

    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (!is_digit(*pa) || !is_digit(*pb)) {
            if (*pa != *pb)
                return (int)*pa - (int)*pb;
            ++pa; ++pb;
        } else {
            // skip leading zeros
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;

            // skip matching digits
            while (is_digit(*pa) && *pa == *pb)
                pa++, pb++;

            // Now mismatching, so see which ends the number sooner
            int diff = (int)*pa - (int)*pb;
            while (is_digit(*pa) && is_digit(*pb))
                pa++, pb++;

            if (is_digit(*pa))
                return  1; // pa still going, so larger
            else if (is_digit(*pb))
                return -1; // pb still going, so larger
            else if (diff)
                return diff; // same length, so earlier diff
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}



// Prints to the provided buffer a nice number of bytes (KB, MB, GB, etc)
//adapted from https://www.mbeckler.org/blog/?p=114
void print_size(const char* name, uint64_t bytes)
{
    const char* suffixes[7];
    suffixes[0] = "B";
    suffixes[1] = "KB";
    suffixes[2] = "MB";
    suffixes[3] = "GB";
    suffixes[4] = "TB";
    suffixes[5] = "PB";
    suffixes[6] = "EB";
    uint64_t s = 0; // which suffix to use
    double count = bytes;
    while (count >= 1024 && s < 7){
        s++;
        count /= 1024;
    }
    if (count - floor(count) == 0.0)
        fprintf(stderr, "[%s] %s : %d %s\n", __func__ , name, (int)count, suffixes[s]);
    else
        fprintf(stderr, "[%s] %s : %.1f %s\n", __func__, name, count, suffixes[s]);
}
