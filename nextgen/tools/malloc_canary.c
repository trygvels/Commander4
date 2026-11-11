/* Malloc canary (E043): finds writes past the end of large heap blocks made by any library in a process.

   Every allocation of at least THRESH bytes gets PAD extra bytes filled with 0xA5; when the block is freed
   (and at exit) the pad is scanned and any change is reported on stderr as a line starting with OVERFLOW.
   Unlike libgmalloc it does not stop at the first byte, so it measures how far the overflow goes.

   Build and use (macOS):
       cc -O1 -dynamiclib -o libcanary.dylib nextgen/tools/malloc_canary.c
       DYLD_INSERT_LIBRARIES=$PWD/libcanary.dylib .venv/bin/python -I script.py ... 2> canary.err
       grep OVERFLOW canary.err

   Found with it on macOS 27.0 (26A5388g), Apple M4: Accelerate dpotrf with uplo = "U" writes up to 768 kB past
   a workspace for matrix orders of about 1500 to 6000; uplo = "L" does not. */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <pthread.h>

#define DYLD_INTERPOSE(_replacement, _replacee) \
    __attribute__((used)) static struct { const void *replacement; const void *replacee; } _interpose_##_replacee \
    __attribute__((section("__DATA,__interpose"))) = {(const void *)(unsigned long)&_replacement, (const void *)(unsigned long)&_replacee};

#define THRESH ((size_t)64 * 1024)
#define PAD ((size_t)8 * 1024 * 1024)
#define NTRACK 8192
static struct { void *p; size_t size; } tab[NTRACK];
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;

static void track(void *p, size_t size) {
    memset((char *)p + size, 0xA5, PAD);
    pthread_mutex_lock(&mu);
    for (int i = 0; i < NTRACK; i++) if (!tab[i].p) { tab[i].p = p; tab[i].size = size; break; }
    pthread_mutex_unlock(&mu);
}
static int untrack(void *p, size_t *size) {
    int found = 0;
    pthread_mutex_lock(&mu);
    for (int i = 0; i < NTRACK; i++) if (tab[i].p == p) { *size = tab[i].size; tab[i].p = NULL; found = 1; break; }
    pthread_mutex_unlock(&mu);
    return found;
}
static void check(void *p, size_t size, const char *when) {
    unsigned char *t = (unsigned char *)p + size;
    size_t first = PAD, last = 0, cnt = 0;
    for (size_t i = 0; i < PAD; i++) if (t[i] != 0xA5) { if (first == PAD) first = i; last = i; cnt++; }
    if (cnt) {
        char buf[400];
        double v[3] = {0, 0, 0};
        memcpy(v, t + (first & ~(size_t)7), sizeof v);
        int n = snprintf(buf, sizeof buf, "OVERFLOW (%s) block size=%zu bytes: pad bytes changed=%zu, first offset=%zu, last offset=%zu; "
                         "first doubles there: %g %g %g\n", when, size, cnt, first, last, v[0], v[1], v[2]);
        write(2, buf, n);
    }
}
static void *my_malloc(size_t size) {
    if (size >= THRESH) { void *p = malloc(size + PAD); if (p) track(p, size); return p; }
    return malloc(size);
}
static void *my_valloc(size_t size) {
    if (size >= THRESH) { void *p = valloc(size + PAD); if (p) track(p, size); return p; }
    return valloc(size);
}
static void *my_aligned_alloc(size_t al, size_t size) {
    if (size >= THRESH) { void *p = NULL; if (posix_memalign(&p, al < sizeof(void *) ? sizeof(void *) : al, size + PAD)) return NULL; track(p, size); return p; }
    return aligned_alloc(al, size);
}
static int my_posix_memalign(void **out, size_t al, size_t size) {
    if (size >= THRESH) { int r = posix_memalign(out, al, size + PAD); if (!r) track(*out, size); return r; }
    return posix_memalign(out, al, size);
}
static void *my_realloc(void *p, size_t size) {
    size_t old;
    if (p && untrack(p, &old)) check(p, old, "realloc");
    if (size >= THRESH) { void *q = realloc(p, size + PAD); if (q) track(q, size); return q; }
    return realloc(p, size);
}
static void my_free(void *p) {
    size_t old;
    if (p && untrack(p, &old)) check(p, old, "free");
    free(p);
}
__attribute__((destructor)) static void at_exit(void) {
    for (int i = 0; i < NTRACK; i++) if (tab[i].p) check(tab[i].p, tab[i].size, "exit");
}
DYLD_INTERPOSE(my_malloc, malloc)
DYLD_INTERPOSE(my_valloc, valloc)
DYLD_INTERPOSE(my_aligned_alloc, aligned_alloc)
DYLD_INTERPOSE(my_posix_memalign, posix_memalign)
DYLD_INTERPOSE(my_realloc, realloc)
DYLD_INTERPOSE(my_free, free)
