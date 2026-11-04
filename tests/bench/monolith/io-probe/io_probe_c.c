/* C helpers of the I/O probe: glibc's strtod at an offset, fma, a conversion in x87 extended precision, the
   process's CPU times, the calling thread's user-mode core cycles, the CPU a thread runs on, and dropping a file's
   pages from the page cache. */
#define _GNU_SOURCE
#include <fcntl.h>
#include <linux/perf_event.h>
#include <math.h>
#include <sched.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <unistd.h>

/* strtod at s[pos] (0-based); *end is the 0-based offset of the first character it did not take. */
double probe_strtod(const char *s, int64_t pos, int64_t *end)
{
  char *e;
  double x = strtod(s + pos, &e);
  *end = (int64_t)(e - s);
  return x;
}

double probe_fma(double a, double b, double c) { return fma(a, b, c); }

/* m times or divided by 10^|e| in x87 extended precision (exact operands: m < 2^63, |e| <= 27), then rounded to
   double: the double rounding of a runtime that converts through long double. */
static const long double P10L[28] = {1e0L, 1e1L, 1e2L, 1e3L, 1e4L, 1e5L, 1e6L, 1e7L, 1e8L, 1e9L, 1e10L, 1e11L, 1e12L,
  1e13L, 1e14L, 1e15L, 1e16L, 1e17L, 1e18L, 1e19L, 1e20L, 1e21L, 1e22L, 1e23L, 1e24L, 1e25L, 1e26L, 1e27L};

double probe_x87(int64_t m, int e)
{
  volatile long double x = (long double)m;
  if (e >= 0) x = x * P10L[e];
  else x = x / P10L[-e];
  return (double)x;
}

/* User and system CPU time of the process, all threads [s]. The barrier makes the compiler read the structure
   getrusage filled (icx at -O2 otherwise reads it as never written). */
void probe_times(double *user, double *sys)
{
  struct rusage u;
  memset(&u, 0, sizeof u);
  getrusage(RUSAGE_SELF, &u);
  __asm__ __volatile__("" : : "r"(&u) : "memory");
  *user = (double)u.ru_utime.tv_sec + 1e-6 * (double)u.ru_utime.tv_usec;
  *sys = (double)u.ru_stime.tv_sec + 1e-6 * (double)u.ru_stime.tv_usec;
}

/* One user-mode cycle counter for the calling thread; -1 when the host refuses it. */
static __thread int cyc_fd = -1;

int probe_cyc_open(void)
{
  struct perf_event_attr a;
  if (cyc_fd >= 0) return cyc_fd;
  memset(&a, 0, sizeof a);
  a.type = PERF_TYPE_HARDWARE;
  a.size = sizeof a;
  a.config = PERF_COUNT_HW_CPU_CYCLES;
  a.exclude_kernel = 1;
  a.exclude_hv = 1;
  cyc_fd = (int)syscall(__NR_perf_event_open, &a, 0, -1, -1, 0);
  return cyc_fd;
}

int64_t probe_cyc_read(void)
{
  int64_t v = -1;
  if (cyc_fd < 0 || read(cyc_fd, &v, sizeof v) != (ssize_t)sizeof v) return -1;
  return v;
}

int probe_cpu(void) { return sched_getcpu(); }

/* Drop the clean pages of a file from the page cache (no privilege needed); 0 on success. */
int probe_evict(const char *path)
{
  int r, fd = open(path, O_RDONLY);
  if (fd < 0) return -1;
  r = posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);
  close(fd);
  return r;
}
