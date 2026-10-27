/* User-mode core cycles of every OpenMP thread, for the [IGLOO-timing] report (IGLOO_Mod_Timing).
 *
 * One counter per thread, opened by the thread itself from inside a parallel region: an inherited counter only
 * follows tasks created after it is enabled, and the OpenMP pool exists before, so it would count the master
 * alone. perf_event_open on the calling thread is allowed at perf_event_paranoid <= 2 once kernel and hypervisor
 * counting are excluded. A slot that cannot open keeps its first errno for the report; nothing else fails.
 */
#define _GNU_SOURCE
#include <errno.h>
#include <stdint.h>
#include <string.h>

#if defined(__linux__)
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

#define IGLOO_CYC_MAX 1024

/* fd + 1 and errno + 1 per slot: static zero already means "unset", with no init pass to race against. */
static int fds1[IGLOO_CYC_MAX];
static int err1[IGLOO_CYC_MAX];
static int n_open = 0;

/* Opens a cycle counter on the CALLING thread, recorded in `slot` (once per OpenMP thread, inside a region). */
void igloo_cyc_open_thread(const int *slot)
{
    int s = *slot;
    if (s < 0 || s >= IGLOO_CYC_MAX || fds1[s] != 0) return;
#if defined(__linux__)
    {
        struct perf_event_attr a;
        int fd;
        memset(&a, 0, sizeof(a));
        a.type           = PERF_TYPE_HARDWARE;
        a.size           = sizeof(a);
        a.config         = PERF_COUNT_HW_CPU_CYCLES;
        a.disabled       = 1;
        a.exclude_kernel = 1;
        a.exclude_hv     = 1;
        /* pid 0 = the calling thread, cpu -1 = wherever it runs */
        fd = (int)syscall(__NR_perf_event_open, &a, 0, -1, -1, 0);
        if (fd < 0) {
            if (err1[s] == 0) err1[s] = errno + 1;
            return;
        }
        ioctl(fd, PERF_EVENT_IOC_RESET, 0);
        ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
        fds1[s] = fd + 1;
        __sync_fetch_and_add(&n_open, 1);
    }
#else
    if (err1[s] == 0) err1[s] = ENOSYS + 1;
#endif
}

/* How many threads carry a counter (0: cycle counting unavailable). */
void igloo_cyc_count(int *n)
{
    *n = n_open;
}

/* Why `slot` has no counter: 0 it has one (or was never tried), 1 refused (EACCES/EPERM),
   2 a limit (EMFILE/ENFILE/ENOMEM), 3 another errno (no such counter, not Linux, ...). */
void igloo_cyc_slot_error(const int *slot, int *cls)
{
    int s = *slot, e;
    *cls = 0;
    if (s < 0 || s >= IGLOO_CYC_MAX || fds1[s] != 0 || err1[s] == 0) return;
    e = err1[s] - 1;
    if (e == EACCES || e == EPERM)                     *cls = 1;
    else if (e == EMFILE || e == ENFILE || e == ENOMEM) *cls = 2;
    else                                               *cls = 3;
}

/* Cycles since the counters were opened, summed over every thread that carries one. */
void igloo_cyc_read_all(int64_t *total)
{
    int64_t sum = 0;
#if defined(__linux__)
    int i;
    for (i = 0; i < IGLOO_CYC_MAX; i++) {
        uint64_t v = 0;
        if (fds1[i] == 0) continue;
        if (read(fds1[i] - 1, &v, sizeof(v)) == (ssize_t)sizeof(v)) sum += (int64_t)v;
    }
#endif
    *total = sum;
}
