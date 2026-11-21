#if !defined(_WIN32) && defined(__linux__)
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif

#include <limits.h>
#include <stdlib.h>

#include "n5_parallel.h"

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <process.h>

/*
 * Windows parallel-for 的共享调度状态。
 *
 * next_task 从 -1 开始，InterlockedIncrement 后得到下一个 index。
 */
typedef struct
{
    volatile LONG next_task;
    unsigned int task_count;
    N5ParallelTaskFunction function;
    void* user_context;
} N5ParallelShared;

/*
 * Windows worker body。
 */
static unsigned __stdcall n5_parallel_worker_windows(
    void* raw_context
)
{
    N5ParallelShared* shared;

    shared = (N5ParallelShared*)raw_context;

    for (;;)
    {
        LONG index;

        index = InterlockedIncrement(&shared->next_task);

        if (index < 0 ||
            (unsigned long)index >= (unsigned long)shared->task_count)
        {
            break;
        }

        shared->function(
            (unsigned int)index,
            shared->user_context
        );
    }

    return 0;
}

#else

#include <pthread.h>
#include <time.h>
#include <unistd.h>
#if defined(__linux__)
#include <sched.h>
#endif

/*
 * POSIX parallel-for 的共享调度状态。
 */
typedef struct
{
    unsigned int next_task;
    unsigned int task_count;
    N5ParallelTaskFunction function;
    void* user_context;
    pthread_mutex_t mutex;
} N5ParallelShared;

/*
 * POSIX worker body。
 */
static void* n5_parallel_worker_posix(
    void* raw_context
)
{
    N5ParallelShared* shared;

    shared = (N5ParallelShared*)raw_context;

    for (;;)
    {
        unsigned int index;
        int have_task;

        pthread_mutex_lock(&shared->mutex);

        if (shared->next_task < shared->task_count)
        {
            index = shared->next_task;
            ++shared->next_task;
            have_task = 1;
        }
        else
        {
            index = 0;
            have_task = 0;
        }

        pthread_mutex_unlock(&shared->mutex);

        if (!have_task)
        {
            break;
        }

        shared->function(index, shared->user_context);
    }

    return NULL;
}

#endif

/*
 * 返回当前机器可用的 logical processor 数。
 */
unsigned int n5_parallel_hardware_threads(void)
{
    unsigned int count;

#ifdef _WIN32
    DWORD system_count;

    system_count = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);

    if (system_count == 0)
    {
        SYSTEM_INFO info;

        GetSystemInfo(&info);
        system_count = info.dwNumberOfProcessors;
    }

    count = (unsigned int)system_count;
#else
#if defined(__linux__)
    /*
     * Linux: prefer the process CPU-affinity mask rather than the machine-wide
     * online CPU count. This respects taskset/cpuset/container restrictions.
     */
    {
        cpu_set_t allowed_cpus;
        int allowed_count;

        CPU_ZERO(&allowed_cpus);
        if (sched_getaffinity(0, sizeof(allowed_cpus), &allowed_cpus) == 0)
        {
            allowed_count = CPU_COUNT(&allowed_cpus);
            count = allowed_count > 0 ? (unsigned int)allowed_count : 1u;
        }
        else
        {
            long system_count;
            system_count = sysconf(_SC_NPROCESSORS_ONLN);
            count = system_count > 0 && (unsigned long)system_count <= UINT_MAX
                ? (unsigned int)system_count
                : 1u;
        }
    }
#else
    {
        long system_count;
        system_count = sysconf(_SC_NPROCESSORS_ONLN);
        count = system_count > 0 && (unsigned long)system_count <= UINT_MAX
            ? (unsigned int)system_count
            : 1u;
    }
#endif
#endif

    return count == 0u ? 1u : count;
}

/*
 * 选择 production worker 数。
 */
unsigned int n5_parallel_resolve_thread_count(
    unsigned int requested_threads,
    unsigned int task_count
)
{
    unsigned int count;
    unsigned int hardware_count;

    hardware_count = n5_parallel_hardware_threads();
    count = requested_threads == 0u ? hardware_count : requested_threads;

    if (count == 0u)
    {
        count = 1u;
    }

    if (count > hardware_count)
    {
        count = hardware_count;
    }

    if (task_count > 0 && count > task_count)
    {
        count = task_count;
    }

    return count;
}

/*
 * 执行 deterministic-index parallel-for。
 */
int n5_parallel_for(
    unsigned int task_count,
    unsigned int thread_count,
    N5ParallelTaskFunction function,
    void* user_context
)
{
    unsigned int resolved;

    if (function == NULL)
    {
        return 0;
    }

    if (task_count == 0)
    {
        return 1;
    }

    resolved = n5_parallel_resolve_thread_count(
        thread_count,
        task_count
    );

    if (resolved <= 1)
    {
        unsigned int index;

        for (index = 0; index < task_count; ++index)
        {
            function(index, user_context);
        }

        return 1;
    }

#ifdef _WIN32
    {
        N5ParallelShared shared;
        HANDLE* handles;
        unsigned int worker_count;
        unsigned int worker;
        int success;

        shared.next_task = -1;
        shared.task_count = task_count;
        shared.function = function;
        shared.user_context = user_context;

        worker_count = resolved - 1u;

        handles = (HANDLE*)calloc(worker_count, sizeof(HANDLE));

        if (handles == NULL)
        {
            return 0;
        }

        success = 1;

        for (worker = 0; worker < worker_count; ++worker)
        {
            uintptr_t handle_value;

            handle_value = _beginthreadex(
                NULL,
                0,
                n5_parallel_worker_windows,
                &shared,
                0,
                NULL
            );

            if (handle_value == 0)
            {
                success = 0;
                break;
            }

            handles[worker] = (HANDLE)handle_value;
        }

        /* 主线程也参与处理剩余 tasks。 */
        n5_parallel_worker_windows(&shared);

        {
            unsigned int created;

            created = success ? worker_count : worker;

            for (worker = 0; worker < created; ++worker)
            {
                if (WaitForSingleObject(handles[worker], INFINITE) !=
                    WAIT_OBJECT_0)
                {
                    success = 0;
                }

                CloseHandle(handles[worker]);
            }
        }

        free(handles);
        return success;
    }
#else
    {
        N5ParallelShared shared;
        pthread_t* threads;
        unsigned int worker_count;
        unsigned int worker;
        unsigned int created;
        int success;

        shared.next_task = 0;
        shared.task_count = task_count;
        shared.function = function;
        shared.user_context = user_context;

        if (pthread_mutex_init(&shared.mutex, NULL) != 0)
        {
            return 0;
        }

        worker_count = resolved - 1u;

        threads = (pthread_t*)calloc(worker_count, sizeof(pthread_t));

        if (threads == NULL)
        {
            pthread_mutex_destroy(&shared.mutex);
            return 0;
        }

        success = 1;
        created = 0;

        for (worker = 0; worker < worker_count; ++worker)
        {
            if (pthread_create(
                    &threads[worker],
                    NULL,
                    n5_parallel_worker_posix,
                    &shared
                ) != 0)
            {
                success = 0;
                break;
            }

            ++created;
        }

        /* 主线程也参与。 */
        n5_parallel_worker_posix(&shared);

        for (worker = 0; worker < created; ++worker)
        {
            if (pthread_join(threads[worker], NULL) != 0)
            {
                success = 0;
            }
        }

        free(threads);
        pthread_mutex_destroy(&shared.mutex);

        return success;
    }
#endif
}

/*
 * monotonic wall-clock seconds。
 */
double n5_parallel_wall_time_seconds(void)
{
#ifdef _WIN32
    LARGE_INTEGER frequency;
    LARGE_INTEGER counter;

    if (!QueryPerformanceFrequency(&frequency) ||
        !QueryPerformanceCounter(&counter) ||
        frequency.QuadPart == 0)
    {
        return 0.0;
    }

    return (double)counter.QuadPart / (double)frequency.QuadPart;
#else
    struct timespec now;

    if (clock_gettime(CLOCK_MONOTONIC, &now) != 0)
    {
        return 0.0;
    }

    return (double)now.tv_sec + (double)now.tv_nsec / 1000000000.0;
#endif
}
