#ifndef N5_PARALLEL_H
#define N5_PARALLEL_H

/*
 * parallel-for 中一个独立 task 的 callback。
 *
 * task_index 在 0..task_count-1 中唯一出现一次。
 * callback 必须自行保证：
 *   不同 task 写入不同 result slots，或使用显式同步。
 */
typedef void (*N5ParallelTaskFunction)(
    unsigned int task_index,
    void* user_context
);

/*
 * 返回当前机器可用的 logical processor 数。
 *
 * Windows 使用跨 processor groups 的 active processor count；
 * Linux 优先使用当前进程 CPU affinity/cpuset；
 * 其他 POSIX 系统使用 _SC_NPROCESSORS_ONLN。
 *
 * 返回值至少为 1；不使用人为设置的固定线程上限。
 */
unsigned int n5_parallel_hardware_threads(void);

/*
 * 选择 production 默认 worker 数。
 *
 * requested_threads：
 *   0 表示自动使用 hardware logical processor 数；
 *   >0 表示用户显式请求；最终不会超过当前机器可用逻辑处理器数。
 *
 * task_count：
 *   worker 数不会超过 task 数。
 */
unsigned int n5_parallel_resolve_thread_count(
    unsigned int requested_threads,
    unsigned int task_count
);

/*
 * 执行一个确定性的 parallel-for。
 *
 * callback 执行顺序不保证；
 * 但每个 task_index 恰执行一次。
 *
 * 主调用线程也参与执行 tasks，避免额外浪费一个 CPU core。
 *
 * 返回：
 *   1 表示线程创建/回收成功；
 *   0 表示系统线程 API 失败。
 */
int n5_parallel_for(
    unsigned int task_count,
    unsigned int thread_count,
    N5ParallelTaskFunction function,
    void* user_context
);

/*
 * monotonic wall-clock seconds。
 *
 * 用于 production 总耗时与 per-orbit wall time；
 * 不使用 process CPU clock，因此多线程时仍表示真实经过时间。
 */
double n5_parallel_wall_time_seconds(void);

#endif
