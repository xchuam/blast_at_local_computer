"""Task execution helpers for blast_at_local_computer."""
from __future__ import annotations

import concurrent.futures
import subprocess
import random
from multiprocessing import Lock, Process, Manager
from pathlib import Path
from time import sleep
from typing import Any, Dict, List, Optional, Tuple

__version__ = "0.1.0"


class BASE_Tasker:
    def __init__(self, op_path: str, V: str = __version__):
        self.version = V
        self.op_path = Path(op_path)
        self.base_cmd_: List[str] = []  # 默认为空列表

    def one_task_cmd_lst_modifier(
        self, task_pack: Tuple[Any, ...], assigned_resource: int
    ) -> List[str]:
        """
        修改基础命令以适配具体任务

        参数:
            task_pack: 任务参数元组
            assigned_resource: 分配的资源数

        返回:
            修改后的命令列表
        """
        # 在子类中应重写此方法
        return self.base_cmd_.copy() + list(task_pack) + [str(assigned_resource)]

    def Psudo_pool_one_sub_excute(
        self,
        task_queue,
        # each element is a task pack. task_pack[0] is the task string, task_pack[1] is the task size
        resource_pool: Dict[str, int],
        resource_lock_: Lock,
        total_size: int,
        available_resource: int,
    ) -> None:
        """动态任务调度器 - 按任务大小比例分配线程"""
        while not task_queue.empty():
            try:
                # 动态获取任务
                with resource_lock_:
                    if task_queue.empty():
                        return
                    task_pack = task_queue.get_nowait()
            except Exception:
                return

            # 计算任务大小比例
            size_ratio = task_pack[1] / total_size

            # 计算所需线程数 (至少1个)
            thread_should_be_given = max(1, int(size_ratio * available_resource))

            # 等待直到有足够资源
            while True:
                with resource_lock_:
                    if resource_pool["available"] >= thread_should_be_given:
                        if resource_pool["tasks_left"] < available_resource:
                            thread_given = max(
                                thread_should_be_given,
                                resource_pool["available"] // resource_pool["tasks_left"],
                            )
                        else:
                            thread_given = thread_should_be_given
                        resource_pool["available"] -= thread_given
                        break
                # wait and check
                sleep(random.randint(1, 2))

            cmd = self.one_task_cmd_lst_modifier(task_pack[0], thread_given)
            # 执行任务
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                error_msg = (
                    f"Error: {cmd} -->\n Err CODE:\n {e.returncode}\n{e.stderr}"
                )
                print(error_msg)

            # 归还资源
            with resource_lock_:
                resource_pool["available"] += thread_given
                resource_pool["tasks_left"] -= 1

            # 标记任务完成
            task_queue.task_done()

    def Psudo_pool_excute(
        self,
        sorted_mission_lst: List[Tuple[str, int]],
        total_size: int,
        available_resource: int,
    ) -> None:
        with Manager() as manager:
            # 创建任务队列
            task_queue = manager.Queue()
            for task_pack in sorted_mission_lst:
                task_queue.put(task_pack)

            # 创建共享资源池
            resource_pool = manager.dict()
            resource_pool["available"] = available_resource  # 初始可用线程数
            resource_pool["tasks_left"] = len(sorted_mission_lst)
            st_lock = manager.Lock()
            # 创建进程 (进程数=min(总任务数, 进程数))
            jobs: List[Process] = []
            actual_procs = min(len(sorted_mission_lst), available_resource)

            for _ in range(actual_procs):
                p = Process(
                    target=self.Psudo_pool_one_sub_excute,
                    args=(
                        task_queue,
                        resource_pool,
                        st_lock,
                        total_size,
                        available_resource,
                    ),
                )
                jobs.append(p)
                p.start()

            # 等待所有任务完成
            task_queue.join()

            # 终止所有进程
            for p in jobs:
                p.join(timeout=1.0)
                if p.is_alive():
                    p.terminate()

    def Pool_worker(
        self,
        task_pack: Tuple[Any, ...],  # 任务参数包
        thread: int = 1,
    ) -> Optional[str]:
        """
        进程池工作函数

        参数:
            task: 任务参数元组
            thread: 分配的线程数

        返回:
            执行结果消息
        """
        cmd = self.one_task_cmd_lst_modifier(task_pack, thread)

        try:
            result = subprocess.run(
                cmd, check=True, capture_output=True, text=True
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            error_msg = f"Error: {cmd} -->\n Err CODE:\n {e.returncode}\n{e.stderr}"
            print(error_msg)
            return None

    def Pool_excute(
        self, tasks: List[Tuple[Any, ...]], processes: int
    ) -> None:
        """
        并行执行任务（使用标准进程池）

        参数:
            tasks: 任务列表
            processes: 使用的进程数
        """
        # 使用进程池执行任务
        with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as pool:
            futures = [pool.submit(self.Pool_worker, *task) for task in tasks]
            for future in concurrent.futures.as_completed(futures):
                future.result()

    def excute(
        self,
        tasks: List[Tuple[Tuple[Any, ...], int]],
        max_workers: Optional[int] = None,
    ) -> None:
        """使用线程池执行任务，适用于 I/O 密集型任务。"""
        if not tasks:
            return

        if max_workers is None:
            max_workers = len(tasks)

        def _run(task: Tuple[Tuple[Any, ...], int]) -> Optional[str]:
            task_pack, thread = task
            return self.Pool_worker(task_pack, thread)

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_run, task) for task in tasks]
            for future in concurrent.futures.as_completed(futures):
                future.result()


class RsyncTasker(BASE_Tasker):
    """简单的 rsync 任务执行器"""

    def __init__(self, op_path: str):
        super().__init__(op_path)
        self.base_cmd_ = ["rsync", "-avzq"]

    def one_task_cmd_lst_modifier(
        self, task_pack: Tuple[Any, ...], assigned_resource: int
    ) -> List[str]:
        if not task_pack:
            raise ValueError("task_pack must include the rsync source address")

        address = str(task_pack[0])
        if len(task_pack) > 1:
            destination = str(task_pack[1])
        else:
            destination = str(self.op_path)

        return self.base_cmd_.copy() + [address, destination]
