# -*- coding: utf-8 -*-

import time
import traceback

from time import sleep
from threading import Thread, Condition, Event, Lock
try:
    from Queue import Queue, Empty, Full
except ImportError:
    from queue import Queue, Empty, Full

from biTK.ngs.concurrent import AbstractLocalThreadPool, ThreadTask, ThreadWorker

__all__ = ['ThreadPool']


class ThreadPool(AbstractLocalThreadPool):
    """
    Pool of threads consuming tasks from a queue
    """
    def __init__(self, num_threads, CheckInterval = 0.5, verbose=False):
        
        super(ThreadPool, self).__init__(num_threads=num_threads,  
                                        WorkerCheckInterval=CheckInterval)
        self.verbose = verbose
        self._next_thread_id = 0
        # call super class's add_threads_if_needed()
        self.add_threads_if_needed()

    def add_worker_thread(self):
        """
        Implemented by derived class and return a thread object
        """
        w = ThreadWorker(self.tasks, 
                         self.dead_thread_queue, 
                         self.shutdown_event,
                         self.WorkerCheckInterval,
                         verbose=self.verbose)

        w.name = "Thread pool #{}".format(self._next_thread_id)
        self._next_thread_id += 1
        return w

    def map(self, func, iterable, chunksize=None, callback=None):
        """
            mapping func 
        """
        if not hasattr(iterable, '__len__'):
            iterable = list(iterable)

        if chunksize is None:
            chunksize, extra = divmod(len(iterable), len(self._pool) * 4)
            if extra:
                chunksize += 1

        if len(iterable) == 0:
            chunksize = 0

        task_batches = self._get_tasks(func, iterable, chunksize)
        #result = MapResult(self._cache, chunksize, len(iterable), callback)
        #for i, x in enumerate(task_batches):
        #    self.add_task("task_{}".format(i), func, *args, **kwargs)

#    def add_task(self, taskid, func, *args, **kwargs):
    def add_task(self, taskid, task, args=(), kwargs={}, callback=None):
        """
        Add a task to the queue
        """
        if not callable(func):
            return False
        else:
            #th_task = ThreadTask(taskid, func, *args, **kwargs)
            th_task = ThreadTask(taskid, func, args=args, kwargs={}, callback=callback)
            self.tasks.put(th_task)
            self.add_threads_if_needed()
        return True

    def shutdown(self):
        """ 
        Shutdown the pool 
        """
        if not self.tasks.empty():
            self.wait_completion()
        self.shutdown_event.set()
        
        # Check the dead threads
        while not self.dead_thread_queue.empty():
            self.dead_thread_queue.get().dismiss()
        
        # Give threads time to die gracefully
        for th in self.threads:
            time.sleep(self.WorkerCheckInterval + 1)
            if isinstance(th, ThreadWorker):
                th.dismiss()
        del self.threads


#    def getNextTask(self):
#        """ 
#        Retrieve the next task from the task queue. For use only by Worker 
#        objects contained in the pool. 
#        """
#        self.__taskLock.acquire()
#        try:
#            if self.tasks.empty():
#                #return (None, None, None, None, None)
#                return (None, None)
#            else:
#                return self.tasks.get(0)
#        finally:
#            self.__taskLock.release()

#    def wait_completion(self, waitForTasks=True, waitForThreads = True):
#        """
#        Wait for completion of all the tasks in the queue
#        """
#        # block task queueing
#        self.__isJoining = True
#        # wait for tasks to finish
#        if waitForTasks:
#            while not self.tasks.empty():
#                #print("Current size of tasks:{}".format(self.tasks.qsize()))
#                sleep(.1)
#        #quit for all threads
#        self.__resizeLock.acquire()
#        try:
#            #Wait until all threads have exited
#            if waitForThreads:
#                for th in self.__threads:
#                    th.dismiss()
#                for th in self.__threads:
#                    th.join()
#                    del th
#            self.__setThreadCountNolock(0)
#            # Reset the pool for potential reuse
#            self.__isJoining = False
#        finally:
#            self.__resizeLock.release()
