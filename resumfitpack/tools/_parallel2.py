#!/usr/bin/env python
import sys,os
import time
import random
import zmq
import multiprocessing
import numpy as np
import cPickle as pickle
import psutil
import mkl
import cProfile
import distutils.util

class PARALLEL2:

    def __init__(self):
        pass


    def setup_master(self):

        self.coordinate = distutils.util.strtobool(os.getenv('PARALLEL_COORDINATE','False'))
        print 'parallel.coordinate=%r'%self.coordinate
        self.cores = [i for i in range(psutil.cpu_count(logical=True))]
        if 'PARALLEL_CORES' in os.environ: self.cores = [int(s) for s in os.getenv('PARALLEL_CORES').split(',')]
        print 'parallel.cores=%r'%','.join(str(i) for i in self.cores)
        self.profile_workers = distutils.util.strtobool(os.getenv('PARALLEL_PROFILE_WORKERS','False'))
        print 'parallel.profile_workers=%r'%self.profile_workers
        self.master_core_lock = int(os.getenv('PARALLEL_MASTER_CORE',-1))

        self.nworkers = len(self.cores)
        if self.master_core_lock >= 0 and self.master_core_lock in self.cores:
            self.nworkers -= 1

        self.context = zmq.Context()
        self.assign_sock = self.context.socket(zmq.REP)
        self.assign_port = self.assign_sock.bind_to_random_port("tcp://*", min_port=49152, max_port=65536, max_tries=100)
        self.result_sock = self.context.socket(zmq.REP)
        self.result_port = self.result_sock.bind_to_random_port("tcp://*", min_port=49152, max_port=65536, max_tries=100)
        if self.coordinate:
            self.acquire_sock = self.context.socket(zmq.REQ)
            self.acquire_sock.connect("tcp://localhost:55550")
            self.release_sock = self.context.socket(zmq.REQ)
            self.release_sock.connect("tcp://localhost:55551")

        self.set_procname('Master')


    def acquire_mutex(self):

        if self.coordinate:
            self.acquire_sock.send("")
            self.acquire_sock.recv()


    def release_mutex(self):

        if self.coordinate:
            self.release_sock.send("")
            self.release_sock.recv()


    def setup_workers(self):

        self.workers = []
        for idx in range(self.nworkers): 
            if self.profile_workers:
                process = multiprocessing.Process(target=self.profile,name='w'+str(idx),args=(idx,))
            else:
                process = multiprocessing.Process(target=self.worker,name='w'+str(idx),args=(idx,))
            process.start()
            self.workers.append(process)
            #psutil.Process(process.pid).cpu_affinity([self.cores[idx]])

        if self.master_core_lock >= 0:
            psutil.Process(os.getpid()).cpu_affinity([self.master_core_lock])
        else:
            psutil.Process(os.getpid()).cpu_affinity(self.cores)
        mkl.set_num_threads(1)


    def profile(self,idx):
        cProfile.runctx('self.worker(idx)',globals(),locals(),'worker-%d.prof' %idx)


    def set_procname(self,name):
        from ctypes import cdll, create_string_buffer, byref
        libc = cdll.LoadLibrary('libc.so.6')
        buff = create_string_buffer(len(name)+1)
        buff.value = name;
        libc.prctl(15, byref(buff), 0, 0, 0)


    def worker(self,idx):
        mkl.set_num_threads(1)
        psutil.Process(os.getpid()).cpu_affinity([self.cores[idx]])

        context = zmq.Context()
        assign_sock = context.socket(zmq.REQ)
        assign_sock.connect("tcp://localhost:%d"%self.assign_port)
        result_sock = context.socket(zmq.REQ)
        result_sock.connect("tcp://localhost:%d"%self.result_port)

        self.set_procname('Worker %d'%idx)
        while True:
            assign_sock.send('')
            idx, state = assign_sock.recv_multipart()
            idx = int(idx)
            if idx < 0: break
            state = pickle.loads(state)
            self.set_state(state)
            result = self.task(idx)
            result_sock.send(pickle.dumps(result,pickle.HIGHEST_PROTOCOL))
            result_sock.recv()

        assign_sock.close()
        result_sock.close()
        context.destroy()


    def update_workers(self,state):

        state_pickle = pickle.dumps(state,pickle.HIGHEST_PROTOCOL)
        self.acquire_mutex()
        for i in range(self.nworkers):
            self.assign_sock.recv()
            self.assign_sock.send_multipart([repr(i),state_pickle])


    def send_tasks(self,requests):

        chunks=[]
        for i in range(self.nworkers):
            msg = self.result_sock.recv()
            self.result_sock.send('')
            results = pickle.loads(msg)
            chunks.append(results)
        self.release_mutex()
        return chunks


    def stop_workers(self):

        for i in range(self.nworkers):
            self.assign_sock.recv()
            self.assign_sock.send_multipart([repr(-1),''])

        self.assign_sock.close()
        self.result_sock.close()
        if self.coordinate:
            self.acquire_sock.close()
            self.acquire_sock.close()
        for w in self.workers:
            w.join()
        self.context.destroy()

def orchestrate():
    context = zmq.Context()
    give = context.socket(zmq.REP)
    give.bind("tcp://*:55550")
    take = context.socket(zmq.REP)
    take.bind("tcp://*:55551")
    while True:
        print 'await give'
        give.recv()
        give.send("")
        print 'await take'
        take.recv()
        take.send("")

if __name__ == "__main__":
    psutil.Process(os.getpid()).cpu_affinity([psutil.cpu_count(logical=True)-1])
    orchestrate()








