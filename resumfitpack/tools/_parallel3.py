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


def set_procname(self,name):
    from ctypes import cdll, create_string_buffer, byref
    libc = cdll.LoadLibrary('libc.so.6')
    buff = create_string_buffer(len(name)+1)
    buff.value = name;
    libc.prctl(15, byref(buff), 0, 0, 0)


class ZmqBroker(object):

  def __init__(self,partition):
    if not partition:
        partition = 8
        # if 'ncpu' in conf and conf['ncpu']!=None:
        #   partition = conf['ncpu']
        # else:
        #   partition = psutil.cpu_count(logical=True)
    self.partition = partition
    self.process = None
      

  def run(self):
    context = zmq.Context(2)
  
    # Socket facing clients
    frontend = context.socket(ZMQ_ROUTER) 
    frontend.bind("tcp://*:55550")
  
    # Socket facing services
    backend  = context.socket(ZMQ_ROUTER) 
    backend.bind("tcp://*:55551")
    backend.router_mandatory = 1
    #backend.setsockopt(ZMQ_ROUTER_MANDATORY, 1)
    
    # Idle workers
    ready_workers = deque()
    
    # FIFO Queue of unassigned work
    pending_work = deque()
    
    # Accumulation of results from workers, per-client
    results = {}
    
    poller = zmq.Poller()
    poller.register(backend,ZMQ_POLLIN)
    poller.register(frontend,ZMQ_POLLIN)
    
    # Used for debugging
    serial = 0

    workers = {}
    set_procname('Broker')

    while True:
      sockets = dict(poller.poll())
      
      # Incoming from worker
      if backend in sockets:
        request = backend.recv_multipart()
        worker, _, client = request[:3]
        ready_workers.append(worker)
        
        # If this is a result message as opposed to a ready message, it will
        # contain more than three parts
        if len(request) > 3:
          _, sn, result = request[3:]
          results[client].append(result)
          #logging.debug('backend %s serial %d results[%s] %d',workers[worker],int(sn),client,len(results[client]))
          
          # If all workers have generated results, we can pass them to the client
          if len(results[client]) == self.partition:
            msg = [client,'']
            msg.extend(results[client])
            frontend.send_multipart(msg)
            
        else:
            #logging.debug('worker %s joins (%s)',client,worker)
            workers[worker] = client
      
      # Incoming from master
      if frontend in sockets:
        client, _, request = frontend.recv_multipart()
        #logging.debug('frontend %s',client)
        results[client] = []
        pending_work.append([client,request,0])
          
      # Can we assign some work?
      while pending_work and ready_workers:
        worker = ready_workers.popleft()
        client,request,offset = pending_work[0]
        serial += 1
        try:
            backend.send_multipart([worker,'',client,'',repr(serial),repr(offset),request],copy=True)
        except zmq.ZMQError as e:
            #logging.warn('Failed send to %s: %s',workers[worker],e)
            continue
        #logging.debug('send %s serial %d offset %d',workers[worker],serial,offset)
        offset += 1
        if offset < self.partition:
          pending_work[0][2] = offset
        else:
          pending_work.popleft()


  def run_subprocess(self):
      def run_runner():
        #logfile = '%s/broker.log'%conf['outputdir']
        #logging.info('Redirecting broker log to %s',logfile)
        #for handler in logging.root.handlers:
        #  logging.root.removeHandler(handler)
        #level=logging.DEBUG
        #level=logging.NOTSET
        #logging.basicConfig(filename=logfile,filemode='w',format='%(asctime)s - %(levelname)s - %(message)s',level=level)
        #console = logging.StreamHandler()
        #console.setLevel(logging.WARN)
        #formatter = logging.Formatter('%(asctime)s - broker - %(levelname)s - %(message)s')
        #console.setFormatter(formatter)
        #logging.root.addHandler(console)
        #os.dup2(log.fileno(),1)
        #os.dup2(log.fileno(),2)
        self.run()
      self.process = Process(target=run_runner,name='broker')
      self.process.start()
      

  def stop(self):
      if self.process:
          self.process.terminate()
          self.process.join(3)



class PARALLEL3:

    def __init__(self):
        self.cores = [i for i in range(psutil.cpu_count(logical=True))]
        if 'PARALLEL_CORES' in os.environ: self.cores = [int(s) for s in os.getenv('PARALLEL_CORES').split(',')]
        print 'parallel.cores=%r'%','.join(str(i) for i in self.cores)
        self.profile_workers = distutils.util.strtobool(os.getenv('PARALLEL_PROFILE_WORKERS','False'))
        print 'parallel.profile_workers=%r'%self.profile_workers
        self.partition = int(os.getenv('PARALLEL_PARTITION',8));
        self.nworkers = len(self.cores)-1


    def run_workforce(self):
        broker = ZmqBroker()
        self.workers = []
        for idx in range(self.nworkers): 
            if self.profile_workers:
                process = multiprocessing.Process(target=self.profile,name='w'+str(idx),args=(idx,))
            else:
                process = multiprocessing.Process(target=self.worker,name='w'+str(idx),args=(idx,))
            process.start()
            self.workers.append(process)
        broker.run()
        for w in self.workers:
            w.terminate()
            w.join()


    def setup_master(self):
        self.context = zmq.Context()
        self.worksock = self.context.socket(ZMQ_REQ)
        self.worksock.connect('tcp://localhost:55550')
        set_procname('Master')
        psutil.Process(os.getpid()).cpu_affinity([self.cores[-1]])
        mkl.set_num_threads(1)


    def setup_workers(self):
        pass


    def profile(self,idx):
        cProfile.runctx('self.worker(idx)',globals(),locals(),'worker-%d.prof' %idx)


    def worker(self,idx):
        mkl.set_num_threads(1)
        psutil.Process(os.getpid()).cpu_affinity([self.cores[idx]])

        context = zmq.Context()
        worksock = context.socket(ZMQ_REQ)
        worksock.connect('tcp://localhost:55551')

        name = 'Worker %d'%idx
        worksock.send(name)
        set_procname('Worker %d'%idx)

        while True:
            address,_,serial,offset,blob = worksock.recv_multipart()
            serial = int(serial)
            #logging.debug('%s recv %d', name,serial)
            if offset == None: break
            offset = int(offset)
            state = pickle.loads(blob)
            self.set_state(state)
            result = self.task(offset)
            result = pickle.dumps(result,pickle.HIGHEST_PROTOCOL)
            #logging.debug('%s send %d', name,serial)
            worksock.send_multipart([address,'',repr(serial),result])
            #logging.debug('%s sent %d', name,serial)
      
        #logging.debug('%s done',name)
        worksock.close()
        context.destroy()


    def update_workers(self,state):
        state_pickle = pickle.dumps(state,pickle.HIGHEST_PROTOCOL)
        self.worksock.send(state_pickle)


    def send_tasks(self,requests):
        data = self.worksock.recv_multipart()
        return [ pickle.loads(part) for part in data[:] ]


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
