import multiprocessing
import multiprocessing.connection
from multiprocessing.connection import Connection
from multiprocessing.context import Process
from typing import Any, Callable, Dict, Iterable, List

"""
Note, AWS Lambda doesn't support shared memory locks because /dev/shm
does not exist. Therefore we use pipes here instead of queues.
"""

SENTINEL = "SENTINEL"


class SingleProcessGeneratorPool:
    def __init__(self, *args, **kwargs):
        self.pending = []

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def __len__(self):
        return 1

    def submit(self, func, kwargss, batchsize=1024):
        self.pending.append((func, kwargss))

    def results(self):
        while self.pending:
            func, kwargss = self.pending.pop()
            for kwargs in kwargss:
                results = func(**kwargs)
                for result in results:
                    yield kwargs, result


class MultiprocessGeneratorPool:
    subprocs: List[Process]
    pipesparent: List[Connection]
    pipeschild: List[Connection]

    def __init__(self, ncpus: int = multiprocessing.cpu_count()):
        self.subprocs = []
        self.pipesparent = []
        self.pipeschild = []
        for i in range(ncpus):
            pipeparent, pipechild = multiprocessing.Pipe(duplex=True)
            subproc = multiprocessing.Process(
                target=self._multiprocess_generator_pool_child,
                args=(pipechild,),
            )

            self.pipesparent.append(pipeparent)
            self.pipeschild.append(pipechild)
            self.subprocs.append(subproc)
        # start them all
        for subproc in self.subprocs:
            subproc.start()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for i in range(len(self.pipesparent)):
            self.pipesparent[i].send(SENTINEL)
        for i in range(len(self.pipesparent)):
            self.pipesparent[i].close()
        for i in range(len(self.subprocs)):
            self.subprocs[i].join(1)

    def __len__(self):
        return len(self.subprocs)

    def submit(
        self, func: Callable, kwargss: Iterable[Dict[Any, Any]], batchsize: int = 1024
    ):
        i = 0
        for kwargs in kwargss:
            # send this set of kwarguments
            self.pipesparent[i].send([func, batchsize, kwargs])
            # point to the next pipe, wrapping if necessary
            i = i + 1 if i + 1 < len(self.subprocs) else 0

    def results(self):
        # wait until we have enough sentinels
        sentinelcount = 0
        while sentinelcount < len(self.subprocs):
            # check which pipes have data in them
            # process each result pipe in turn
            for pipe in multiprocessing.connection.wait(self.pipesparent):
                result = pipe.recv()
                # we recieved a sentinel value to say that a chunk is complete
                if result == SENTINEL:
                    sentinelcount += 1
                elif result[2]:
                    # an exception was raised in a worker
                    # reraise it in the parent
                    # pool as context manager will handle cleanup
                    raise result[2] from result[2]
                else:
                    # note this will be out of order between subprocesses
                    # so we include the fkwargs for disambiguation by the caller if necessary
                    assert len(result) == 3, f"expected 3 got {result}"
                    kwargs, batch, _ = result
                    for item in batch:
                        yield kwargs, item

    @staticmethod
    def _multiprocess_generator_pool_child(pipe: Connection) -> None:
        while True:
            # wait for a message
            pipe.poll(None)
            # read the message
            msg = pipe.recv()
            if msg == SENTINEL:
                # received a sentinel message to terminate self
                break
            else:
                # start of a new function invocation
                func, batchsize, kwargs = msg
                # start a fresh batch of results
                batch = []
                try:
                    for result in func(**kwargs):
                        batch.append(result)
                        # batch is full, send it and start a new one
                        if len(batch) >= batchsize:
                            pipe.send([kwargs, batch, None])
                            batch = []
                except Exception as e:
                    # if an error happened send it up
                    pipe.send([kwargs, batch, e])
                    # continue to the next arg
                else:
                    # no error was thrown
                    # send any leftover lines smaller than a batch
                    pipe.send([kwargs, batch, None])
                # send a sentinel to say we've finished an arg
                pipe.send(SENTINEL)
