try:
    # this supports e.g. partial funcs, lambdas
    import multiprocessing_on_dill as multiprocessing
except ImportError:
    import multiprocessing


"""
Note, AWS Lambda doesn't support shared memory locks because /dev/shm
does not exist. Therefore we use pipes here instead of queues.
"""

SENTINEL = "SENTINEL"


def _multiprocess_generator_child(f, fkwargs, q, batchsize):
    try:
        batch = []
        for item in f(**fkwargs):
            batch.append(item)
            # batch is full, send it and start a new one
            if len(batch) >= batchsize:
                print("spam", fkwargs)
                q.send([fkwargs, batch])
                batch = []
        # send any leftover lines smaller than a batch
        q.send([fkwargs, batch])
    finally:
        # send a sentinel
        q.send(SENTINEL)
        # close the child-side end of the pipe
        q.close()


def from_multiprocess_generator(
    genfunc,
    genfunckwargs=[{}],
    batchsize=8,  # number of results in a batch
    queuesize=2,  # number of batches queued per subprocess
):
    """
    A generator that uses multiprocessing under the hood
    to distribute generation of results, and gethers then back via a queue.

    genfunc generator function to be invoked to create generators in the subprocesses

    genfunckwargs iterator of keyword arguments to be called in the subprocess generators
    Each item in genfunckwargs is a dict of arguments to genfunc in a subprocess, and
    all subprocesses will be started at the same time e.g. genfunc(**genfunckwargs[0]).

    batchsize is the number of generator results collected to be returned through the queue
    this improves performance by reducing time spent locking and unlocking the queue

    queuesize is the number of batches * number of subprocess that are in the queue, so total
    number of items in the queue is queuesize * batchsize * len(genfunckwargs)
    """

    # multiprocessing.Queue isn't support on AWS Lambda
    # so instead of a shared queue, need one pipe per process
    piperesults = []
    # create the subprocesses
    # TODO wrap in conext manager for cleanup
    subprocs = []
    for i in range(len(genfunckwargs)):
        pipeparent, pipechild = multiprocessing.Pipe(
            duplex=True
        )  # need duplex to flow from child to parent
        subproc = multiprocessing.Process(
            target=_multiprocess_generator_child,
            args=(
                genfunc,
                genfunckwargs[i],
                pipechild,
                batchsize,
            ),
        )
        subprocs.append(subproc)
        piperesults.append(pipeparent)

    # start all the subprocesses
    try:
        for subproc in subprocs:
            subproc.start()

        # as long as we are doing something
        while len(subprocs):
            # check which pipes have data in them
            # process each result pipe in turn
            for piperesult in multiprocessing.connection.wait(piperesults):
                result = piperesult.recv()

                # we recieved a sentinel value to say that a chunk is complete
                if result == SENTINEL:
                    # this pipe can be closed now
                    piperesult.close()
                    i = piperesults.index(piperesult)
                    piperesults.remove(piperesult)
                    # also close the subproc
                    subproc = subprocs[i]
                    subproc.join(1)
                    subprocs.remove(subproc)
                else:
                    # note this will be out of order between subprocesses
                    # so we include the fkwargs for disambiguation by the caller if necessary
                    assert len(result) == 2, f"expected 2 got {len(result)}"
                    fkwargs, batch = result
                    for item in batch:
                        yield fkwargs, item

        assert not len(piperesults), "unclosed pipes remaining"
        assert not len(subprocs), "unterminated subprocesses remaining"
    finally:
        # ensure subprocs are cleaned up by joining
        for subproc in subprocs:
            subproc.join(1)
