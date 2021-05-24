try:
    # this supports e.g. partial funcs, lambdas
    import multiprocessing_on_dill as multiprocessing
except ImportError:
    import multiprocessing

SENTINEL = "SENTINEL"


def _multiprocess_generator_child(f, fkwargs, q, batchsize):
    batch = []
    for item in f(**fkwargs):
        batch.append(item)
        # batch is full, send it and start a new one
        if len(batch) >= batchsize:
            q.put((fkwargs, batch))
            batch = []
    # send any leftover lines smaller than a batch
    q.put((fkwargs, batch))
    # send a sentinel
    q.put(SENTINEL)


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

    # this is the queue the results will be returned from subprocs via
    q = multiprocessing.Queue(len(genfunckwargs) * queuesize)

    # create the subprocesses
    # TODO wrap in conext manager for cleanup
    subprocs = []
    for i in range(len(genfunckwargs)):
        subproc = multiprocessing.Process(
            target=_multiprocess_generator_child,
            args=(genfunc, genfunckwargs[i], q, batchsize,),
        )
        subprocs.append(subproc)

    # start all the subprocesses
    for subproc in subprocs:
        subproc.start()

    # to keep track of the number of active subprocs
    active = len(subprocs)
    while active > 0:
        result = q.get()
        # we recieved a sentinel value to say that this chunk is complete
        if result == SENTINEL:
            active -= 1
        else:
            # note this will be out of order between subprocesses
            # so we include the fkwargs for disambiguation by the caller if necessary
            assert len(result) == 2, f"expected 2 got {len(result)}"
            fkwargs, batch = result
            for item in batch:
                yield fkwargs, item

    # make sure all the subprocesses terminate tidily
    # at this point we have recieved all the sentinels, so we should be fine
    for subproc in subprocs:
        subproc.join(1)
