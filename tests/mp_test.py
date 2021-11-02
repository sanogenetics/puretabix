import pytest

from puretabix.mp import MultiprocessGeneratorPool


def myrange(stop):
    for i in range(stop):
        yield i
    if stop >= 5:
        # raise an exception
        1 / 0


class TestMultiprocessing:
    def test_success(self):
        with MultiprocessGeneratorPool(2) as pool:
            pool.submit(myrange, [{"stop": 2}, {"stop": 2}])
            results = tuple(sorted((i for _, i in pool.results())))
            assert results == (
                0,
                0,
                1,
                1,
            )
            pool.submit(myrange, [{"stop": 3}, {"stop": 4}])
            results = tuple(sorted((i for _, i in pool.results())))
            assert results == (
                0,
                0,
                1,
                1,
                2,
                2,
                3,
            )

    def test_error(self):
        with pytest.raises(Exception):
            with MultiprocessGeneratorPool(2) as pool:
                pool.submit(myrange, [{"stop": 2}, {"stop": 5}])
                results = tuple(sorted((i for _, i in pool.results())))
                assert results == (
                    0,
                    0,
                    1,
                    1,
                    2,
                    3,
                    4,
                )
