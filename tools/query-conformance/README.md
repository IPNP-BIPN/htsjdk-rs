# QueryInterval and Chunk conformance harness

Dumps `QueryInterval`'s ordering, overlap, abuttal and `optimizeIntervals`, and `Chunk`'s
comparison, overlap, adjacency and `optimizeChunkList`.

`Chunk.overlaps`, `isAdjacentTo` and `optimizeChunkList` are package-private and are reached by
reflection. The alternative is to exercise them through a real index, which would measure the index
rather than the arithmetic, and the arithmetic is what a query's chunk list is made of.

## Run

```sh
python3 ../conformance/run_suite.py --suites query
```
