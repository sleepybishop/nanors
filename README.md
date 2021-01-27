# nanors

nanors is a tiny, performant implementation of reed solomon codes capable of reaching multi-gigabit speeds on a single core.

## Performance
Good performance is dependent on CPU acceleration, to this end a stripped down and amalgamated version of https://github.com/sleepybishop/oblas is included. Automatic CPU capabilities is attempted via compiler macro defines.

## Benchmark
![](graph.png)

## Use Cases
Applications with small data volumes and low latency requirements, it was designed to be compatible with the reed solomon scheme used in NVIDIA GameStream.
