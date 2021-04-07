# pangraph

> a bioinformatic toolkit to align large sets of genome assemblies into a graph data structure

### example run

verbose logging:
julia src/PanGraph.jl build [your sequence file(s) here] >output.json

quiet logging:
julia src/PanGraph.jl build [your sequence file(s) here] 2>/dev/null 1>output.json

#### to come:
julia src/PanGraph.jl (generate)

### license

[MIT License](LICENSE)
