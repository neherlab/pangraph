data := data/synthetic/accuracy
d := script

compare:
	$(julia) $(d)/ancestors.jl -N 100 -L 200000 | pangraph generate -m 0 -t 75 -r 5e-2 -o $(data)/known.json | pangraph build --circular -m 0 -b 0 > $(data)/guess.json
