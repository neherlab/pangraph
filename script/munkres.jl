# ------------------------------------------------------------------------
# munkres assignment algorithm
# adaption of Yi Cao's MATLAB algorithm
# https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm

Mark = (
	Star  = 1,
	Prime = 2,

	NotStar = 510
)

function assign(data)
	star(x)  = (x & Mark.Star)  > 0
	prime(x) = (x & Mark.Prime) > 0

	mark = zeros(UInt8, size(data))

	# step 1: subtract minimum from each row
	cost = data .- minimum(data, dims=2)

	# step 2: find zeros of cost. if no zeros in column or row, star the zero
	hits = cost .== 0
	@label step2 # START LOOP =====
	hit = findfirst(hits)
	hit !== nothing || @goto end2

		mark[I] |= Mark.Star
		hits[hit.I[1],:] = false
		hits[:,hit.I[2]] = false

	@goto step2
	@label end2 # END LOOP =====

	# step 3: cover each column with a starred zero
	@label step3 # START LOOP ======
	cover = (
		col = any(mark .== Mark.Star,dims=1),
		row = zeros(Bool, size(cost,2))
	)
	any(cover.col) || @goto end3
	# step 4: find a noncovered zero and prime it.
	# if no starred zero, go to step 5.
	# otherwise cover row and uncover the column of starred zero.
	# continue until no uncovered zeros left.
	# go to step 6
		@label step4 # START LOOP ======
		hits .= false
		index = (
			row = .!cover.row,
			col = .!cover.col,
		)
		hits[index.row, index.col] = cost[index.row, index.col] .== 0

		@label kernel
		hit = findfirst(hits)
		hit !== nothing || @goto endkernel

			mark[hit] |= Mark.Prime
			stars = star.(mark[hit.I[1],:])

			any(stars) || @goto end4

			cover.row[hit.I[1]] = true
			cover.col[hit.I[2]] = true

			hits[hit.I[1],:] = false
			hits[.!cover.row,stars] = (cost[.!cover.row,stars] .== 0)

		@goto kernel
		@label endkernel

		# step 6: add the minimum uncovered value to every element of each covered row
		index = (
			row = .!cover.row,
			col = .!cover.col,
		)
		m = minimum(cost[index.row, index.col])
		cost[cover.row,cover.col] .+= m
		cost[index.row,index.col] .-= m

		@goto step4
		@label end4 # END LOOP =======

		# step 5
		col = hit.I[2]
		row = star.(mark[:,col])
		mark[hit] .|= Mark.Star

		@label step5
		hit = findfirst(row)
		hit !== nothing || @goto end5

			mark[hit,col] .&= Mark.NotStar
			j, i = prime.(mark[row,:]), row
			row  = star.(mark[:,j])
			mark[i,j] .|= Mark.Star

		@goto step5
		@label end5

	@label end3 # END LOOP ======
end
