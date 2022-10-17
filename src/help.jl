Help = Command(
	"help",
	"pangraph help <command>",
	"prints usage information for a given command",
	"one command name",
	Arg[],
	function(args)
		# print version
		run(Version,[])

		command = parse(Help, args)
		command === nothing  && return usage(Dispatch)
		length(command) == 0 && return usage(Dispatch)
		length(command)  > 1 && return usage(Help)

		for sub in Dispatch.sub
			if sub.cmd == command[1]
				usage(sub)
				return 0
			end
		end

		return usage(Help)
	end
)
