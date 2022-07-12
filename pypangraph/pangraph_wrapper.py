# file containing utilities to run pangraph from python code

import subprocess

def __check_pangraph_command():
    rc = subprocess.call(['which', 'pangraph'])
    assert rc == 0, "pangraph command not present on path."

def pangraph_command(
    input_file,
    output_file,
    circular=False,
    min_block_len=None,
    junction_cost=None,
    diversity_cost=None,
    verbose=True,
):
    """Python wrapper for the pangraph pipeline. It executes pangraph on the
    specified input file, and saves the results on the output file. It requires
    the presence of the `pangraph` binary on the path.


    Args:
    - input_file (str) : path to the fasta (optionally gzipped) file to be taken
        as input by PanGraph
    - output_file (str) : name of the file where results are to be stored. Must
        be in `.json` format
    - circular (bool, defalut False) : whether the genomes considered are circular
    - min_block_len (optional int) : the minimum length of blocks
    - junction_cost (optional int) : nergy cost for introducing junction due to
        alignment merger
    - diversity_cost (optional int) : energy cost for interblock diversity due
        to alignment merger

    Returns:
    - the output of subprocess.run. From it one can extract .stderr and .returncode
    """

    __check_pangraph_command()

    # initial command
    command = [
        "pangraph",
        "build",
    ]

    # optional arguments
    if circular:
        command.append("--circular")
    if min_block_len is not None:
        command += ["--len", f"{int(min_block_len)}"]
    if junction_cost is not None:
        command += ["--mu", f"{int(junction_cost)}"]
    if diversity_cost is not None:
        command += ["--beta", f"{int(diversity_cost)}"]

    # input file
    command.append(f"{input_file}")

    def printv(*args):
        if verbose:
            print(*args)
        else:
            return None

    printv("executing pangraph with command:")
    printv(" ".join(command))
    printv(f"saving output in:")
    printv(output_file)
    # run the command and save the results in the specified output
    with open(output_file, "w") as f_out:
        p = subprocess.run(command, stdout=f_out, stderr=subprocess.PIPE)

    print(p.stderr.decode())

    if p.returncode != 0:
        print("Warning: pangraph run with command")
        print(" ".join(command))
        print("saved in", output_file)
        print("failed with exit code", p.returncode)

    return p
