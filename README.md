# Data Availability Sampling Efficiency Scripts
This repository contains Python scripts to compute efficiency metrics of data availability sampling schemes. The code given here has been used in [this](https://eprint.iacr.org/2023/1079.pdf) and [this](https://eprint.iacr.org/2024/248.pdf) paper.

Note: The parameters that this script selects internally have not been audited. The script has not been audited. Therefore, parameters picked by this script should not be used in production without further checking.

## Generating Tables
First make sure that the `tabulate` module is installed.
Then, you can generate tables by running
```
uv run table.py <data size in MB>
```
For example, if you are interested in efficiency metrics for a data size of 32 Megabytes, run
```
uv run table.py 32
```
If you want to print the table as LaTeX, add the option `-l`.

## Generating Plots
Run
```
uv run graphs.py
```
As a result, you will find csv files in `./csvdata/`
The script will generate csv files for all schemes and the parameter ranges specified in `graphs.py`.
For each scheme, `./csvdata/` will contain separate csv files for the commitment size, communication per query, total communication, and encoding size.
You can then plot this data, e.g., using LaTeX. 
Here is an example of how to plot the encoding size:
```tex
    \begin{tikzpicture}
      \begin{axis}[
          at={(ax1.south east)},
          xshift=2cm,
          width=0.45\linewidth,
          height=0.3\linewidth,
          grid=major,
          grid style={dashed,gray!30},
          xlabel=\small Data,
          y unit= GB,
          x unit= MB,
          ylabel=\small Encoding,
        ]
        \addplot[red, mark =x]
        table[col sep=comma] {/csvdata/tensor_encoding.csv};
        \addplot[black!70!green, mark =triangle]
        table[col sep=comma] {/csvdata/hash_encoding.csv};
        \addplot[blue,thick]
        table[col sep=comma] {/csvdata/fri_encoding.csv};
      \end{axis}
    \end{tikzpicture}
```

## Codes
Erasure codes (or rather their parameters) are modelled as a dataclass, see `codes.py`.
A code maps a message to codeword. In this script, a code is therefore specified by the following parameters:
- `size_msg_symbol` the size of a message symbol in bits.
- `msg_len` the number of symbols the message contains. 
- `size_code_symbol` the size of a codeword symbol in bits.
- `codeword_len` the number of symbols the codeword contains. 
- `reception` if that many symbols from the codeword are known, then the message can be reconstructed
- `samples` the number of random samples to reconstruct the message with high probability

The last parameter deserves some explanation: assume we randomly sample positions of the codeword and collect the corresponding codeword symbols. Then `samples` models how many random samples guarantee that we collect enough symbols to reconstruct the data. 

When you want to add a new code, you can implement a function that defines and returns the code.
For example, the following snippet returns the trivial code that maps a message to itself:
```python
def makeTrivialCode(symbolsize, msg_len):
    '''
    Identity Code, mapping a message of msg_len many symbos to itself.
    Each symbol is of size symbolsize bits
    '''
    return Code(
        size_msg_symbol=symbolsize,
        msg_len=msg_len,
        size_code_symbol=symbolsize,
        codeword_len=msg_len,
        reception=msg_len,
        samples=samples_from_reception(SECPAR_SOUND, msg_len, msg_len)
    )
```

## Schemes
A data availability sampling scheme (see `schemes.py`) is specified by the following:
- `code` the erasure code that is used.
- `com_size` the size of the commitment in bits.
- `opening_overhead` the size of an opening proof, i.e., the overhead of opening a symbol in the encoding

For example, we can consider a (bad) data availability sampling scheme in which we do not encode the data, and commit to the data using a Merkle tree. Then this can be written as
```python
def makeMerkleScheme(datasize, chunksize=1024):
    '''
    Merkle scheme:
    Take a merkle tree and the identity code
    Every leaf contains chunksize bits of the data
    '''
    k = math.ceil(datasize / chunksize)
    return Scheme(
        code=makeTrivialCode(chunksize, k),
        com_size=HASH_SIZE,
        opening_overhead=math.ceil(math.log(k, 2))*HASH_SIZE
    )
```

## License
MIT License.
