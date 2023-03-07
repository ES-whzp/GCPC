**GCPC**(GC content and PAMs count Calculation), is a tool for calculating the GC content and PAM count in different strains.

The tool supports six types of PAMs, including 'TTN', 'YTN', 'KYTV', 'KTTV', 'ATTC', and 'NGG'.


**Usage**:

```bash
 python GCPC.py [-h] [-d] [-o]
 -h --help, show this help message and exit.
 -d --genomeDir, The directory that stores the fasta-format genome files.
 -o --output, The output fileName (default name is 'GCPC_output'). The output is saved under the genome directory.
```

**Example**:

python /path/to/GCPC.py -d /path/to/genomeDirectory -o outputName
