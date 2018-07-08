module GAFF {
  // Global Alignment with Scoring Matrix and Affine Gap Penalty
  // http://rosalind.info/problems/gaff/
  // 
  // Useful background info: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf   
  config const defaultfilename = "none";
  config const infile = defaultfilename;
  config const scorefile = "blosum.txt";
  config const gap_open = -10;
  config const gap_extend = -1;

  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
    // Iterate over fasta sequences
    var sequence = "", seqlabel = "";
    var started = false;
    var line = "";

    while (inchannel.readline(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..).strip();
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }
    if started then
      yield (seqlabel, sequence);
  }

  proc read_table(infile) throws {
    // Read substitution score matrix, 
    // return it as a matrix of integers,
    // and a table mapping bases to integers
    try! {
      var file = open(infile, iomode.r);
      var channel = file.reader();

      var line: string;
      channel.readline(line);

      var bases: domain(string);
      var baseorder: [bases] int;
      for (i, base) in zip(1.., line.split()) {
        bases.add(base);
        baseorder[base] = i;
      }

      const n = bases.size;
      var table : [1..n, 1..n] int;

      while (channel.readline(line)) {
        const items = line.split();
        const base = items[1];
        const i = baseorder[base];
        for (j, valstr) in zip(1.., items[2..]) {
          const val = valstr:int;
          table[i, j] = val;
        }
      }
      return (baseorder, table);
    }
  }

  proc alignment(seq1, seq2, gap, extend, bases, blosum) {
    const n = seq1.length;
    const m = seq2.length;

    // Three score matrices:
    //     - M: best score for aligning seq1[1..i-1], seq2[1..j-1]
    //          with sequence ending in a match
    //     - X: best score with sequence ending in a gap in seq1
    //     - Y: best score with sequence ending in a gap in seq2
    // Three dimensions (e.g., M_ijk).  Third dimension:
    //    - i,j,0: score
    //    - i,j,1: i from previous move
    //    - i,j,2: j from previous move
    //    - i,j,3: matrix of previous move (0=m, 1=x, 2=y)
    var scores : [0..3, 0..n, 0..m, 0..3] int;
    const m_mtx = 0;
    const x_mtx = 1;
    const y_mtx = 2;

    // boundary conditions:
    //   - no alignment of either length 0 that ends in a match,
    //     or gap in that sequence
    for i in 1..n {
      scores[m_mtx, i, 0, 0] = -2147483648 / 2;
      scores[x_mtx, i, 0, 0] = gap + i*extend;
      scores[y_mtx, i, 0, 0] = -2147483648 / 2;
    }
    for j in 1..m {
      scores[m_mtx, 0, j, 0] = -2147483648 / 2;
      scores[x_mtx, 0, j, 0] = -2147483648 / 2;
      scores[y_mtx, 0, j, 0] = gap + j*extend;
    }

    for i in 1..n {
      const basei = seq1(i);
      const basei_idx = bases[basei];
      for j in 1..m {
        const basej = seq2(j);
        const basej_idx = bases[basej];

        // match
        var options = [(scores[m_mtx, i-1, j-1, 0], i-1, j-1, m_mtx),
                       (scores[x_mtx, i-1, j-1, 0], i-1, j-1, x_mtx),
                       (scores[y_mtx,i-1, j-1, 0], i-1, j-1, y_mtx)];
        scores[m_mtx, i, j, 0..3] = max reduce options;
        scores[m_mtx, i, j, 0] = scores[m_mtx, i, j, 0] + blosum[basei_idx, basej_idx];

        // match gap in seq 1
        options = [(scores[m_mtx, i, j-1, 0] + gap + extend, i, j-1, m_mtx),
                   (scores[x_mtx, i, j-1, 0] + extend, i, j-1, x_mtx),
                   (scores[y_mtx, i, j-1, 0] + gap + extend, i, j-1, y_mtx)];
        scores[x_mtx, i, j, 0..3] = max reduce options;

        // gap in seq 2
        options = [(scores[m_mtx, i-1, j, 0] + gap + extend, i-1, j, m_mtx),
                   (scores[x_mtx, i-1, j, 0] + gap + extend, i-1, j, x_mtx),
                   (scores[y_mtx, i-1, j, 0] + extend, i-1, j, y_mtx)];
        scores[y_mtx, i, j, 0..3] = max reduce options;
      }
    }

    // backtrack
    var (dist, mtx) = maxloc reduce zip(scores[0..2, n, m, 0], 0..2);
    var i = n;
    var j = m;
    var aug1, aug2: string;

    while i > 0 || j > 0 {
        var move_i = scores[mtx, i, j, 1];
        var move_j = scores[mtx, i, j, 2];
        var move_mtx = scores[mtx, i, j, 3];
        if mtx == m_mtx {
            aug1 = seq1[i] + aug1;
            aug2 = seq2[j] + aug2;
        } else if mtx == x_mtx {
            aug1 = '-' + aug1;
            aug2 = seq2[j] + aug2;
        } else {
            aug1 = seq1[i] + aug1;
            aug2 = '-' + aug2;
        }

        (mtx, i, j) = (move_mtx, move_i, move_j);
    }

    return (dist, aug1, aug2);
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      const (bases, blosum) = read_table(scorefile);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var (n, aug1, aug2) = alignment(seqs[2], seqs[1], gap_open, gap_extend, bases, blosum);
      writeln(n);
      writeln(aug2);
      writeln(aug1);
    }
  }
}
