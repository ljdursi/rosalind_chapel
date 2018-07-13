module LAFF {
  // Local Alignment with Scoring Matrix and Affine Gap Penalty
  // http://rosalind.info/problems/laff/
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
    //    - i,j,1: matrix of previous move (0=m, 1=x, 2=y)
    var moves : [0..n, 0..m, 0..3, 0..1] int;
    const m_mtx = 0;
    const x_mtx = 1;
    const y_mtx = 2;
    const start = 3;

    // boundary conditions:
    //   - no alignment of either length 0 that ends in a match,
    //     or gap in that sequence
    for i in 1..n {
      moves[i, 0, m_mtx, 0] = -1073741824;
      moves[i, 0, x_mtx, 0] = 0; // gap + i*extend;
      moves[i, 0, y_mtx, 0] = -1073741824;
    }
    for j in 1..m {
      moves[0, j, m_mtx, 0] = -1073741824;
      moves[0, j, x_mtx, 0] = -1073741824;
      moves[0, j, y_mtx, 0] = 0; // gap + j*extend;
    }

    for i in 1..n {
      const basei = seq1(i);
      const basei_idx = bases[basei];

      // var m_mtx_score, xmtx_score, y_mtx_score;
      var scores : [0..3] int;

      for j in 1..m {
        const basej = seq2(j);
        const basej_idx = bases[basej];

        // gap in seq 2
        scores[m_mtx] = moves[i-1, j, m_mtx, 0] + gap + extend;
        scores[x_mtx] = moves[i-1, j, x_mtx, 0] + gap + extend;
        scores[y_mtx] = moves[i-1, j, y_mtx, 0] + extend;
        scores[start] = 0;

        moves[i, j, y_mtx, 0] = scores[m_mtx];
        moves[i, j, y_mtx, 1] = m_mtx;

        for mtx in 1..3 do
          if (scores[mtx] > moves[i, j, y_mtx, 0]) {
            moves[i, j, y_mtx, 0] = scores[mtx];
            moves[i, j, y_mtx, 1] = mtx;
          }

        // match
        scores[m_mtx] = moves[i-1, j-1, m_mtx, 0] + blosum[basei_idx, basej_idx];
        scores[x_mtx] = moves[i-1, j-1, x_mtx, 0] + blosum[basei_idx, basej_idx];
        scores[y_mtx] = moves[i-1, j-1, y_mtx, 0] + blosum[basei_idx, basej_idx];
        scores[start] = 0;

        moves[i, j, m_mtx, 0] = scores[m_mtx];
        moves[i, j, m_mtx, 1] = m_mtx;

        for mtx in 1..3 do
          if (scores[mtx] > moves[i, j, m_mtx, 0]) {
            moves[i, j, m_mtx, 0] = scores[mtx];
            moves[i, j, m_mtx, 1] = mtx;
          }

        // gap in seq 1
        scores[m_mtx] = moves[i, j-1, m_mtx, 0] + gap + extend;
        scores[x_mtx] = moves[i, j-1, x_mtx, 0] + extend;
        scores[y_mtx] = moves[i, j-1, y_mtx, 0] + gap + extend;
        scores[start] = 0;

        moves[i, j, x_mtx, 0] = scores[m_mtx];
        moves[i, j, x_mtx, 1] = m_mtx;

        for mtx in 1..3 do
          if (scores[mtx] > moves[i, j, m_mtx, 0]) {
            moves[i, j, x_mtx, 0] = scores[mtx];
            moves[i, j, x_mtx, 1] = mtx;
          }
      }
    }

    // writeln(moves[1..n, 1..m, m_mtx, 1]);
    // writeln('--');
    // writeln(moves[1..n, 1..m, x_mtx, 1]);
    // writeln('--');
    // writeln(moves[1..n, 1..m, y_mtx, 1]);

    // backtrack
    var i = n;
    var j = m;
    var mtx = m_mtx;
    var dist = moves[i, j, mtx, 0];
    for ii in 1..n do
      for jj in 1..m do
        for kk in 0..2 do
          if moves[ii, jj, kk, 0] > dist {
            dist = moves[ii, jj, kk, 0];
            mtx = kk;
            i = ii;
            j = jj;
          }
    
    var sub1, sub2: string;
    while i > 0 && j > 0 && moves[i, j, mtx, 0] > 0 && moves[i, j, mtx, 1] != start {
        var move_i, move_j : int;
        var move_mtx = moves[i, j, mtx, 1];
        if mtx == m_mtx {
            sub1 = seq1[i] + sub1;
            sub2 = seq2[j] + sub2;
            move_i = i-1;
            move_j = j-1;
        } else if mtx == x_mtx {
            sub2 = seq2[j] + sub2;
            move_i = i;
            move_j = j-1;
        } else {
            sub1 = seq1[i] + sub1;
            move_i = i-1;
            move_j = j;
        }

        (mtx, i, j) = (move_mtx, move_i, move_j);
    }

    return (dist, sub1, sub2);
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      const (bases, blosum) = read_table(scorefile);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var (n, sub1, sub2) = alignment(seqs[2], seqs[1], gap_open, gap_extend, bases, blosum);
      writeln(n);
      writeln(sub2);
      writeln(sub1);
    }
  }
}
