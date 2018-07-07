module OSYM {
  // Isolating Symbols in Alignments
  // http://rosalind.info/problems/osym/
  const defaultfilename="none";
  config const infile=defaultfilename;
  enum move {diag=0, ins=1, del=2, best=3};

  proc input_channel(filename: string, defaultfilename: string) throws {
    // return open file from filename filename if provided,
    // or stdin
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
    // Iterate over sequnces in a FASTA
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

  proc matchscore(seq1, i, seq2, j) {
    if (seq1[i] == seq2[j]) then
      return +1;
    return -1;
  }

  proc rev(seq) {
    var rseq = "";
    for base in seq[1..seq.length by -1] do
      rseq += base;
    return rseq;
  }

  proc glob_align_scores(seq1, seq2) {
    // Create scoring matrix from global alignment of two sequences
    // 3rd dimension of matrix = operation
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m, move.diag..move.best] int;

    for k in move.diag..move.best {
      score[1..n, 0, k] = -n..-1 by -1;
      score[0, 1..m, k] = -m..-1 by -1;
    }

    for i in 1..n {
      for j in 1..m {
        const bp_score = matchscore(seq1, i, seq2, j);
        score[i, j, move.diag] = score[i-1, j-1, move.best] + bp_score;
        score[i, j, move.ins] = score[i-1, j, move.best] - 1;
        score[i, j, move.del] = score[i, j-1, move.best] - 1;
        score[i, j, move.best] = max(score[i, j, move.diag], score[i, j, move.ins], score[i, j, move.del]);
      }
    }
    return score;
  }

  proc isolation_matrix(seq1, seq2) {
    // Find matrix of scores M_i,j with individual base pairs (i, j)
    // forced into alignment
    const n = seq1.length;
    const m = seq2.length;

    // This is the best score so far involving a match at position i, j
    // which is the global alignment score for the left portion of
    // the strings up to and including positions i, j..
    const left = glob_align_scores(seq1, seq2);

    // This is the best score so far involving a match at position i, j
    // for the reverse strings
    // which is the global alignment score for the right portion of
    // the strings up to and including positions i, j..
    const right = glob_align_scores(rev(seq1), rev(seq2));

    // The score is the sum of the two (reversed for the right portion)
    // minus the score of the actual pair of bases themselves, which are
    // double counted
    var score : [1..n, 1..m] int;
    for i in 1..n do
      for j in 1..m {
        score[i, j] = left[i, j, move.diag];
        score[i, j] -= matchscore(seq1, i, seq2, j);
        score[i, j] += right[n-i+1, m-j+1, move.diag];
      }

    return score;
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());
      
      var score = isolation_matrix(seqs[2], seqs[1]);

      writeln(max reduce score);
      writeln(+ reduce score);
    }
  }
}
