module LCIS {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var n: int;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.read(n);
    var sequence : [1..n] int;
    success = channel.read(sequence);
    return sequence;
  }

  proc lcs(seq1, seq2) {
    const n = seq1.size;
    const m = seq2.size;
    var length : [0..n, 0..m] int;
    var i, j : int;

    // calculate subsequence lengths
    length = 0;
    for i in 1..n {
      for j in 1..m {
        if seq1[i] == seq2[j] then
          length[i, j] = length[i-1, j-1] + 1;
        else 
          length[i, j] = max(length[i-1, j], length[i, j-1]);
      }
    }

    // backtrack
    const l = length[n, m];
    var subsequence : [1..l] int;
    var position = l;
    i = n;
    j = m;

    while (i > 0 && j > 0) {
        var up = length[i-1, j];
        var left = length[i, j-1];
        if up < position && left < position {
          subsequence[position] = seq1[i];
          position -= 1;
          i -= 1;
          j -= 1;
        } else {
          if left > up then
            j -= 1;
          else
            i -= 1;
        }
    }
    return subsequence;
  }

  proc main() {
    try! {
      var sequence = params_from_file(infile, defaultfilename);
      var sortedseq : [1..sequence.size] int;
      var i : int;

      for i in 1..sequence.size do
        sortedseq[i] = i;
      const result = lcs(sequence, sortedseq);
      writeln(result);

      for i in 1..sequence.size do
        sortedseq[i] = sequence.size - i + 1;
      const result2 = lcs(sequence, sortedseq);
      writeln(result2);
    }
  }
}
