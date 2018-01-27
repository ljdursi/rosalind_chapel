module EDTA {
  const defaultfilename="none";
  config const infile=defaultfilename;
  enum move {ins, del, diag};

  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
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

  proc edit_alignment(seq1, seq2) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var moves : [0..n, 0..m] move;

    score[0, 1..m] = 1..m;
    score[1..n, 0] = 1..n;
    moves = move.ins;

    var i, j: int;
    for i in 1..n {
      for j in 1..m {
        const match = score[i-1, j-1];
        const mismatch = score[i-1, j-1] + 1;
        const ins = score[i-1, j] + 1;
        const del = score[i, j-1] + 1;

        if ascii(seq1[i]) == ascii(seq2[j]) then
          (score[i, j], moves[i,j]) = min((ins, move.ins), (del, move.del), (match, move.diag));
        else 
          (score[i, j], moves[i,j]) = min((ins, move.ins), (del, move.del), (mismatch, move.diag));
      }
    }

    // backtrack
    const editdist = score[n, m];
    var aug1, aug2: string;
    (i, j) = (n, m);

    while i > 0 || j > 0 {
      if moves[i, j] == move.diag {
        aug1 = seq1[i] + aug1;
        aug2 = seq2[j] + aug2;
        i -= 1;
        j -= 1;
      } else if moves[i, j] == move.del {
        aug1 = "-" + aug1;
        aug2 = seq2[j] + aug2;
        j -= 1;
      } else {
        aug1 = seq1[i] + aug1;
        aug2 = "-" + aug2;
        i -= 1;
      }
    }

    return (editdist, aug1, aug2);
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var (n, aug1, aug2) = edit_alignment(seqs[1], seqs[2]);
      writeln(n);
      writeln(aug1);
      writeln(aug2);
    }
  }
}
