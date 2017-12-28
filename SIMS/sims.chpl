module SIMS {
  const defaultfilename="none";
  config const infile=defaultfilename;
  enum move {ins, del, diag, start};

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

  proc motif_alignment(motif, seq2) {
    const n = motif.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var moves : [0..n, 0..m] move;

    score[1..n, 0] = -n..-1 by -1;
    moves = move.start;

    for i in 1..n {
      for j in 1..m {
        const match = score[i-1, j-1] + 1;
        const mismatch = score[i-1, j-1] - 1;
        const ins = score[i-1, j] - 1;
        const del = score[i, j-1] - 1;

        if motif[i] == seq2[j] then
          (score[i, j], moves[i,j]) = max((ins, move.ins), (del, move.del), (match, move.diag));
        else 
          (score[i, j], moves[i,j]) = max((ins, move.ins), (del, move.del), (mismatch, move.diag));
      }
    }

    // backtrack
    var (dist, j) = maxloc reduce zip(score[n, 0..m], 0..m);
    var i = n;
    var aug1, aug2: string;

    while i > 0 && j > 0 {
      if moves[i, j] == move.diag {
        aug1 = motif[i] + aug1;
        aug2 = seq2[j] + aug2;
        i -= 1;
        j -= 1;
      } else if moves[i, j] == move.del {
        aug1 = "-" + aug1;
        aug2 = seq2[j] + aug2;
        j -= 1;
      } else {
        aug1 = motif[i] + aug1;
        aug2 = "-" + aug2;
        i -= 1;
      }
    }

    return (dist, aug1, aug2);
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var (n, aug1, aug2) = motif_alignment(seqs[2], seqs[1]);
      writeln(n);
      writeln(aug2);
      writeln(aug1);
    }
  }
}
