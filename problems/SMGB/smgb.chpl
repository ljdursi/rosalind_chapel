module SMGB {
  config const defaultfilename = "none";
  config const infile = defaultfilename;
  config const gap_score = -1;
  config const match_score = +1;
  config const mismatch_score = -1;

  enum move {ins=1, del=2, diag=4};

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

  proc semiglobal_alignment(seq1, seq2, match_score, mismatch_score, gap_score) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var moves : [0..n, 0..m] move;
    moves[0..n, 0] = move.ins;
    moves[0, 0..m] = move.del;

    score = 0;

    for i in 1..n {
      const del_gap = if i < n then gap_score else 0;

      for j in 1..m {
        const match_delta = if seq1(i) == seq2(j) then match_score else mismatch_score;
        const ins_gap = if j < m then gap_score else 0;

        const match = score[i-1, j-1] + match_delta;
        const ins = score[i-1, j] + ins_gap;
        const del = score[i, j-1] + del_gap;

        (score[i, j], moves[i,j]) = max((ins, move.ins), (del, move.del),
                                        (match, move.diag));
      }
    }

    // backtrack
    var j = m;
    var i = n;
    const dist = score[n, m];

    var aug1, aug2: string;

    while i > 0 || j > 0 {
      if moves[i, j] == move.diag {
        aug1 = seq1[i] + aug1;
        aug2 = seq2[j] + aug2;
        i -= 1;
        j -= 1;
      } else if moves[i, j] == move.del {
        aug2 = seq2[j] + aug2;
        aug1 = "-" + aug1;
        j -= 1;
      } else if moves[i, j] == move.ins {
        aug1 = seq1[i] + aug1;
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

      const score = semiglobal_alignment(seqs[2], seqs[1], match_score, mismatch_score, gap_score);
      writeln(score[1]);
      writeln(score[3]);
      writeln(score[2]);
    }
  }
}
