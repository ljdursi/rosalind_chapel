module OAP {
  config const defaultfilename = "none";
  config const infile = defaultfilename;
  config const gap_score = -2;
  config const match_score = +1;
  config const mismatch_score = -2;

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

  proc overlap_alignment(seq1, seq2, match_score, mismatch_score, gap_score) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var moves : [0..n, 0..m] move;

    moves = move.start;
    for i in 1..n {
      score[i, 0] = i*gap_score;
      moves[i, 0] = move.ins;
    }

    for i in 1..n {
      const basei = ascii(seq1(i));
      for j in 1..m {
        var match_delta = mismatch_score;
        if basei == ascii(seq2(j)) then
          match_delta = match_score;

        const match = score[i-1, j-1] + match_delta;
        const ins = score[i-1, j] + gap_score;
        const del = score[i, j-1] + gap_score;

        (score[i, j], moves[i,j]) = max((ins, move.ins), (del, move.del),
                                        (match, move.diag));
      }
    }

    // backtrack
    var j = m;
    const dist = max reduce score[0..n, j];

    // highest i that has maximum score
    var i = 0;
    for maxi in 0..n do 
      if score[maxi, j] == dist then
        i = maxi;

    var aug1, aug2: string;

    while true {
      if i <= 0 || j <= 0 then
        break;

      if moves[i, j] == move.start then
        break;

      if moves[i, j] == move.diag {
        aug1 = seq1[i] + aug1;
        aug2 = seq2[j] + aug2;
        i -= 1;
        j -= 1;
      } else if moves[i, j] == move.del {
        aug2 = seq2[j] + aug2;
        aug1 = "-" + aug1;
        j -= 1;
      } else {
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

      const score = overlap_alignment(seqs[2], seqs[1], match_score, mismatch_score, gap_score);
      writeln(score[1]);
      writeln(score[3]);
      writeln(score[2]);
    }
  }
}
