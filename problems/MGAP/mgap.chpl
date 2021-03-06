module MGAP {
  // Maximum gap alignment
  // http://rosalind.info/problems/mgap
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const match=+10;
  config const mismatch=-3;
  config const gap=-1;

  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
    // simple fasta parser
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

  proc edit_dist(seq1, seq2, match_score, mismatch_score, gap_score) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var ngaps : [0..n, 0..m] int;
    var i, j : int;

    score[0, 1..m] = (1..m)*gap_score;
    score[1..n, 0] = (1..n)*gap_score;
    ngaps[0, 1..m] = 1..m;
    ngaps[1..n, 0] = 1..n;

    // calculate subsequence lengths
    for i in 1..n {
      for j in 1..m {
        const diag = score[i-1, j-1] + if ascii(seq1[i]) == ascii(seq2[j]) then match_score else mismatch_score;
        const insertion = score[i-1, j] + gap_score;
        const deletion = score[i, j-1] + gap_score;
        (score[i, j], ngaps[i, j]) = max((insertion, ngaps[i-1, j] + 1),
                                         (deletion, ngaps[i, j-1] + 1),
                                         (diag, ngaps[i-1, j-1]));
      }
    }

    return (score[n, m], ngaps[n, m]);
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
      var sequences : [1..0] string;
      for (seqlabel, sequence) in fasta_iterator(inchannel) {
        sequences.push_back(sequence);
      }
      var (score, gaps) = edit_dist(sequences[1], sequences[2], match, mismatch, gap);
      writeln(gaps);
    }
  }
}
