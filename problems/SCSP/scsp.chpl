module SCSP {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc seqs_from_file(filename: string, defaultfilename: string) throws {
    var seq1, seq2: string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var success = channel.readln(seq1);
    success = channel.read(seq2);
    return (seq1.strip(), seq2.strip());
  }

  proc scs(seq1, seq2) {
    const n = seq1.length;
    const m = seq2.length;
    var shortest : [-1..n, -1..m] int;
    var i, j : int;

    shortest = 0;
    shortest[0, 1..m] = 1..m;
    shortest[1..n ,0] = 1..n;

    for i in 1..n do
      for j in 1..m do
        if seq1[i] == seq2[j] then
          shortest[i, j] = shortest[i-1, j-1] + 1;
        else 
          shortest[i, j] = min(shortest[i-1, j] + 1, shortest[i, j-1] + 1);

    // backtrack
    const l = shortest[n, m];
    var supersequence = "";
    var position = l;
    i = n;
    j = m;

    while (i > 0 || j > 0) {
      var up = shortest[i-1, j];
      var left = shortest[i, j-1];
      var diag = shortest[i-1, j-1];
      if (diag + 1 == position) && (i >0 && j > 0) && (seq1[i] == seq2[j]) {
        supersequence = seq1[i] + supersequence;
        i -= 1;
        j -= 1;
      } else if (left <= up) && (j > 0) {
        supersequence = seq2[j] + supersequence;
        j -= 1;
      } else {
        supersequence = seq1[i] + supersequence;
        i -= 1;
      }
      position -= 1;
    }
    return supersequence;
  }

  proc main() {
    try! {
      var (seq1, seq2) = seqs_from_file(infile, defaultfilename);
      writeln(scs(seq1, seq2));
    }
  }
}
