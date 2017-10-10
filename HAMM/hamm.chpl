module HAMM {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
  proc params_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    var line1="", line2="";

    var success = channel.read(line1);
    success = channel.read(line2);
    return (line1, line2);
  }

  proc hamming_distance(seq1: string, seq2: string): int {
    var dist = 0;
    for (base1, base2) in zip(seq1, seq2) do
      if base1 != base2 then
        dist += 1;

    return dist;
  }

  proc main() {
    var seq1, seq2: string;
    try! {
      (seq1, seq2) = params_from_file(infile, defaultfilename);
      var result = hamming_distance(seq1, seq2);
      writeln(result);
    }
  }
}
