module GC {
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
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

    if (started) then
      yield (seqlabel, sequence);
  }

  proc sequence_gc(sequence: string) {
    const n = sequence.len;
    var gc = 0.0;
    for base in sequence do
      if (base == "G" || base == "C") then
        gc += 1.0;

    return 100.0*gc/n;
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var maxlabel = "none", maxgc = -1.0;

      for (seqlabel, sequence) in fasta_iterator(inchannel) {
        var gc = sequence_gc(sequence);
        if (gc > maxgc) then
           (maxlabel, maxgc) = (seqlabel, gc);
      }
      writeln(maxlabel);
      writeln(maxgc);
    }
  }
}
