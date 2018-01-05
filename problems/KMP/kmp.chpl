module KMP {
  param defaultfilename="none";
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

  proc failure_array(sequence) {
    const n = sequence.length;
    var fa : [1..#n] int;

    var prefixlen = 0;
    var prefixchar = sequence(prefixlen+1);

    for (i, base) in zip(2..n, sequence(2..n)) {
      while true {
        if base == prefixchar {
          prefixlen += 1;
          prefixchar = sequence(prefixlen+1);
          break;
        } 
        if prefixlen == 0 then
          break;
        prefixlen = fa[prefixlen];
        prefixchar = sequence(prefixlen+1);
      }
      fa[i] = prefixlen;
    }

    return fa;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        for (seqlabel, sequence) in fasta_iterator(inchannel) do
          writeln(failure_array(sequence));
    }
  }
}
