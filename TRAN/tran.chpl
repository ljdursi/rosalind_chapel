module TRAN {
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

  proc transition_transversion(seq1, seq2) {
    const basetype = [ "A" => "purine", "G" => "purine", 
                       "C" => "pyramidine", "T" => "pyramidine" ];
    var transitions = 0;
    var transversions = 0;

    for (base1, base2) in zip(seq1, seq2) {
        if base1 != base2 {
            if basetype[base1] != basetype[base2] then
                transversions += 1;
            else
                transitions += 1;
        }
    }
    return (transitions, transversions);
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        var sequences : [1..0] string;
        for (seqlabel, sequence) in fasta_iterator(inchannel) do
          sequences.push_back(sequence);

        var tran = transition_transversion(sequences[1], sequences[2]);
        writeln(1.0*tran[1]/tran[2]);
    }
  }
}
