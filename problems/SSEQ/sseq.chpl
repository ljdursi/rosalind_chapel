module SSEQ {
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

  proc subsequence_indices(sequence, subsequence) {
    var seq_cursor = 1;
    var subseq_positions : [1..subsequence.length] int;

    for i in 1..subsequence.length {
        while sequence[seq_cursor] != subsequence[i] do
            seq_cursor += 1;
        subseq_positions[i] = seq_cursor;
        seq_cursor += 1;
    }

    return subseq_positions;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        var sequences : [1..0] string;
        for (seqlabel, sequence) in fasta_iterator(inchannel) do
          sequences.push_back(sequence);

        var indices = subsequence_indices(sequences[1], sequences[2]);
        writeln(indices);
    }
  }
}
