module LCSM {
  use Sort;
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

    while (inchannel.read(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..);
        sequence = "";
      } else {
        sequence += line;
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }


  iter substr_iterator(sequence) {
    const len = sequence.length;
    for sublen in 1..len by -1 {
      for start in 1..(len-sublen+1) do
        yield sequence(start..#sublen);
    }
  }

  proc naive_lcsm(sequences): string {
    var sorted_seqs = sequences;
    sort(sorted_seqs);
    const shortest = sorted_seqs[1];
    const others = sorted_seqs[2..];

    for substr in substr_iterator(shortest) {
       var allpresent = true;
       for otherstr in others {
         if otherstr.find(substr) == 0 {
            allpresent = false;
            break;
         }
       }
       if allpresent then
         return substr;
    }
    return "";
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    var sequences : [1..0] string;
    try! {
        for (seqlabel, sequence) in fasta_iterator(inchannel) do
          sequences.push_back(sequence);

        const lcsm = naive_lcsm(sequences);
        writeln(lcsm);
    }
  }
}
