module MPRT {
  use Regexp;
  use Curl;
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const motif="N[^P][ST][^P]";

  iter fasta_iterator(inchannel) throws {
    var sequence = "", seqlabel = "";
    var started = false;
    var line = "";

    while (inchannel.readline(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..);
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  // read file into an array of strings
  // use stdin if filename == defaultfilename
  proc strings_from_file(filename: string, defaultfilename: string) throws {
    var text: string = "";
    var line: string;
    var lines: [1..0] string;

    var channel = stdin;
    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    while (channel.read(line)) {
      lines.push_back(line);
    }

    return lines;
  }

  proc sequence_from_uniprot(uniprotid) throws {
    const uniproturl="http://www.uniprot.org/uniprot/" + uniprotid + ".fasta";
    var c = open(url=uniproturl, mode=iomode.r);
    c.setopt(curlopt_followlocation, 1);
    var reader = c.reader();

    var sequences: [1..0] string;
    for (seqlabel, sequence) in fasta_iterator(reader) do
      sequences.push_back(sequence);

    return sequences[1];
  }

  proc motif_locations(sequence, motif_re) {
    var posns : [1..0] int;
    var startidx = 1;

    // Note: cant use re.matches because that skips
    // overlapping matches
    var match = motif_re.search(sequence(startidx..));
    while match {
      const idx = startidx + match.offset;
      posns.push_back(idx);
      startidx = idx + 1;
      match = motif_re.search(sequence(startidx..));
    }
    return posns;
  }

  proc main() {
    try! {
        var uniprot_ids = strings_from_file(infile, defaultfilename);
        for uniprot_id in uniprot_ids {
          const seq = sequence_from_uniprot(uniprot_id);
          var motif_re = compile(motif);
          const locs = motif_locations(seq, motif_re);
          if locs.size != 0 {
            writeln(uniprot_id);
            for loc in locs do 
              write(loc, " ");
            writeln("");
          }
        }
    }
  }
}
