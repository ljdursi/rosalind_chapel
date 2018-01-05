module CORR {
  const defaultfilename="none";
  config const infile=defaultfilename;
  use Sort;

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


  proc reverse_complement(seq) {
    const matchings = ['A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A'];
    var rc = "";
    
    for base in seq[1..seq.length by -1] do
      rc += matchings[base];

    return rc;
  }

  proc normalize(seq) {
    const rc = reverse_complement(seq);
    if rc < seq then
      return rc;
    else
      return seq;
  }

  proc known_good_seqs(seqs) {
    var lastseq = "";
    var known_good : domain(string);
    var normedseqs : [1..seqs.size] string;
    for i in 1..seqs.size do
      normedseqs[i] = normalize(seqs[i]);

    for seq in sorted(normedseqs) {
      if lastseq != "" then
        if seq == lastseq {
          known_good += seq;
          known_good += reverse_complement(seq);
        }
      lastseq = seq;
    }
    return known_good;
  }

  proc find_neighbor(seq, known_good) {
    for good in known_good {
      var edit_dist = 0;
      for (base1, base2) in zip(seq, good) {
        if base1 != base2 then
          edit_dist += 1;
        if edit_dist > 1 then
          break;
      }
      if edit_dist == 1 then
        return good;
    }
    return "";
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
      var sequences : [1..0] string;
      for (seqlabel, sequence) in fasta_iterator(inchannel) do
        sequences.push_back(sequence);

      var known_good = known_good_seqs(sequences);
      for seq in sequences do
        if !known_good.member(seq) {
          var fixed = find_neighbor(seq, known_good);
          writeln(seq, "->", fixed);
        }
    }
  }
}
